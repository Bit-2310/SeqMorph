# src/outputs/manager.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from typing import Any, Iterable, Optional
from contextlib import contextmanager
import os
import json
import shutil


# Root directory for all outputs; can be overridden via environment
DEFAULT_ROOT = Path(os.environ.get("SEQMORPH_OUTPUTS_DIR", "outputs")).resolve()


@dataclass(frozen=True, slots=True)
class RunPaths:
    """
    Structured paths for a single run. Use these to write artifacts/reports/logs.
    """
    root: Path
    artifacts: Path
    reports: Path
    logs: Path
    manifest: Path

    def exists(self) -> bool:
        return self.root.exists()


class OutputsManager:
    """
    Small helper to standardize where SeqMorph writes outputs.

    Layout:
        <root>/runs/<run_id>/
            artifacts/
            reports/
            logs/
            manifest.json
        <root>/index.json

    Features:
        - Per-run folder creation with timestamped run_id.
        - Global index append with each new run.
        - Simple list/get/delete helpers for GUI/API use.
        - Optional retention: keep only the most-recent N runs (FIFO by run_id timestamp).
        - Safe file resolution within a run to avoid directory traversal.

    Notes:
        - JSON files are written with UTF-8 and indent=2.
        - All methods are synchronous and filesystem-only.
    """

    def __init__(
        self,
        root: Optional[Path | str] = None,
        *,
        retention: Optional[int] = None,
    ) -> None:
        """
        Parameters
        ----------
        root : Path | str | None
            Outputs root. Defaults to DEFAULT_ROOT (env SEQMORPH_OUTPUTS_DIR or ./outputs).
        retention : int | None
            If set, keep only the most recent N runs. Older ones are removed after new_run().
        """
        self.root: Path = Path(root).resolve() if root is not None else DEFAULT_ROOT
        self.runs: Path = self.root / "runs"
        self.index: Path = self.root / "index.json"
        self.retention: Optional[int] = retention

        self.runs.mkdir(parents=True, exist_ok=True)

        if not self.index.exists():
            # initialize empty index for convenience
            self._write_json(self.index, [])

    # --------------------------------------------------------------------- #
    # Creation
    # --------------------------------------------------------------------- #

    def new_run(self, label: Optional[str] = None, meta: Optional[dict[str, Any]] = None) -> RunPaths:
        """
        Create a new run directory tree and manifest; append to global index.

        Parameters
        ----------
        label : str | None
            Optional label to suffix the run_id (e.g., "mutate_and_analyze").
        meta : dict | None
            Extra metadata to include in the manifest (merged under "meta").

        Returns
        -------
        RunPaths
            Paths object for the newly created run.
        """
        ts = self._timestamp()
        run_id = f"{ts}{f'_{label}' if label else ''}"

        base = self.runs / run_id
        artifacts = base / "artifacts"
        reports = base / "reports"
        logs = base / "logs"

        for p in (artifacts, reports, logs):
            p.mkdir(parents=True, exist_ok=True)

        manifest = base / "manifest.json"
        manifest_obj: dict[str, Any] = {
            "run_id": run_id,
            "created": ts,
            "label": label,
        }
        if meta:
            manifest_obj["meta"] = meta

        self._write_json(manifest, manifest_obj)

        # Append to index
        idx_entry = {"run_id": run_id, "created": ts, "label": label}
        self._append_index(idx_entry)

        # Enforce retention policy
        if self.retention is not None and self.retention > 0:
            self._prune_retention(self.retention)

        return RunPaths(base, artifacts, reports, logs, manifest)

    # --------------------------------------------------------------------- #
    # Listing / lookup / deletion
    # --------------------------------------------------------------------- #

    def list_runs(self) -> list[dict[str, Any]]:
        """
        Return the global index as a list (newest last). Empty list if none.
        """
        try:
            return self._read_json(self.index)  # type: ignore[return-value]
        except Exception:
            return []

    def get_run_paths(self, run_id: str) -> RunPaths:
        """
        Build RunPaths for an existing run_id (does not validate content).
        Raises FileNotFoundError if the run directory is missing.
        """
        base = self.runs / run_id
        if not base.exists():
            raise FileNotFoundError(f"Run not found: {run_id}")
        return RunPaths(
            root=base,
            artifacts=base / "artifacts",
            reports=base / "reports",
            logs=base / "logs",
            manifest=base / "manifest.json",
        )

    def read_manifest(self, run_id: str) -> dict[str, Any]:
        """
        Load manifest.json for a run.
        """
        rp = self.get_run_paths(run_id)
        if not rp.manifest.exists():
            raise FileNotFoundError(f"Manifest missing for run: {run_id}")
        return self._read_json(rp.manifest)

    def update_manifest(self, run_id: str, update: dict[str, Any]) -> None:
        """
        Merge `update` into existing manifest (shallow merge).
        """
        rp = self.get_run_paths(run_id)
        obj = {}
        if rp.manifest.exists():
            obj = self._read_json(rp.manifest)
        obj.update(update)
        self._write_json(rp.manifest, obj)

    def delete_run(self, run_id: str) -> bool:
        """
        Delete a run directory and remove it from index.
        Returns True on success, False if the run was missing.
        """
        base = self.runs / run_id
        if not base.exists():
            return False

        shutil.rmtree(base, ignore_errors=True)

        # Remove from index (if present)
        try:
            entries: list[dict[str, Any]] = self._read_json(self.index)
            entries = [e for e in entries if e.get("run_id") != run_id]
            self._write_json(self.index, entries)
        except Exception:
            # If index is unreadable, re-scan runs dir to rebuild
            self._rebuild_index()

        return True

    # --------------------------------------------------------------------- #
    # Safe file IO within a run
    # --------------------------------------------------------------------- #

    def safe_path(self, run_id: str, rel_path: str | Path) -> Path:
        """
        Resolve a path inside the run directory, preventing directory traversal.
        """
        base = self.get_run_paths(run_id).root.resolve()
        target = (base / Path(rel_path)).resolve()
        if not str(target).startswith(str(base)):
            raise ValueError("Attempted path escape outside the run directory.")
        return target

    def write_text(self, run_id: str, rel_path: str | Path, content: str, *, encoding: str = "utf-8") -> Path:
        """
        Write a text file inside the run (creates parent dirs). Returns the written path.
        """
        p = self.safe_path(run_id, rel_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding=encoding)
        return p

    def read_text(self, run_id: str, rel_path: str | Path, *, encoding: str = "utf-8") -> str:
        """
        Read a text file from inside the run.
        """
        p = self.safe_path(run_id, rel_path)
        return p.read_text(encoding=encoding)

    def write_json(self, run_id: str, rel_path: str | Path, obj: Any) -> Path:
        """
        Write a JSON file inside the run with indent=2.
        """
        p = self.safe_path(run_id, rel_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self._write_json(p, obj)
        return p

    def read_json(self, run_id: str, rel_path: str | Path) -> Any:
        """
        Read JSON from a file inside the run.
        """
        p = self.safe_path(run_id, rel_path)
        return self._read_json(p)

    # --------------------------------------------------------------------- #
    # Context manager (ergonomic "session")
    # --------------------------------------------------------------------- #

    @contextmanager
    def session(self, label: Optional[str] = None, meta: Optional[dict[str, Any]] = None) -> Iterable[RunPaths]:
        """
        Context manager to create a run and optionally finalize the manifest.

        Example:
            with om.session(label="mutate_and_analyze", meta={"seed": 42}) as run:
                (run.logs / "step.log").write_text("started")
                ...
        """
        run = self.new_run(label=label, meta=meta)
        try:
            yield run
        finally:
            # Add a simple "completed" timestamp; callers can update_manifest() with richer info.
            if run.exists():
                manifest = self._safe_read_json(run.manifest) or {}
                manifest["completed"] = self._timestamp()
                self._write_json(run.manifest, manifest)

    # --------------------------------------------------------------------- #
    # Internal helpers
    # --------------------------------------------------------------------- #

    def _append_index(self, entry: dict[str, Any]) -> None:
        try:
            entries: list[dict[str, Any]] = self._read_json(self.index)
        except Exception:
            entries = []
        entries.append(entry)
        self._write_json(self.index, entries)

    def _rebuild_index(self) -> None:
        """
        Rebuild index.json from the runs directory if it's missing or corrupt.
        """
        items: list[dict[str, Any]] = []
        if self.runs.exists():
            for child in sorted(self.runs.iterdir()):
                if not child.is_dir():
                    continue
                manifest = child / "manifest.json"
                if manifest.exists():
                    obj = self._safe_read_json(manifest) or {}
                    items.append({
                        "run_id": obj.get("run_id", child.name),
                        "created": obj.get("created"),
                        "label": obj.get("label"),
                    })
                else:
                    items.append({
                        "run_id": child.name,
                        "created": None,
                        "label": None,
                    })
        self._write_json(self.index, items)

    def _prune_retention(self, keep_last: int) -> None:
        """
        Keep only the most recent `keep_last` runs (by run_id lexical order: timestamp prefix).
        """
        # Determine current runs by directory names
        names = sorted([p.name for p in self.runs.iterdir() if p.is_dir()])
        if keep_last <= 0 or len(names) <= keep_last:
            return
        to_delete = names[: len(names) - keep_last]
        for name in to_delete:
            shutil.rmtree(self.runs / name, ignore_errors=True)

        # Rebuild index to reflect deletions
        self._rebuild_index()

    @staticmethod
    def _timestamp() -> str:
        # ISO-like but filesystem-friendly: YYYYMMDD_HHMMSS
        return datetime.now().strftime("%Y%m%d_%H%M%S")

    @staticmethod
    def _write_json(path: Path, obj: Any) -> None:
        path.write_text(json.dumps(obj, indent=2), encoding="utf-8")

    @staticmethod
    def _read_json(path: Path) -> Any:
        return json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def _safe_read_json(path: Path) -> Optional[dict[str, Any]]:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
