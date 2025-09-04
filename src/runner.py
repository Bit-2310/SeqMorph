# src/runner.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Callable, Any
from time import perf_counter
import json

from outputs.manager import OutputsManager, RunPaths
from mutation.engine import MutationEngine
from mutation.context import ContextChooserConfig
from mutation.types import ChooseBase
from sequence_store.string_store import StringStore
from sequence_store.base import BaseStore

@dataclass(slots=True)
class RunSummary:
    """Concise summary of a single mutate+analyze run."""
    run_id: str
    sequence_id: str
    length: int
    gc: float
    point_mutations: int
    structural_events: int
    status: str = "success"
    message: str = ""


class Runner:
    """
    Orchestrates mutate+analyze flows and writes standardized outputs via OutputsManager.
    Designed for reuse by CLI, API/GUI, or notebooks.

    Typical usage:
        r = Runner(retention=50)
        r.add_sequence("s1", "ACGCGTACGT")
        summary = r.mutate_and_analyze("s1", point_positions=[2,5,8], seed=123)
    """

    def __init__(
        self,
        *,
        outputs_root: Optional[Path | str] = None,
        retention: Optional[int] = None,
        progress: Optional[Callable[[str, dict[str, Any]], None]] = None,
    ) -> None:
        """
        Parameters
        ----------
        outputs_root : Path | str | None
            Root directory for outputs (defaults to SEQMORPH_OUTPUTS_DIR or ./outputs).
        retention : int | None
            If set, keep only the most recent N runs (FIFO by timestamped run_id).
        progress : callable | None
            Optional callback receiving (event_name, context_dict) at key milestones.
        """
        self.om = OutputsManager(root=outputs_root, retention=retention)
        self._sequences: dict[str, str] = {}
        self._progress = progress

    # ------------------------------------------------------------------ #
    # Sequence registry (minimal, in-memory for now)
    # ------------------------------------------------------------------ #

    def add_sequence(self, sid: str, sequence: str) -> None:
        """
        Register a sequence by id. RNA U is mapped to T. Allowed alphabet: A,C,G,T,N.
        """
        up = sequence.upper().replace("U", "T")
        if not set(up).issubset({"A", "C", "G", "T", "N"}):
            raise ValueError("Invalid DNA: allowed symbols are A,C,G,T,N (U is mapped to T).")
        self._sequences[sid] = up

    def get_sequence(self, sid: str) -> str:
        """
        Retrieve a registered sequence by id.
        """
        if sid not in self._sequences:
            raise ValueError(f"Unknown sequence_id: {sid}")
        return self._sequences[sid]

    # ------------------------------------------------------------------ #
    # Core flow: mutate + analyze (single job)
    # ------------------------------------------------------------------ #

    def mutate_and_analyze(
        self,
        sequence_id: str,
        *,
        point_positions: Optional[List[int]] = None,
        struct_window: Optional[tuple[int, int]] = None,
        rate_pct: float = 0.0,
        mean_seg_len: int = 500,
        mix: Optional[Dict[str, float]] = None,
        seed: Optional[int] = 42,
        chooser_cfg: Optional[ContextChooserConfig] = None,
        choose_base: Optional[ChooseBase] = None,
        label: str = "mutate_and_analyze",
        model_path: Optional[str] = None,
        pre_hooks: Optional[List[Callable[[RunPaths, dict[str, Any]], None]]] = None,
        post_hooks: Optional[List[Callable[[RunPaths, dict[str, Any]], None]]] = None,
    ) -> RunSummary:
        """
        Perform point and/or structural mutations, compute simple metrics, and
        write a complete output bundle (artifacts/reports/logs + manifest).

        Returns
        -------
        RunSummary
            A small, serializable result for CLI/API printing or return payloads.
        """
        seq = self.get_sequence(sequence_id)
        self._validate_params(len(seq), point_positions, struct_window, rate_pct)

        # Create run and record provenance
        meta = {
            "seed": seed,
            "model_path": model_path,
            "point_positions": point_positions,
            "struct_window": struct_window,
            "rate_pct": rate_pct,
            "mean_seg_len": mean_seg_len,
            "mix": mix,
        }
        run = self.om.new_run(label=label, meta=meta)
        self._emit("run_created", {"run_id": run.root.name, "sequence_id": sequence_id})

        # Write input artifacts & initial log
        self._write_input_artifacts(run, sequence_id, seq)

        # Pre-hooks for custom behaviors (e.g., extra input validation or metrics)
        ctx: dict[str, Any] = {"sequence_id": sequence_id}
        for h in pre_hooks or []:
            h(run, ctx)

        start = perf_counter()
        status = "success"
        message = ""

        try:
            # Engine + store
            store = StringStore(seq)
            engine = MutationEngine(
                seed=seed,
                chooser_cfg=chooser_cfg or ContextChooserConfig(),
                default_choose_base=choose_base,
                model_path=model_path,
            )
            # Record engine/chooser provenance
            self.om.update_manifest(run.root.name, {
                "engine": "MutationEngine",
                "chooser_cfg": {
                    "cpg_enabled": (chooser_cfg.cpg_enabled if chooser_cfg else True),
                    "cpg_strength": (chooser_cfg.cpg_strength if chooser_cfg else 3.0),
                },
                "store": "StringStore",
            })

            # Mutations
            pm = 0
            if point_positions:
                pm = engine.mutate_points(store, point_positions, choose_base=choose_base)

            structural_events = []
            if struct_window and rate_pct > 0.0:
                s, e = struct_window
                structural_events = engine.mutate_structural(
                    store, win_start=s, win_end=e, rate_pct=rate_pct,
                    mean_seg_len=mean_seg_len, mix=mix
                )

            mutated = store.to_string()

            # Simple analysis (extendable)
            gc = (mutated.count("G") + mutated.count("C")) / max(1, len(mutated))

            # Write outputs
            self._write_outputs(run, sequence_id, mutated, gc, pm, len(structural_events))

            # Final manifest enrichment
            self.om.update_manifest(run.root.name, {
                "inputs": {"sequence_id": sequence_id},
                "result": {
                    "length": len(mutated),
                    "gc": gc,
                    "point_mutations": pm,
                    "structural_events": len(structural_events),
                },
            })

            # Post-hooks (e.g., extra analyses into reports/)
            for h in post_hooks or []:
                h(run, {"gc": gc, "length": len(mutated)})

            self._emit("run_completed", {"run_id": run.root.name})

            elapsed_ms = int((perf_counter() - start) * 1000)
            self.om.update_manifest(run.root.name, {"elapsed_ms": elapsed_ms, "status": status})

            return RunSummary(
                run_id=run.root.name,
                sequence_id=sequence_id,
                length=len(mutated),
                gc=gc,
                point_mutations=pm,
                structural_events=len(structural_events),
                status=status,
                message=message,
            )

        except Exception as exc:
            status = "failed"
            message = str(exc)
            elapsed_ms = int((perf_counter() - start) * 1000)
            self.om.update_manifest(run.root.name, {"elapsed_ms": elapsed_ms, "status": status, "error": message})
            self._emit("run_failed", {"run_id": run.root.name, "error": message})
            # Still return a summary (useful to show in UI)
            return RunSummary(
                run_id=run.root.name,
                sequence_id=sequence_id,
                length=0,
                gc=0.0,
                point_mutations=0,
                structural_events=0,
                status=status,
                message=message,
            )

    # ------------------------------------------------------------------ #
    # Batch mode: run many jobs (each gets its own run folder)
    # ------------------------------------------------------------------ #

    def batch_mutate_and_analyze(self, jobs: List[dict[str, Any]]) -> List[RunSummary]:
        """
        Run multiple mutate+analyze jobs. Each item is a dict that must include
        'sequence_id' and may include any kwargs accepted by mutate_and_analyze.
        Returns a list of RunSummary (successful or failed).
        """
        results: List[RunSummary] = []
        for job in jobs:
            try:
                sid = job.get("sequence_id")
                if not sid:
                    raise ValueError("Job missing 'sequence_id'.")
                # Shallow copy to avoid mutating caller-provided dict
                params = {k: v for k, v in job.items() if k != "sequence_id"}
                res = self.mutate_and_analyze(sid, **params)
                results.append(res)
            except Exception as exc:
                results.append(RunSummary(
                    run_id="ERROR",
                    sequence_id=str(job.get("sequence_id", "?")),
                    length=0,
                    gc=0.0,
                    point_mutations=0,
                    structural_events=0,
                    status="failed",
                    message=str(exc),
                ))
        return results

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    def _validate_params(
        self,
        seq_len: int,
        point_positions: Optional[List[int]],
        struct_window: Optional[tuple[int, int]],
        rate_pct: float,
    ) -> None:
        """
        Validate inputs; raise ValueError on invalid parameters.
        """
        if point_positions:
            if any((p is None) or (p < 0) or (p >= seq_len) for p in point_positions):
                raise ValueError("point_positions contain out-of-bounds indices.")
            if len(point_positions) != len(set(point_positions)):
                raise ValueError("point_positions must be unique.")
        if struct_window:
            s, e = struct_window
            if not (0 <= s < e <= seq_len):
                raise ValueError("struct_window must be within sequence bounds and s < e.")
            if not (0.0 <= rate_pct <= 100.0):
                raise ValueError("rate_pct must be between 0 and 100.")

    def _write_input_artifacts(self, run: RunPaths, sid: str, seq: str) -> None:
        """
        Write input FASTA and start log.
        """
        (run.artifacts / "input.fasta").write_text(f">{sid}\n{seq}\n", encoding="utf-8")
        (run.logs / "runner.log").write_text("Run started\n", encoding="utf-8")

    def _write_outputs(
        self,
        run: RunPaths,
        sid: str,
        mutated: str,
        gc: float,
        pm: int,
        se_count: int,
    ) -> None:
        """
        Write mutated FASTA, summary report, and append to log.
        """
        (run.artifacts / "mutated.fasta").write_text(f">{sid}\n{mutated}\n", encoding="utf-8")
        summary = {
            "sequence_id": sid,
            "length": len(mutated),
            "gc": gc,
            "point_mutations": pm,
            "structural_events": se_count,
        }
        (run.reports / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        with (run.logs / "runner.log").open("a", encoding="utf-8") as fh:
            fh.write("Run completed\n")

    def _emit(self, event: str, payload: dict[str, Any]) -> None:
        """
        Invoke progress callback if provided.
        """
        if self._progress:
            try:
                self._progress(event, payload)
            except Exception:
                # Progress callbacks must not crash the runner
                pass
