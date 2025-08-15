# src/cli.py
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
import json
from .sequence_store.base import BaseStore
from .runner import Runner
from .mutation.context import ContextChooserConfig


# -------------------------- Banner / Colors --------------------------- #

def _supports_color() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()

def _c(code: str) -> str:
    return code if _supports_color() else ""

RESET = _c("\033[0m")
BOLD = _c("\033[1m")
YELLOW = _c("\033[93m")
CYAN = _c("\033[96m")
MAGENTA = _c("\033[95m")
GREEN = _c("\033[92m")
BLUE = _c("\033[94m")
WHITE = _c("\033[97m")

BORDER = "·······································································"

def print_logo(multicolor: bool = True) -> None:
    lines = [
        BORDER,
        ":███████╗███████╗ ██████╗ ███╗   ███╗ ██████╗ ██████╗ ██████╗ ██╗  ██╗:",
        ":██╔════╝██╔════╝██╔═══██╗████╗ ████║██╔═══██╗██╔══██╗██╔══██╗██║  ██║:",
        ":███████╗█████╗  ██║   ██║██╔████╔██║██║   ██║██████╔╝██████╔╝███████║:",
        ":╚════██║██╔══╝  ██║▄▄ ██║██║╚██╔╝██║██║   ██║██╔══██╗██╔═══╝ ██╔══██║:",
        ":███████║███████╗╚██████╔╝██║ ╚═╝ ██║╚██████╔╝██║  ██║██║     ██║  ██║:",
        ":╚══════╝╚══════╝ ╚══▀▀═╝ ╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝:",
        BORDER,
    ]
    if not _supports_color() or not multicolor:
        print("\n" + "\n".join(lines) + "\n")
        return

    palette = [YELLOW, CYAN, MAGENTA, BLUE, GREEN, MAGENTA, YELLOW]
    colored = []
    for i, line in enumerate(lines):
        color = palette[i % len(palette)]
        if line == BORDER:
            colored.append(f"{YELLOW}{line}{RESET}")
        else:
            # Add a subtle gradient effect by tinting the left token on each line
            l = line
            l = l.replace("███████╗", f"{CYAN}███████╗{RESET}")
            l = l.replace("██╔════╝", f"{MAGENTA}██╔════╝{RESET}")
            colored.append(color + l + RESET)
    print("\n" + "\n".join(colored) + "\n")


# -------------------------- IO helpers -------------------------------- #

def printf(msg: str = "", *, style: str | None = None) -> None:
    if style == "ok":
        print(f"{GREEN}{msg}{RESET}" if _supports_color() else msg)
    elif style == "err":
        print(f"{MAGENTA}{msg}{RESET}" if _supports_color() else msg)
    elif style == "bold":
        print(f"{BOLD}{msg}{RESET}" if _supports_color() else msg)
    else:
        print(msg)

def ask(prompt: str, default: Optional[str] = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        try:
            val = input(f"{prompt}{suffix}: ").strip()
        except EOFError:
            val = ""
        if val:
            return val
        if default is not None:
            return default
        printf("Please enter a value.", style="err")

def ask_yes_no(prompt: str, default: bool = True) -> bool:
    d = "Y/n" if default else "y/N"
    while True:
        val = ask(f"{prompt} ({d})", "").lower()
        if val == "": return default
        if val in ("y", "yes"): return True
        if val in ("n", "no"): return False
        printf("Please answer y/n.", style="err")

def ask_int(prompt: str, default: Optional[int] = None, *, min_val: Optional[int] = None, max_val: Optional[int] = None) -> int:
    while True:
        raw = ask(prompt, str(default) if default is not None else None)
        try:
            v = int(raw)
            if min_val is not None and v < min_val:
                printf(f"Value must be >= {min_val}", style="err"); continue
            if max_val is not None and v > max_val:
                printf(f"Value must be <= {max_val}", style="err"); continue
            return v
        except ValueError:
            printf("Please enter an integer.", style="err")

def ask_float(prompt: str, default: Optional[float] = None, *, min_val: Optional[float] = None, max_val: Optional[float] = None) -> float:
    while True:
        raw = ask(prompt, str(default) if default is not None else None)
        try:
            v = float(raw)
            if min_val is not None and v < min_val:
                printf(f"Value must be >= {min_val}", style="err"); continue
            if max_val is not None and v > max_val:
                printf(f"Value must be <= {max_val}", style="err"); continue
            return v
        except ValueError:
            printf("Please enter a number.", style="err")

def parse_positions(raw: str) -> List[int]:
    parts = [p.strip() for p in raw.replace(",", " ").split() if p.strip()]
    out: List[int] = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            raise ValueError(f"Invalid position: {p}")
    return out

def read_fasta_first_record(path: Path) -> Tuple[str, str]:
    text = path.read_text(encoding="utf-8").splitlines()
    if not text or not text[0].startswith(">"):
        raise ValueError("Invalid FASTA: first line must start with '>'")
    header = text[0][1:].strip() or path.stem
    seq_lines: List[str] = []
    for ln in text[1:]:
        if ln.startswith(">"):
            break
        s = ln.strip()
        if s:
            seq_lines.append(s)
    seq = "".join(seq_lines).upper().replace("U", "T")
    return header, seq


# -------------------------- Interactive flow -------------------------- #

def choose_sequence_source() -> Tuple[str, str]:
    printf("Choose sequence source:", style="bold")
    printf("  [1] FASTA file path")
    printf("  [2] Database (stub): enter ID and sequence manually")
    printf("  [3] Manual entry")

    while True:
        choice = ask("Select 1/2/3", "1")
        if choice == "1":
            p = Path(ask("FASTA file path"))
            if not p.exists():
                printf("File not found, try again.", style="err")
                continue
            sid, seq = read_fasta_first_record(p)
            sid = ask("Sequence ID", sid)
            return sid, seq
        elif choice == "2":
            sid = ask("Database sequence ID")
            seq = ask("Sequence (ACGT…; U will be mapped to T)").upper().replace("U", "T")
            return sid, seq
        elif choice == "3":
            sid = ask("Sequence ID", "seq1")
            seq = ask("Sequence (ACGT…; U will be mapped to T)").upper().replace("U", "T")
            return sid, seq
        else:
            printf("Please choose 1, 2, or 3.", style="err")


def progress_printer(event: str, payload: Dict[str, Any]) -> None:
    if event == "run_created":
        printf(f"[+] Run created: {payload.get('run_id')}", style="ok")
    elif event == "run_completed":
        printf(f"[✓] Run completed: {payload.get('run_id')}", style="ok")
    elif event == "run_failed":
        printf(f"[x] Run failed: {payload.get('run_id')} -> {payload.get('error')}", style="err")


def main() -> int:
    # 1) Banner
    print_logo(multicolor=True)

    # 2) Outputs root (optional override)
    printf("Output directory can be overridden via SEQMORPH_OUTPUTS_DIR or entered below.")
    out_root_in = ask("Outputs root (blank to use default ./outputs)", "")
    outputs_root = out_root_in or None

    # 3) Sequence
    sid, seq = choose_sequence_source()

    # 4) Parameters
    printf("\nParameters", style="bold")
    seed = ask_int("Seed", 42)
    use_cpg = ask_yes_no("Apply CpG C->T bias?", True)
    cpg_strength = ask_float("CpG strength multiplier", 3.0, min_val=0.01)

    pos_raw = ask("Point mutation positions (e.g., '1,3,5' or blank for none)", "")
    positions: Optional[List[int]] = None
    if pos_raw.strip():
        try:
            positions = parse_positions(pos_raw)
        except ValueError as e:
            printf(str(e), style="err")
            return 1

    do_struct = ask_yes_no("Configure structural window?", False)
    struct_window: Optional[Tuple[int, int]] = None
    rate_pct = 0.0
    mean_seg_len = 500
    if do_struct:
        start = ask_int("Window start (0-based index)", 0, min_val=0)
        end = ask_int("Window end (exclusive)", max(1, len(seq)), min_val=1)
        if end <= start:
            printf("End must be greater than start.", style="err")
            return 1
        struct_window = (start, end)
        rate_pct = ask_float("Structural rate % (0..100)", 1.0, min_val=0.0, max_val=100.0)
        mean_seg_len = ask_int("Mean segment length", 500, min_val=1)

    model_path = ask("Model JSON path (blank for built-in baseline+CpG)", "")
    model_path = model_path or None

    label = ask("Run label", "mutate_and_analyze")

    # 5) Preview & confirm
    printf("\nSummary", style="bold")
    summary_preview = {
        "sequence_id": sid, "length": len(seq),
        "seed": seed, "cpg_enabled": use_cpg, "cpg_strength": cpg_strength,
        "positions": positions, "struct_window": struct_window,
        "rate_pct": rate_pct, "mean_seg_len": mean_seg_len,
        "model_path": model_path, "label": label,
        "outputs_root": outputs_root or os.environ.get("SEQMORPH_OUTPUTS_DIR", "./outputs"),
    }
    print(json.dumps(summary_preview, indent=2))
    if not ask_yes_no("Proceed?", True):
        printf("Aborted.")
        return 0

    # 6) Run
    chooser_cfg = ContextChooserConfig(cpg_enabled=use_cpg, cpg_strength=cpg_strength)
    runner = Runner(outputs_root=outputs_root, progress=progress_printer)

    try:
        runner.add_sequence(sid, seq)
    except ValueError as e:
        printf(f"Invalid sequence: {e}", style="err")
        return 1

    printf("\nRunning...", style="bold")
    res = runner.mutate_and_analyze(
        sid,
        point_positions=positions,
        struct_window=struct_window,
        rate_pct=rate_pct,
        mean_seg_len=mean_seg_len,
        seed=seed,
        chooser_cfg=chooser_cfg,
        model_path=model_path,
        label=label,
    )

    # 7) Output summary + where files are saved
    out_dir = Path(runner.om.root) / "runs" / res.run_id
    printf("\nResult", style="bold")
    print(json.dumps({
        "run_id": res.run_id,
        "sequence_id": res.sequence_id,
        "length": res.length,
        "gc": res.gc,
        "point_mutations": res.point_mutations,
        "structural_events": res.structural_events,
        "status": res.status,
        "message": res.message,
        "outputs_path": str(out_dir),
    }, indent=2))

    log_path = out_dir / "logs" / "runner.log"
    if log_path.exists():
        printf("\n--- runner.log ---", style="bold")
        try:
            print(log_path.read_text(encoding="utf-8"))
        except Exception:
            pass

    printf(f"Saved outputs to: {out_dir}", style="ok")
    return 0 if res.status == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
