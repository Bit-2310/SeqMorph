from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, List, Optional

from mutation.types import StructuralEvent, OpName, ChooseBase
from mutation.structural import sample_structural_events, apply_structural_events_on_store
from mutation.point import apply_point_mutations_on_store
from mutation.context import ContextAwareChooser, ContextChooserConfig
# Optional: load JSON k-mer profiles (if you use models/)
try:
    from mutation.model_loader import load_model  # noqa: F401
    _HAS_MODEL_LOADER = True
except Exception:
    _HAS_MODEL_LOADER = False

# IMPORTANT: use a relative import for BaseStore to avoid top-level shadowing
from sequence_store.base import BaseStore


@dataclass
class MutationReport:
    """Minimal mutation report; extend as needed."""
    point_mutations: int
    structural_events: List[StructuralEvent]


class MutationEngine:
    """
    Store-first mutation engine.

    Design
    ------
    - Uses a single RNG (seedable) for reproducibility across all ops.
    - Structural ops operate directly on BaseStore backends (ChunkedStore recommended).
    - Point mutations use an injected `ChooseBase` strategy; by default we use a
      context-aware chooser with a CpG bias toggle + strength.
    - Optionally, you can pass a JSON model path (3-mer/5-mer tables). If not provided,
      a uniform baseline (with optional CpG bump) is used.

    Parameters
    ----------
    rng : random.Random | None
        Inject a pre-seeded RNG. If provided, `seed` is ignored.
    seed : int | None
        Seed for deterministic behavior when `rng` is not provided.
    chooser_cfg : ContextChooserConfig | None
        Configuration for the default context-aware chooser (CpG toggle + strength).
    default_choose_base : ChooseBase | None
        Custom default chooser strategy. If provided, overrides everything else.
    model_path : str | None
        Optional JSON model file with baseline + 3-mer/5-mer weights. Ignored if
        `default_choose_base` is provided or if model loader is unavailable.
    """

    def __init__(
        self,
        rng: Optional[random.Random] = None,
        seed: Optional[int] = None,
        chooser_cfg: Optional[ContextChooserConfig] = None,
        default_choose_base: Optional[ChooseBase] = None,
        model_path: Optional[str] = None,
    ):
        # RNG setup
        self.rng: random.Random = rng or random.Random()
        if rng is None and seed is not None:
            self.rng.seed(seed)

        # Default chooser config
        self._chooser_cfg = chooser_cfg or ContextChooserConfig()

        # Resolve default chooser
        if default_choose_base is not None:
            self._default_choose_base: ChooseBase = default_choose_base
        else:
            self._default_choose_base = self._build_default_chooser(model_path)

    # ---- internal ---------------------------------------------------------

    def _build_default_chooser(self, model_path: Optional[str]) -> ChooseBase:
        """
        Build a ContextAwareChooser either from a JSON model (if provided)
        or from a uniform baseline (still CpG-aware via chooser_cfg).
        """
        # If a model is requested but model_loader isn't available, fall back safely.
        if model_path and _HAS_MODEL_LOADER:
            mdl = load_model(model_path)  # type: ignore[name-defined]
            return ContextAwareChooser(
                self._chooser_cfg,
                mdl.ctx5 if mdl.k >= 5 else None,
                mdl.ctx3 if mdl.k >= 3 else None,
                mdl.baseline,
            )  # type: ignore[return-value]

        # Uniform baseline fallback (A,C,G,T uniform; N uniform to A/C/G/T)
        baseline = {
            "A": {"C": 1.0, "G": 1.0, "T": 1.0},
            "C": {"A": 1.0, "G": 1.0, "T": 1.0},
            "G": {"A": 1.0, "C": 1.0, "T": 1.0},
            "T": {"A": 1.0, "C": 1.0, "G": 1.0},
            "N": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
        }
        return ContextAwareChooser(self._chooser_cfg, None, None, baseline)  # type: ignore[return-value]

    # ---- point mutations --------------------------------------------------

    def mutate_points(
        self,
        store: BaseStore,
        positions: List[int],
        choose_base: Optional[ChooseBase] = None,
    ) -> int:
        """
        Mutate specific positions using the provided chooser or the engine's default.

        Parameters
        ----------
        store : BaseStore
            Sequence storage backend implementing get/set/len.
        positions : list[int]
            Zero-based positions to mutate.
        choose_base : ChooseBase | None
            Strategy for choosing the mutated base. If None, uses the engine default.

        Returns
        -------
        int
            Number of point mutations attempted (same as len(positions)).
        """
        policy: ChooseBase = choose_base or self._default_choose_base
        return apply_point_mutations_on_store(self.rng, store, positions, policy)

    # ---- structural mutations --------------------------------------------

    def mutate_structural(
        self,
        store: BaseStore,
        win_start: int,
        win_end: int,
        rate_pct: float,
        mean_seg_len: int = 500,
        mix: Optional[Dict[OpName, float]] = None,
    ) -> List[StructuralEvent]:
        """
        Sample and apply structural events within [win_start, win_end).

        Returns the list of applied events.
        """
        events = sample_structural_events(
            self.rng, win_start, win_end, rate_pct, mean_seg_len, mix  # type: ignore[arg-type]
        )
        apply_structural_events_on_store(store, events)
        return events

    # ---- combined convenience --------------------------------------------

    def mutate(
        self,
        store: BaseStore,
        *,
        point_positions: Optional[List[int]] = None,
        choose_base: Optional[ChooseBase] = None,
        struct_window: Optional[tuple[int, int]] = None,
        rate_pct: float = 0.0,
        mean_seg_len: int = 500,
        mix: Optional[Dict[OpName, float]] = None,
    ) -> MutationReport:
        """
        Convenience API to run point and/or structural mutations in one call.

        - If `point_positions` is provided, runs point mutations using `choose_base`
          or the engine's default context-aware chooser.
        - If `struct_window` is provided with `rate_pct > 0`, runs structural ops
          within the specified window.
        """
        pm = 0
        if point_positions:
            pm = self.mutate_points(store, point_positions, choose_base=choose_base)

        se: List[StructuralEvent] = []
        if struct_window and rate_pct > 0.0:
            s, e = struct_window
            se = self.mutate_structural(store, s, e, rate_pct, mean_seg_len, mix)

        return MutationReport(point_mutations=pm, structural_events=se)

    # ---- runtime tuning ---------------------------------------------------

    def set_default_chooser(self, chooser: ChooseBase) -> None:
        """Replace the engine's default point-mutation chooser."""
        self._default_choose_base = chooser

    def set_cpg_bias(self, enabled: bool, strength: Optional[float] = None) -> None:
        """
        Toggle CpG bias and optionally update strength on the default chooser config.
        Creates a new default chooser instance with the updated config.
        """
        cfg = self._chooser_cfg
        cfg = ContextChooserConfig(
            cpg_enabled=enabled,
            cpg_strength=strength if strength is not None else cfg.cpg_strength,
            alphabet=cfg.alphabet,
        )
        self._chooser_cfg = cfg
        # Rebuild default chooser using the previous model (if any) by passing None; the
        # chooser will reuse uniform baseline unless you created a custom chooser.
        self._default_choose_base = self._build_default_chooser(model_path=None)
