from __future__ import annotations

import csv
import math
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


@dataclass(slots=True)
class SelectionObservation:
    sequence: str
    round_index: int
    target_count: int
    counter_count: int
    enrichment: float
    selectivity: float


def simulate_rapid_selection(sequences: Iterable[str], target_score: Callable[[str], float],
    counter_score: Callable[[str], float] | None = None, rounds: int = 6,
    reads_per_round: int = 20_000, seed: int = 17) -> list[SelectionObservation]:
    """Research simulator for RaPID-like iterative selection, not wet-lab PDPS."""
    pool = list(dict.fromkeys(sequences))
    if not pool:
        return []
    counter_score = counter_score or (lambda _sequence: 0.0)
    rng = random.Random(seed)
    abundance = {sequence: 1.0 / len(pool) for sequence in pool}
    previous = Counter({sequence: 1 for sequence in pool})
    observations = []
    for round_index in range(1, rounds + 1):
        target_weight = {s: abundance[s] * math.exp(max(-8.0, min(8.0, target_score(s)))) for s in pool}
        counter_weight = {s: abundance[s] * math.exp(max(-8.0, min(8.0, counter_score(s)))) for s in pool}
        target_reads = Counter(rng.choices(pool, weights=target_weight.values(), k=reads_per_round))
        counter_reads = Counter(rng.choices(pool, weights=counter_weight.values(), k=max(1000, reads_per_round // 4)))
        for s in pool:
            observations.append(SelectionObservation(s, round_index, target_reads[s], counter_reads[s],
                (target_reads[s] + 1) / (previous[s] + 1), (target_reads[s] + 1) / (counter_reads[s] + 1)))
        retained = sorted(pool, key=lambda s: (target_reads[s] + 1) / (counter_reads[s] + 1), reverse=True)
        pool = retained[:max(32, len(retained) // 2)]
        total = sum(target_reads[s] + 1 for s in pool)
        abundance = {s: (target_reads[s] + 1) / total for s in pool}
        previous = target_reads
    return observations


def write_selection_csv(path: str | Path, rows: Iterable[SelectionObservation]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence", "round", "target_count", "counter_count", "enrichment", "selectivity"])
        writer.writeheader()
        for row in rows:
            writer.writerow({"sequence": row.sequence, "round": row.round_index, "target_count": row.target_count,
                "counter_count": row.counter_count, "enrichment": f"{row.enrichment:.6g}", "selectivity": f"{row.selectivity:.6g}"})
