from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .generator import Generator
from .scoring import Scorer


@dataclass(slots=True)
class RankedCandidate:
    sequence: str
    scores: dict[str, float]
    aggregate: float


class DesignPipeline:
    def __init__(self, generator: Generator, scorers: dict[str, Scorer], weights: dict[str, float] | None = None):
        self.generator = generator
        self.scorers = scorers
        self.weights = weights or {name: 1.0 for name in scorers}

    def design(self, target_id: str, n: int = 20, min_len: int = 8, max_len: int = 20) -> list[RankedCandidate]:
        sequences = self.generator.generate(target_id, n, min_len, max_len)
        matrices = {name: scorer.score(target_id, sequences) for name, scorer in self.scorers.items()}
        denominator = sum(abs(self.weights.get(name, 1.0)) for name in matrices) or 1.0
        ranked = []
        for index, sequence in enumerate(sequences):
            scores = {name: values[index] for name, values in matrices.items()}
            aggregate = sum(scores[name] * self.weights.get(name, 1.0) for name in scores) / denominator
            ranked.append(RankedCandidate(sequence, scores, aggregate))
        return sorted(ranked, key=lambda item: item.aggregate, reverse=True)

