from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
import json
import math
from pathlib import Path
import subprocess
from typing import Sequence

from .schema import PeptideRecord


def _features(sequence: str) -> dict[str, float]:
    seq = sequence.upper()
    values: dict[str, float] = {"bias": 1.0, "length": len(seq) / 20.0}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        values[f"aa:{aa}"] = seq.count(aa) / max(1, len(seq))
    for a, b in zip(seq, seq[1:]):
        values[f"di:{a}{b}"] = values.get(f"di:{a}{b}", 0.0) + 1.0 / max(1, len(seq) - 1)
    return values


class Scorer(ABC):
    @abstractmethod
    def score(self, target_id: str, sequences: Sequence[str]) -> list[float]:
        raise NotImplementedError


class LinearStudent(Scorer):
    """Interpretable baseline student; replace with PepExplainer through the same API."""

    def __init__(self):
        self.weights: dict[str, float] = defaultdict(float)

    def fit(self, records: Sequence[PeptideRecord], epochs: int = 100, lr: float = 0.15) -> None:
        train = [r for r in records if r.assay.value is not None and r.assay.evidence in {"measured", "reported"}]
        if not train:
            raise ValueError("No measured/reported assay rows with numeric values")
        raw = [float(r.assay.value) for r in train]
        midpoint = sorted(raw)[len(raw) // 2]
        for _ in range(epochs):
            for record in train:
                # Lower concentration-like values are treated as stronger activity.
                label = 1.0 if float(record.assay.value) <= midpoint else 0.0
                feats = _features(record.sequence)
                z = sum(self.weights[k] * v for k, v in feats.items())
                pred = 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z))))
                for key, value in feats.items():
                    self.weights[key] += lr * (label - pred) * value

    def score(self, target_id: str, sequences: Sequence[str]) -> list[float]:
        del target_id
        output = []
        for sequence in sequences:
            z = sum(self.weights[k] * v for k, v in _features(sequence).items())
            output.append(1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z)))))
        return output

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(dict(self.weights), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "LinearStudent":
        model = cls()
        model.weights.update(json.loads(Path(path).read_text(encoding="utf-8")))
        return model


class JsonCommandScorer(Scorer):
    """Adapter for PepExplainer, EvoBind, docking, or a structure service."""

    def __init__(self, command: list[str]):
        self.command = command

    def score(self, target_id: str, sequences: Sequence[str]) -> list[float]:
        payload = {"target_id": target_id, "sequences": list(sequences)}
        result = subprocess.run(
            self.command, input=json.dumps(payload), text=True,
            capture_output=True, check=True,
        )
        values = json.loads(result.stdout)
        return [float(v["score"] if isinstance(v, dict) else v) for v in values]


class DevelopabilityScorer(Scorer):
    def score(self, target_id: str, sequences: Sequence[str]) -> list[float]:
        del target_id
        scores = []
        for seq in sequences:
            length_term = max(0.0, 1.0 - abs(len(seq) - 14) / 14)
            hydrophobic = sum(seq.count(x) for x in "AILMFWVY") / max(1, len(seq))
            charge = abs(sum(seq.count(x) for x in "KR") - sum(seq.count(x) for x in "DE"))
            scores.append(max(0.0, min(1.0, 0.5 * length_term + 0.4 * (1 - abs(hydrophobic - 0.4)) - 0.03 * charge)))
        return scores

