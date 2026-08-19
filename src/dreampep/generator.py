from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
import json
import random
import subprocess
from typing import Sequence

from .schema import PeptideRecord


CANONICAL = "ACDEFGHIKLMNPQRSTVWY"


class Generator(ABC):
    @abstractmethod
    def generate(self, target_id: str, n: int, min_len: int, max_len: int) -> list[str]:
        raise NotImplementedError


class MarkovGenerator(Generator):
    """Small deterministic baseline for pipeline validation, not a discovery model."""

    def __init__(self, seed: int = 7):
        self.rng = random.Random(seed)
        self.transitions: dict[str, list[str]] = defaultdict(list)

    def fit(self, records: Sequence[PeptideRecord]) -> None:
        for record in records:
            seq = "".join(x for x in record.sequence.upper() if x in CANONICAL)
            for left, right in zip("^" + seq, seq + "$"):
                self.transitions[left].append(right)

    def generate(self, target_id: str, n: int, min_len: int, max_len: int) -> list[str]:
        del target_id
        results: list[str] = []
        for _ in range(n * 10):
            seq, token = "", "^"
            while len(seq) < max_len:
                choices = self.transitions.get(token) or list(CANONICAL)
                nxt = self.rng.choice(choices)
                if nxt == "$":
                    if len(seq) >= min_len:
                        break
                    continue
                seq += nxt
                token = nxt
            if min_len <= len(seq) <= max_len and seq not in results:
                results.append(seq)
            if len(results) == n:
                break
        return results


class JsonCommandGenerator(Generator):
    """Adapter for PepINVENT or another generator exposed as a JSON subprocess."""

    def __init__(self, command: list[str]):
        self.command = command

    def generate(self, target_id: str, n: int, min_len: int, max_len: int) -> list[str]:
        payload = {"target_id": target_id, "n": n, "min_len": min_len, "max_len": max_len}
        result = subprocess.run(
            self.command,
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        output = json.loads(result.stdout)
        return [item["sequence"] if isinstance(item, dict) else item for item in output]
