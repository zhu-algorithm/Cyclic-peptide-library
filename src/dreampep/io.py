from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from .schema import PeptideRecord, record_from_dict


def read_jsonl(path: str | Path) -> list[PeptideRecord]:
    records: list[PeptideRecord] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                records.append(record_from_dict(json.loads(line)))
            except Exception as exc:
                raise ValueError(f"{path}:{line_no}: {exc}") from exc
    return records


def write_jsonl(path: str | Path, records: Iterable[PeptideRecord]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for record in records:
            record.validate()
            handle.write(json.dumps(record.to_dict(), ensure_ascii=False) + "\n")

