from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from .schema import Assay, PeptideRecord, SourceRef


REQUIRED_COLUMNS = {"sequence", "target_id", "source_identifier", "source_url"}


def _float_or_none(value: str) -> float | None:
    return float(value) if value.strip() else None


def ingest_reviewed_csv(path: str | Path) -> list[PeptideRecord]:
    """Import human-reviewed patent/paper rows; never scrapes or infers measurements."""
    records: list[PeptideRecord] = []
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
        for line_no, row in enumerate(reader, 2):
            seed = "|".join(
                [row["source_identifier"], row.get("location", ""), row["sequence"]]
            )
            record_id = row.get("record_id") or hashlib.sha256(seed.encode()).hexdigest()[:16]
            evidence = row.get("evidence", "reported").strip().lower()
            if evidence == "measured" and not row.get("assay_value", "").strip():
                raise ValueError(f"Line {line_no}: measured evidence requires assay_value")
            record = PeptideRecord(
                record_id=record_id,
                sequence=row["sequence"].strip(),
                target_id=row["target_id"].strip(),
                topology=row.get("topology", "linear").strip() or "linear",
                helm=row.get("helm", "").strip(),
                smiles=row.get("smiles", "").strip(),
                assay=Assay(
                    name=row.get("assay_name", "").strip(),
                    value=_float_or_none(row.get("assay_value", "")),
                    unit=row.get("assay_unit", "").strip(),
                    relation=row.get("assay_relation", "=").strip() or "=",
                    evidence=evidence,
                ),
                source=SourceRef(
                    source_type=row.get("source_type", "patent").strip(),
                    identifier=row["source_identifier"].strip(),
                    url=row["source_url"].strip(),
                    license_or_terms=row.get("license_or_terms", "unknown").strip(),
                    location=row.get("location", "").strip(),
                ),
                split_group=row.get("patent_family", "").strip() or row["source_identifier"].strip(),
                metadata={"review_status": "human-reviewed", "input_line": line_no},
            )
            record.validate()
            records.append(record)
    return records

