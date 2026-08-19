from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PUBLIC_SOURCE_FIELDS = [
    "source_name", "source_kind", "record_id", "sequence", "target_id", "topology",
    "assay_name", "assay_value", "assay_unit", "evidence", "source_identifier",
    "source_url", "source_location", "license_or_terms", "access_date", "notes",
]


def normalize_public_catalog(input_csv: str | Path, output_csv: str | Path) -> int:
    """Normalize manually exported public rows without bypassing site controls."""
    count = 0
    with Path(input_csv).open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        missing = {"source_name", "source_url"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
        output = Path(output_csv)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8-sig") as target:
            writer = csv.DictWriter(target, fieldnames=PUBLIC_SOURCE_FIELDS)
            writer.writeheader()
            for line_no, row in enumerate(reader, 2):
                normalized = {field: (row.get(field) or "").strip() for field in PUBLIC_SOURCE_FIELDS}
                seed = "|".join([normalized["source_name"], normalized["source_identifier"], normalized["sequence"]])
                normalized["record_id"] = normalized["record_id"] or hashlib.sha256(seed.encode()).hexdigest()[:16]
                normalized["evidence"] = normalized["evidence"] or "catalog"
                normalized["source_location"] = normalized["source_location"] or f"CSV line {line_no}"
                writer.writerow(normalized)
                count += 1
    return count


def training_eligible(row: dict[str, str]) -> bool:
    """Only sequence-level reported/measured data are eligible for supervised training."""
    return bool(row.get("sequence", "").strip()) and row.get("evidence", "").lower() in {"reported", "measured"}
