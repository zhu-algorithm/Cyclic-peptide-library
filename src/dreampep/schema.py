from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
from typing import Any


EVIDENCE_LEVELS = {"measured", "reported", "claimed", "inferred", "generated"}


@dataclass(slots=True)
class SourceRef:
    source_type: str
    identifier: str
    url: str
    retrieved_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )
    license_or_terms: str = "unknown"
    location: str = ""


@dataclass(slots=True)
class Assay:
    name: str = ""
    value: float | None = None
    unit: str = ""
    relation: str = "="
    evidence: str = "reported"
    conditions: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if self.evidence not in EVIDENCE_LEVELS:
            raise ValueError(f"Unsupported evidence level: {self.evidence}")


@dataclass(slots=True)
class PeptideRecord:
    record_id: str
    sequence: str
    target_id: str
    topology: str = "linear"
    helm: str = ""
    smiles: str = ""
    cyclization_bonds: list[dict[str, Any]] = field(default_factory=list)
    assay: Assay = field(default_factory=Assay)
    source: SourceRef | None = None
    split_group: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.record_id or not self.sequence or not self.target_id:
            raise ValueError("record_id, sequence and target_id are required")
        if self.topology not in {"linear", "cyclic", "bicyclic", "branched", "other"}:
            raise ValueError(f"Unsupported topology: {self.topology}")
        self.assay.validate()
        if self.source is None:
            raise ValueError("Every record requires a source")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def fingerprint(self) -> str:
        payload = json.dumps(self.to_dict(), sort_keys=True, ensure_ascii=False)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def record_from_dict(data: dict[str, Any]) -> PeptideRecord:
    data = dict(data)
    data["assay"] = Assay(**data.get("assay", {}))
    source = data.get("source")
    data["source"] = SourceRef(**source) if source else None
    record = PeptideRecord(**data)
    record.validate()
    return record

