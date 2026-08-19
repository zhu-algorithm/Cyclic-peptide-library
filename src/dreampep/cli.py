from __future__ import annotations

import argparse
import json
from pathlib import Path

from .generator import MarkovGenerator
from .ingest import ingest_reviewed_csv
from .io import read_jsonl, write_jsonl
from .library import (
    BillionThioetherLibrary,
    build_disulfide_library,
    write_thioether_library_chunk,
    write_thioether_library_manifest,
    write_vendor_csv,
)
from .pipeline import DesignPipeline
from .scoring import DevelopabilityScorer, LinearStudent


def main() -> None:
    parser = argparse.ArgumentParser(prog="dreampep")
    commands = parser.add_subparsers(dest="command", required=True)
    ingest = commands.add_parser("ingest")
    ingest.add_argument("input_csv")
    ingest.add_argument("output_jsonl")
    train = commands.add_parser("train-student")
    train.add_argument("dataset")
    train.add_argument("model")
    design = commands.add_parser("design")
    design.add_argument("dataset")
    design.add_argument("model")
    design.add_argument("--target", required=True)
    design.add_argument("--n", type=int, default=20)
    design.add_argument("--min-len", type=int, default=8)
    design.add_argument("--max-len", type=int, default=20)
    design.add_argument("--output", required=True)
    library = commands.add_parser("build-library")
    library.add_argument("--target", required=True)
    library.add_argument("--size", type=int, default=384)
    library.add_argument("--ring-length", type=int, default=10)
    library.add_argument("--seed", type=int, default=17)
    library.add_argument("--output", required=True)
    thioether = commands.add_parser("build-thioether-library")
    thioether.add_argument("--manifest", required=True)
    thioether.add_argument("--chunk-output")
    thioether.add_argument("--start", type=int, default=0)
    thioether.add_argument("--stop", type=int, default=0)
    thioether.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()

    if args.command == "ingest":
        records = ingest_reviewed_csv(args.input_csv)
        write_jsonl(args.output_jsonl, records)
        print(json.dumps({"records": len(records), "output": args.output_jsonl}))
    elif args.command == "train-student":
        records = read_jsonl(args.dataset)
        model = LinearStudent()
        model.fit(records)
        model.save(args.model)
        print(json.dumps({"records": len(records), "model": args.model}))
    elif args.command == "design":
        records = read_jsonl(args.dataset)
        generator = MarkovGenerator()
        generator.fit(records)
        pipeline = DesignPipeline(
            generator,
            {"student_activity": LinearStudent.load(args.model), "developability": DevelopabilityScorer()},
            {"student_activity": 0.65, "developability": 0.35},
        )
        candidates = pipeline.design(args.target, args.n, args.min_len, args.max_len)
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(json.dumps({"sequence": c.sequence, "scores": c.scores, "aggregate": c.aggregate}) for c in candidates) + "\n", encoding="utf-8")
        print(json.dumps({"candidates": len(candidates), "output": str(output)}))
    elif args.command == "build-library":
        rows = build_disulfide_library(args.size, args.ring_length, args.seed)
        write_vendor_csv(args.output, rows, args.target)
        print(json.dumps({"designs": len(rows), "chemistry": "disulfide", "output": args.output}))
    elif args.command == "build-thioether-library":
        virtual_library = BillionThioetherLibrary(seed=args.seed)
        write_thioether_library_manifest(args.manifest, virtual_library)
        exported = 0
        if args.chunk_output:
            if args.stop <= args.start:
                raise ValueError("--stop must be greater than --start when exporting a chunk")
            write_thioether_library_chunk(
                args.chunk_output, virtual_library, args.start, args.stop
            )
            exported = args.stop - args.start
        print(json.dumps({
            "designs": len(virtual_library),
            "chemistry": "N-chloroacetyl-to-Cys thioether",
            "manifest": args.manifest,
            "exported_rows": exported,
            "chunk_output": args.chunk_output,
        }))


if __name__ == "__main__":
    main()
