# DreamPep

DreamPep is an independent, traceable research pipeline for target-conditioned
peptide design. It is an engineering scaffold—not a reproduction of
PeptiDream's proprietary PDPS platform and not a validated drug-discovery
system.

## Architecture

1. **Knowledge layer** — human-reviewed patent/paper examples in a strict
   provenance schema. Claimed, inferred, generated, reported and measured
   evidence are kept separate.
2. **Generator** — a JSON adapter is provided for PepINVENT. The included
   Markov generator is only a runnable smoke-test baseline.
3. **Activity student** — a JSON scorer adapter supports PepExplainer. The
   included linear student validates the training pipeline on small data.
4. **Target/structure scoring** — the same scorer adapter supports EvoBind,
   docking or structure services.
5. **Ranking** — weighted multi-objective ranking with separate score channels.

No PeptiDream code, model weights, confidential data or automated website
scraping is included.

## Quick start

Use Python 3.10 or newer. From this directory:

```powershell
$env:PYTHONPATH = "src"
python -m dreampep.cli ingest examples/reviewed_patent_examples.csv work/dataset.jsonl
python -m dreampep.cli train-student work/dataset.jsonl work/student.json
python -m dreampep.cli design work/dataset.jsonl work/student.json --target DEMO_TARGET --n 10 --output work/candidates.jsonl
python -m dreampep.cli build-library --target DEMO_TARGET --size 384 --ring-length 10 --output work/cyclic_library_384.csv
python -m dreampep.cli build-thioether-library --manifest work/thioether_manifest.json --chunk-output work/thioether_first_1000.csv --start 0 --stop 1000
python -m unittest discover -s tests -v
```

The bundled CSV is synthetic and exists only to test the workflow. It must not
be presented as PeptiDream experimental data.

## Billion-member thioether macrocycle library

`BillionThioetherLibrary` defines exactly 1,000,000,000 unique virtual 8--11mer
peptides. Each length receives approximately one quarter of the identifiers.
The C-terminal L-Cys is fixed for intramolecular thioether formation with an
N-terminal chloroacetyl group; the remaining positions use a 53-member default
catalogue spanning canonical L residues, D residues, N-methyl residues and
other representative noncanonical building blocks.

The library is deterministically indexable and should not be materialized as a
billion-row file. Use `build-thioether-library` to write the reproducibility
manifest and optionally export a manageable `[start, stop)` CSV chunk. The
curated catalogue is extensible and is not a claim that every known or
commercially available noncanonical amino acid is included. Confirm monomer,
translation/synthesis, protection, cyclization and analytical compatibility
before physical production. A human-readable workbook is provided at
`docs/dreampep_thioether_library.xlsx`.

## Patent curation rule

Only import a numeric label as `measured` or `reported` when the source gives a
specific assay, value, unit and example/table location. Use `claimed` for claim
scope and Markush enumeration. Never train the activity model on claimed or
generated rows as experimental positives. Keep all members of a patent family
in the same split using `split_group` to prevent leakage.

Required CSV columns:

```text
sequence,target_id,source_identifier,source_url
```

Recommended columns are demonstrated in
`examples/reviewed_patent_examples.csv`. Use HELM or SMILES for noncanonical
residues; the plain sequence field may contain a normalized display string.

## External model contract

`JsonCommandGenerator` sends one JSON object to stdin:

```json
{"target_id":"P04637","n":20,"min_len":8,"max_len":20}
```

It expects a JSON array of sequences or objects containing `sequence`.

`JsonCommandScorer` sends:

```json
{"target_id":"P04637","sequences":["ACDEFGHC"]}
```

It expects a JSON array of floats or objects containing `score`. This isolates
PepINVENT, PepExplainer and EvoBind environments from the orchestration layer.

## Production roadmap

- Curate patent examples with dual review and source snapshots.
- Represent modified/cyclic peptides in HELM plus canonical isomeric SMILES.
- Replace the baseline generator with a pinned PepINVENT environment.
- Train PepExplainer per assay family with patent-family grouped splits.
- Add EvoBind/structure scores, calibration and uncertainty estimates.
- Evaluate novelty against training data and patent claims.
- Require synthesis, binding, selectivity and developability experiments before
  interpreting any generated candidate as a lead.

## License and third-party components

DreamPep's original code is Apache-2.0. External repositories are not bundled.
Pin and record each external repository commit and retain its own license and
notices when integrating it. Patent publication does not grant freedom to
practice the claimed inventions; obtain a formal FTO review before commercial
use.
