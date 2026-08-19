# DreamPep cyclic-peptide library plan

## Recommended starting point

Build a target-conditioned, arrayed pilot library rather than attempting a
proprietary PDPS/RaPID reproduction. The first physical library should contain
384 individually identified disulfide-cyclized peptides in 96-well plates,
followed by a 960-member focused expansion after assay calibration. Maintain a
much larger virtual library in DreamPep and select compounds by diversity,
predicted binding, developability and uncertainty.

This format is deliberately compatible with common contract synthesis and
per-peptide LC-MS workflows. It is also easier to deconvolute and audit than a
pooled library. Disulfide closure is suitable for platform validation; later
rounds should compare thioether, lactam and head-to-tail variants because
disulfides may be reduced in some biological environments.

## Three complementary library tiers

| Tier | Format | Practical diversity | Primary purpose |
|---|---|---:|---|
| Virtual | HELM/SMILES records generated and scored by DreamPep | 10^5–10^7 | Explore sequence, residue and topology space |
| Arrayed physical | Individually synthesized and barcoded wells | 384–5,000 | Quantitative biochemical/cell assays and clean deconvolution |
| Encoded/display | OBOC, phage, mRNA display/RaPID-like | 10^5–10^13+ | De novo hit discovery when target and selection assay are mature |

The tiers should be connected: display hits are resynthesized as arrayed
compounds; arrayed assay data retrain DreamPep; DreamPep proposes the next
focused library.

## Pilot design specification

- Ring size: initially 8, 10 and 12 residues, stratified across plates.
- Baseline chemistry: terminal Cys–Cys disulfide with no internal cysteine.
- C terminus: amide for the default pilot; include matched acid controls where
  assay biology suggests it matters.
- Diversity: balance hydrophobic, polar, charged and turn-promoting residues;
  reject extreme hydrophobicity and charge before synthesis.
- Controls: linear counterparts, scrambled controls, known binder if available,
  no-peptide wells, plate duplicates and a small set of inter-plate standards.
- Identity: unique library ID, sequence, HELM, isomeric SMILES where possible,
  bridge definition, batch, well, source, QC and freeze/thaw history.
- QC: identity MS for every member; analytical LC/UPLC for every member in the
  screening set; record actual yield and purity rather than only pass/fail.
- Assay: establish a concentration-response confirmation step and orthogonal
  counter-screen before treating a primary signal as a hit.

## Design funnel

1. Define target construct, binding site hypothesis, assay modality and desired
   mechanism.
2. Generate a large virtual set with PepINVENT/DreamPep and known motif seeds.
3. Normalize every design to HELM plus explicit cyclization bonds.
4. Remove duplicates, reactive liabilities, aggregation-prone designs and
   close analogues of training/test records.
5. Rank by target/structure score, activity-student score, developability,
   novelty and uncertainty.
6. Select a max-min diverse set within score strata; reserve 10–20% for
   exploration rather than only taking the highest predicted scores.
7. Commission an arrayed pilot and import vendor QC into the provenance store.
8. Screen, confirm, retrain and construct focused SAR libraries around validated
   hits.

## Build a vendor-ready pilot CSV

The included command creates a chemistry-validation design. It is not yet
target-conditioned and should be rescored before ordering:

```powershell
$env:PYTHONPATH = "src"
python -m dreampep.cli build-library --target YOUR_TARGET --size 384 --ring-length 10 --output work/cyclic_library_384.csv
```

## Vendor fit

| Provider | Best fit based on public specifications |
|---|---|
| GenScript | Small-to-medium arrayed libraries, micro-scale or purified options, cyclic modifications and per-peptide QC |
| WuXi TIDES | Larger parallel campaigns, broad cyclization chemistry, unusual residues and progression toward process/CMC work |
| JPT | Plate-based peptide libraries and cyclization libraries; especially mature for immune/epitope applications |
| Biosynth/Pepscan | Arrayed libraries plus CLIPS constrained-peptide phage display and lead optimization |
| Creative Peptides | Custom arrayed cyclic libraries and outsourced screening; specifications should be confirmed in a formal proposal |
| Mimotopes | Multipin/SynPhase-style parallel synthesis heritage; confirm current commercial availability and cyclic-library QC directly |

Do not compare quotations on price alone. Require the same ring definition,
scale, purity basis, per-member MS/LC data, failed-sequence policy, plate map,
solubilization format, residual counterion specification and ownership/data-use
terms.

## IP boundary

PeptiDream states that PDPS combines Flexizyme, FIT and RaPID and is protected by
a portfolio built around Flexizyme-related patents. Public patent disclosure is
useful for understanding the field but is not permission to practice active
claims. The recommended pilot uses conventional arrayed chemical synthesis and
an independently implemented computational design pipeline. Obtain a
jurisdiction-specific freedom-to-operate review before commercial deployment,
especially before implementing Flexizyme/RaPID-like translation, CLIPS, encoded
libraries or proprietary scaffolds.

## Public sources

- PeptiDream PDPS: https://www.peptidream.com/en/science/pdps/
- PeptiDream publications: https://www.peptidream.com/en/science/paper/
- RaPID patent family entry: https://patents.google.com/patent/EP2492344A4/en
- GenScript libraries: https://www.genscript.com/peptide-library.html
- WuXi TIDES discovery synthesis: https://tides.wuxiapptec.com/services-solutions/peptide/discovery/custom-peptide-synthesis/
- JPT library design: https://www.jpt.com/support-contact/resources/peptide-library-pool-design/
- Biosynth libraries: https://www.biosynth.com/peptides/custom-peptides/libraries
- Biosynth constrained-peptide discovery: https://www.biosynth.com/peptides/peptide-manufacturing/lead-discovery-and-optimization
- Creative Peptides cyclic libraries: https://www.creative-peptides.com/services/cyclic-peptide-library-construction.html
- OBOC review: https://pmc.ncbi.nlm.nih.gov/articles/PMC7301614/
- Backbone macrocycle screening methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC7314982/

