# Source registry

This registry records architectural references only. A source appears here does
not mean its data or code is redistributed by DreamPep.

| Component | Source | Intended use | Bundled |
|---|---|---|---|
| Generator | https://github.com/MolecularAI/PepINVENT | External generation adapter | No |
| Activity student | https://github.com/zhaisilong/PepExplainer | External activity-scoring adapter | No |
| Target scoring | https://github.com/patrickbryant1/EvoBind | External target-conditioned scorer | No |
| Peptide utilities | https://github.com/novonordisk-research/pepfunn | Optional analysis utilities | No |
| Flexizyme research code | https://github.com/everyday847/flexizyme | Methodological reference | No |
| PeptiDream publications | https://www.peptidream.com/en/science/paper/ | Literature discovery | No |
| PeptiDream patents | https://patents.justia.com/assignee/peptidream-inc | Patent discovery only; verify against patent offices | No |
| RaPID patent family | https://patents.google.com/patent/EP2492344A4/en | IP/reference review | No |
| Peptide-library patent | https://patents.google.com/patent/US10711268B2/en | IP/reference review | No |

Locally retrieved repository commits and installation results are recorded in
`INSTALL_STATUS.md`.

For each integrated repository, add a lock entry containing repository URL,
commit SHA, retrieval date, license and local environment identifier. For each
curated assay row, retain the publication number, patent family, exact
example/table/paragraph, retrieval date and reviewer identity.
