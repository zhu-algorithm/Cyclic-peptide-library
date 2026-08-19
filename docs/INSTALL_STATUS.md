# Dependency installation status

Checked on 2026-08-19 on Windows with Python 3.12.13 and an NVIDIA Quadro P620
(4 GB). The repository commits are recorded below.

| Component | Commit/version | Status |
|---|---|---|
| PepINVENT | `e976dd18d68a6e0496ec0f70d7a4dc78d473ae6c` / package `0.0.1` | Source downloaded; local package installed in `.packages`; top-level import verified |
| PepExplainer | `3a64091469cce94df9f8f90113693f3598c09bbd` / package `ame 2023.9.20` | Source downloaded; local package installed in `.packages`; top-level `core` import verified |
| EvoBind | `1600488bc4160017054926d9174729e677796c0b` | Source downloaded; not runnable on native Windows because the official workflow assumes Linux, NVIDIA GPU, AlphaFold2 parameters and MSA databases |
| PyTorch | requested CPU wheel | Not installed: package-index traffic is blocked by the execution environment's proxy layer |
| DGL and model dependencies | not installed | Depend on a compatible PyTorch installation; package-index traffic is blocked |

The successful source-package imports do not mean model inference is ready.
PepINVENT requires its chemistry/ML dependencies and model assets.
PepExplainer requires PyTorch, DGL, RDKit, pandas, transformers and related
packages before its model modules can load.

To expose the locally installed source packages in PowerShell:

```powershell
. .\activate-local.ps1
```

For a complete production environment, use WSL2/Linux with separate Conda
environments. Do not combine PepINVENT's historical Python 3.8 / PyTorch 1.8.1
environment with PepExplainer's PyTorch 2.1.2 environment. EvoBind should have
its own environment and data volume.

