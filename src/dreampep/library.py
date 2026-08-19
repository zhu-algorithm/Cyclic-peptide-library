from __future__ import annotations

import csv
import itertools
import json
import math
import random
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


RESIDUE_GROUPS = {
    "hydrophobic": "AILMFVWY",
    "polar": "STNQ",
    "positive": "KRH",
    "negative": "DE",
    "turn": "GP",
}

CANONICAL_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
PEPTIDE_LENGTH_RANGE = range(8, 12)
THIOETHER_LIBRARY_SIZE = 1_000_000_000


@dataclass(frozen=True)
class BuildingBlock:
    """A residue building block represented without lossy one-letter encoding."""

    code: str
    name: str
    residue_class: str
    stereo: str


# The catalogue is intentionally explicit and can be replaced with a vendor list.
# L-Cys is reserved for the C-terminal cyclisation residue and excluded internally.
_NATURAL_NAMES = {
    "A": "alanine", "D": "aspartic acid", "E": "glutamic acid", "F": "phenylalanine",
    "G": "glycine", "H": "histidine", "I": "isoleucine", "K": "lysine",
    "L": "leucine", "M": "methionine", "N": "asparagine", "P": "proline",
    "Q": "glutamine", "R": "arginine", "S": "serine", "T": "threonine",
    "V": "valine", "W": "tryptophan", "Y": "tyrosine",
}


def default_thioether_building_blocks() -> tuple[BuildingBlock, ...]:
    """Natural and representative commercially accessible non-natural residues."""
    natural = tuple(
        BuildingBlock(aa, name, "natural", "L" if aa != "G" else "achiral")
        for aa, name in _NATURAL_NAMES.items()
    )
    d_residues = tuple(
        BuildingBlock(f"d{aa}", f"D-{name}", "non-natural", "D")
        for aa, name in _NATURAL_NAMES.items() if aa != "G"
    )
    non_natural = (
        BuildingBlock("NMeA", "N-methyl-L-alanine", "non-natural", "L"),
        BuildingBlock("NMeF", "N-methyl-L-phenylalanine", "non-natural", "L"),
        BuildingBlock("NMeG", "sarcosine", "non-natural", "achiral"),
        BuildingBlock("Aib", "2-aminoisobutyric acid", "non-natural", "achiral"),
        BuildingBlock("BAl", "beta-alanine", "non-natural", "achiral"),
        BuildingBlock("Orn", "L-ornithine", "non-natural", "L"),
        BuildingBlock("Dab", "L-2,4-diaminobutyric acid", "non-natural", "L"),
        BuildingBlock("Nle", "L-norleucine", "non-natural", "L"),
        BuildingBlock("Cha", "L-cyclohexylalanine", "non-natural", "L"),
        BuildingBlock("Nal1", "L-1-naphthylalanine", "non-natural", "L"),
        BuildingBlock("Nal2", "L-2-naphthylalanine", "non-natural", "L"),
        BuildingBlock("Fpa", "4-fluoro-L-phenylalanine", "non-natural", "L"),
        BuildingBlock("ClF", "4-chloro-L-phenylalanine", "non-natural", "L"),
        BuildingBlock("Bpa", "4-benzoyl-L-phenylalanine", "non-natural", "L"),
        BuildingBlock("Hph", "L-homophenylalanine", "non-natural", "L"),
        BuildingBlock("Cit", "L-citrulline", "non-natural", "L"),
    )
    return natural + d_residues + non_natural


@dataclass(frozen=True)
class ThioetherPeptide:
    library_id: str
    residues: tuple[str, ...]
    length: int
    n_terminal: str = "N-chloroacetyl"
    c_terminal: str = "carboxamide"
    cyclization: str = "N-chloroacetyl-to-C-terminal-Cys thioether"

    @property
    def sequence(self) -> str:
        return "-".join(self.residues)


class BillionThioetherLibrary:
    """Indexable 8--11-mer thioether virtual library with one billion members.

    Only the variable positions are decoded.  The final residue is invariant L-Cys,
    whose thiol displaces chloride from the N-terminal chloroacetyl group.  An affine
    length strata and affine permutations spread the billion IDs evenly across lengths
    and chemical space while preserving a one-to-one mapping and reproducibility.
    """

    def __init__(
        self,
        building_blocks: Sequence[BuildingBlock] | None = None,
        size: int = THIOETHER_LIBRARY_SIZE,
        seed: int = 17,
    ) -> None:
        self.building_blocks = tuple(building_blocks or default_thioether_building_blocks())
        if len(self.building_blocks) < 2:
            raise ValueError("at least two variable building blocks are required")
        codes = [block.code for block in self.building_blocks]
        if len(codes) != len(set(codes)):
            raise ValueError("building-block codes must be unique")
        self.lengths = tuple(PEPTIDE_LENGTH_RANGE)
        self.capacity_by_length = tuple(len(codes) ** (length - 1) for length in self.lengths)
        self.capacity = sum(self.capacity_by_length)
        balanced_capacity = min(self.capacity_by_length) * len(self.lengths)
        if not 1 <= size <= balanced_capacity:
            raise ValueError(f"size must be between 1 and balanced capacity {balanced_capacity}")
        self.size = size
        rng = random.Random(seed)
        multipliers: list[int] = []
        offsets: list[int] = []
        for modulus in self.capacity_by_length:
            multiplier = rng.randrange(2, modulus)
            while math.gcd(multiplier, modulus) != 1:
                multiplier += 1
            multipliers.append(multiplier)
            offsets.append(rng.randrange(modulus))
        self._multipliers = tuple(multipliers)
        self._offsets = tuple(offsets)

    def __len__(self) -> int:
        return self.size

    def _decode_space_index(self, space_index: int) -> tuple[str, ...]:
        for length, count in zip(self.lengths, self.capacity_by_length):
            if space_index < count:
                variable: list[str] = []
                radix = len(self.building_blocks)
                for _ in range(length - 1):
                    space_index, digit = divmod(space_index, radix)
                    variable.append(self.building_blocks[digit].code)
                return tuple(reversed(variable)) + ("C",)
            space_index -= count
        raise IndexError("space index outside library capacity")

    def _decode_length_index(self, length: int, value: int) -> tuple[str, ...]:
        variable: list[str] = []
        radix = len(self.building_blocks)
        for _ in range(length - 1):
            value, digit = divmod(value, radix)
            variable.append(self.building_blocks[digit].code)
        return tuple(reversed(variable)) + ("C",)

    def peptide(self, index: int) -> ThioetherPeptide:
        if not 0 <= index < self.size:
            raise IndexError(f"library index must be in [0, {self.size})")
        length_slot = index % len(self.lengths)
        local_index = index // len(self.lengths)
        modulus = self.capacity_by_length[length_slot]
        space_index = (
            local_index * self._multipliers[length_slot] + self._offsets[length_slot]
        ) % modulus
        residues = self._decode_length_index(self.lengths[length_slot], space_index)
        return ThioetherPeptide(f"DPT-{index + 1:010d}", residues, len(residues))

    def iter_range(self, start: int = 0, stop: int | None = None) -> Iterator[ThioetherPeptide]:
        stop = self.size if stop is None else stop
        if not 0 <= start <= stop <= self.size:
            raise ValueError("invalid library range")
        for index in range(start, stop):
            yield self.peptide(index)

    def sample(self, count: int, seed: int = 17) -> list[ThioetherPeptide]:
        if not 0 <= count <= self.size:
            raise ValueError("sample count outside library size")
        rng = random.Random(seed)
        return [self.peptide(index) for index in rng.sample(range(self.size), count)]

    def manifest(self) -> dict[str, object]:
        return {
            "library_type": "virtual thioether macrocycle",
            "size": self.size,
            "combinatorial_capacity": self.capacity,
            "lengths": list(self.lengths),
            "variable_positions": "positions 1..N-1",
            "fixed_residue": "C-terminal L-Cys",
            "n_terminal": "N-chloroacetyl",
            "c_terminal": "carboxamide",
            "cyclization": "intramolecular thioether",
            "building_blocks": [asdict(block) for block in self.building_blocks],
            "index_mapping": {
                "type": "length-stratified affine permutations",
                "length_assignment": "index modulo 4",
                "multipliers": list(self._multipliers),
                "offsets": list(self._offsets),
                "moduli": list(self.capacity_by_length),
            },
        }


def write_thioether_library_manifest(path: str | Path, library: BillionThioetherLibrary) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(library.manifest(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_thioether_library_chunk(
    path: str | Path,
    library: BillionThioetherLibrary,
    start: int,
    stop: int,
) -> None:
    """Export a manageable slice; exporting all billion rows is intentionally avoided."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "library_id", "sequence", "length", "n_terminal", "c_terminal", "cyclization"
        ))
        writer.writeheader()
        for peptide in library.iter_range(start, stop):
            writer.writerow({
                "library_id": peptide.library_id,
                "sequence": peptide.sequence,
                "length": peptide.length,
                "n_terminal": peptide.n_terminal,
                "c_terminal": peptide.c_terminal,
                "cyclization": peptide.cyclization,
            })


@dataclass(frozen=True)
class BackboneTemplate:
    name: str
    phi: float
    psi: float


EXTENDED_TEMPLATE = BackboneTemplate("extended", -60.0, 120.0)
BETA_TEMPLATE = BackboneTemplate("beta", -135.0, 135.0)


@dataclass(frozen=True)
class WorkflowStep:
    """One externally executable and auditable structure-design step."""

    name: str
    command: tuple[str, ...]
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]


def validate_peptide_sequence(sequence: str, min_length: int = 8, max_length: int = 11) -> str:
    sequence = sequence.strip().upper()
    if not min_length <= len(sequence) <= max_length:
        raise ValueError(f"peptide length must be between {min_length} and {max_length}")
    invalid = sorted(set(sequence) - set(CANONICAL_AMINO_ACIDS))
    if invalid:
        raise ValueError(f"unsupported amino-acid codes: {''.join(invalid)}")
    return sequence


def enumerate_peptides(
    lengths: Iterable[int] = PEPTIDE_LENGTH_RANGE,
    alphabet: str = CANONICAL_AMINO_ACIDS,
    limit: int | None = None,
) -> Iterator[str]:
    """Lazily enumerate a sequence space without materialising it in memory."""
    if limit is not None and limit < 0:
        raise ValueError("limit must be non-negative")
    alphabet = "".join(dict.fromkeys(alphabet.upper()))
    if not alphabet or set(alphabet) - set(CANONICAL_AMINO_ACIDS):
        raise ValueError("alphabet must contain only canonical amino acids")
    produced = 0
    for length in lengths:
        if length not in PEPTIDE_LENGTH_RANGE:
            raise ValueError("all peptide lengths must be between 8 and 11")
        for residues in itertools.product(alphabet, repeat=length):
            if limit is not None and produced >= limit:
                return
            yield "".join(residues)
            produced += 1


def choose_backbone_template(sequence: str, bulky_fraction: float = 0.45) -> BackboneTemplate:
    """Choose the initial template; callers may force beta after a failed build."""
    sequence = validate_peptide_sequence(sequence)
    bulky = sum(residue in "FWYRK" for residue in sequence) / len(sequence)
    return BETA_TEMPLATE if bulky >= bulky_fraction else EXTENDED_TEMPLATE


def write_chimerax_batch(
    sequences: Iterable[str],
    output_dir: str | Path,
    *,
    chimerax_executable: str = "ChimeraX",
    batch_size: int = 1000,
) -> list[Path]:
    """Create ChimeraX 1.8 command files and JSONL metadata for sequence builds.

    Commands create an initial model, set uniform backbone torsions, add/optimise
    hydrogens, and save PDB files.  Failed extended builds can be regenerated by
    calling this function again with sequences selected from the run log; bulky
    sequences are assigned the beta template up front and every choice is recorded.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    root = Path(output_dir)
    scripts = root / "scripts"
    structures = root / "initial_structures"
    scripts.mkdir(parents=True, exist_ok=True)
    structures.mkdir(parents=True, exist_ok=True)
    metadata_path = root / "build_metadata.jsonl"
    script_paths: list[Path] = []
    handle = None
    with metadata_path.open("w", encoding="utf-8") as metadata:
        for index, raw_sequence in enumerate(sequences, 1):
            sequence = validate_peptide_sequence(raw_sequence)
            template = choose_backbone_template(sequence)
            if (index - 1) % batch_size == 0:
                if handle is not None:
                    handle.close()
                script_path = scripts / f"build_{len(script_paths) + 1:05d}.cxc"
                script_paths.append(script_path)
                handle = script_path.open("w", encoding="utf-8")
            peptide_id = f"DP-{index:09d}"
            pdb_path = (structures / f"{peptide_id}.pdb").resolve()
            # One phi,psi pair is reused for every residue by ChimeraX.
            handle.write(
                f'build start peptide "{peptide_id}" {sequence} '
                f"{template.phi},{template.psi} chainId A rotLib Dunbrack\n"
            )
            handle.write("addh\nminimize\n")
            handle.write(f"save \"{pdb_path}\" format pdb\nclose all\n")
            metadata.write(json.dumps({
                "library_id": peptide_id,
                "sequence": sequence,
                "length": len(sequence),
                "template": asdict(template),
                "structure": str(pdb_path),
                "status": "planned",
                "chimerax_version": "1.8",
            }) + "\n")
    if handle is not None:
        handle.close()
    launcher = root / "run_chimerax.jsonl"
    with launcher.open("w", encoding="utf-8") as stream:
        for script in script_paths:
            stream.write(json.dumps({"command": [chimerax_executable, "--nogui", str(script.resolve())]}) + "\n")
    return script_paths


def conformer_sampling_workflow(
    initial_structure_dir: str | Path,
    output_dir: str | Path,
    *,
    conformers_per_peptide: int = 32,
    sampler_command: Sequence[str] = ("python", "sample_macrocycle_conformers.py"),
) -> list[WorkflowStep]:
    """Describe hydrogen-aware macrocycle sampling and low-energy optimisation.

    The sampler is deliberately configurable because ChimeraX prepares/minimises
    the starting structures but does not provide a headless exhaustive macrocycle
    conformer search.  A Rosetta, RDKit, OpenMM, or laboratory sampler wrapper can
    implement the command contract and write ranked SDF/PDB ensembles.
    """
    if conformers_per_peptide < 1:
        raise ValueError("conformers_per_peptide must be positive")
    initial = Path(initial_structure_dir).resolve()
    root = Path(output_dir).resolve()
    sampled, optimised = root / "sampled", root / "low_energy"
    return [
        WorkflowStep("sample_conformers", tuple(sampler_command) + (
            "--input", str(initial), "--output", str(sampled),
            "--num-conformers", str(conformers_per_peptide), "--cyclic",
        ), (str(initial),), (str(sampled),)),
        WorkflowStep("optimise_and_rank", tuple(sampler_command) + (
            "--input", str(sampled), "--output", str(optimised),
            "--optimise", "--rank-by-energy",
        ), (str(sampled),), (str(optimised),)),
    ]


def rfpeptides_workflow(
    receptor_pdb: str | Path,
    output_dir: str | Path,
    *,
    min_length: int = 8,
    max_length: int = 11,
    backbones: int = 100,
    sequences_per_backbone: int = 8,
    rfdiffusion_command: Sequence[str] = ("python", "run_inference.py"),
    proteinmpnn_command: Sequence[str] = ("python", "protein_mpnn_run.py"),
    rosetta_command: Sequence[str] = ("rosetta_scripts",),
) -> list[WorkflowStep]:
    """Build an RFpeptides-style RFdiffusion -> MPNN -> Rosetta execution plan."""
    if min_length not in PEPTIDE_LENGTH_RANGE or max_length not in PEPTIDE_LENGTH_RANGE or min_length > max_length:
        raise ValueError("length interval must fall within 8..11")
    if backbones < 1 or sequences_per_backbone < 1:
        raise ValueError("sampling counts must be positive")
    receptor = Path(receptor_pdb).resolve()
    root = Path(output_dir).resolve()
    backbones_dir, sequences_dir, refined_dir = (root / name for name in ("backbones", "sequences", "refined"))
    return [
        WorkflowStep("rfdiffusion", tuple(rfdiffusion_command) + (
            f"inference.input_pdb={receptor}", f"inference.output_prefix={backbones_dir / 'macrocycle'}",
            f"inference.num_designs={backbones}", f"contigmap.contigs=[{min_length}-{max_length}]",
            "inference.cyclic=True",
        ), (str(receptor),), (str(backbones_dir),)),
        WorkflowStep("proteinmpnn", tuple(proteinmpnn_command) + (
            "--pdb_path", str(backbones_dir), "--out_folder", str(sequences_dir),
            "--num_seq_per_target", str(sequences_per_backbone),
        ), (str(backbones_dir),), (str(sequences_dir),)),
        WorkflowStep("rosetta", tuple(rosetta_command) + (
            "-parser:protocol", str(root / "refine.xml"), "-s", str(sequences_dir),
            "-out:path:all", str(refined_dir),
        ), (str(sequences_dir),), (str(refined_dir),)),
    ]


def write_workflow_manifest(path: str | Path, steps: Sequence[WorkflowStep]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(json.dumps(asdict(step)) for step in steps) + "\n", encoding="utf-8")


def run_workflow(steps: Sequence[WorkflowStep], *, dry_run: bool = True) -> list[dict[str, object]]:
    """Run a workflow serially, stopping at the first failed external program."""
    results: list[dict[str, object]] = []
    for step in steps:
        if dry_run:
            results.append({"name": step.name, "command": list(step.command), "status": "planned"})
            continue
        completed = subprocess.run(step.command, text=True, capture_output=True, check=False)
        result = {"name": step.name, "command": list(step.command), "returncode": completed.returncode,
                  "stdout": completed.stdout, "stderr": completed.stderr,
                  "status": "complete" if completed.returncode == 0 else "failed"}
        results.append(result)
        if completed.returncode != 0:
            raise RuntimeError(f"{step.name} failed with exit code {completed.returncode}: {completed.stderr.strip()}")
    return results


def _descriptor(sequence: str) -> tuple[float, ...]:
    core = sequence[1:-1]
    length = max(1, len(core))
    return tuple(sum(core.count(aa) for aa in residues) / length for residues in RESIDUE_GROUPS.values())


def _distance(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    return sum((a - b) ** 2 for a, b in zip(left, right)) ** 0.5


def build_disulfide_library(size: int, ring_length: int, seed: int = 17) -> list[dict[str, str]]:
    """Create a balanced Cys-Cys pilot design; outputs designs, not synthesized compounds."""
    if ring_length < 6 or ring_length > 18:
        raise ValueError("ring_length must be between 6 and 18 residues including terminal cysteines")
    if size < 1:
        raise ValueError("size must be positive")
    rng = random.Random(seed)
    alphabet = "ADEFGHIKLMNPQRSTVWY"  # internal Cys excluded to avoid ambiguous bridges
    pool: set[str] = set()
    while len(pool) < max(size * 40, 2000):
        core = "".join(rng.choice(alphabet) for _ in range(ring_length - 2))
        # Exclude strongly hydrophobic and highly charged sequences at pilot stage.
        hydrophobic = sum(core.count(x) for x in "AILMFWVY") / len(core)
        net_charge_proxy = abs(sum(core.count(x) for x in "KR") - sum(core.count(x) for x in "DE"))
        if 0.2 <= hydrophobic <= 0.65 and net_charge_proxy <= max(2, len(core) // 3):
            pool.add("C" + core + "C")
    candidates = sorted(pool)
    selected = [candidates[rng.randrange(len(candidates))]]
    candidate_desc = {seq: _descriptor(seq) for seq in candidates}
    remaining = set(candidates) - set(selected)
    min_distances = {
        seq: _distance(candidate_desc[seq], candidate_desc[selected[0]])
        for seq in remaining
    }
    while len(selected) < size:
        best = max(remaining, key=min_distances.__getitem__)
        selected.append(best)
        remaining.remove(best)
        min_distances.pop(best)
        for seq in remaining:
            min_distances[seq] = min(
                min_distances[seq], _distance(candidate_desc[seq], candidate_desc[best])
            )
    rows = []
    for index, sequence in enumerate(selected, 1):
        rows.append({
            "library_id": f"DPC-{index:05d}",
            "sequence": sequence,
            "topology": "cyclic",
            "cyclization": "Cys1-CysN disulfide",
            "n_terminal": "free",
            "c_terminal": "amide",
            "evidence": "generated",
            "design_method": "descriptor max-min diversity; internal Cys excluded",
        })
    return rows


def write_vendor_csv(path: str | Path, rows: list[dict[str, str]], target_id: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["library_id", "target_id", "sequence", "topology", "cyclization", "n_terminal", "c_terminal", "evidence", "design_method"]
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "target_id": target_id})
