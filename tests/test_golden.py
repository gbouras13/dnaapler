"""
Golden-output regression tests.

Each test runs a dnaapler command and compares its output files against a stored
snapshot in tests/test_data/golden/. The snapshot records, for each output file, a
compact but exact fingerprint - reoriented sequences are captured as lengths + sha256
hashes (so large genomes don't bloat the repo), and small text files (e.g. the
no-reorientation summary) are stored verbatim. Any change to dnaapler's reoriented
output is therefore caught.

Only the primary reoriented outputs (the FASTA / GFA, and the fully deterministic
no-circular summary) are goldened. Summary *statistics* for real reorientations
(coverage %, identity %, etc.) are intentionally not goldened, as those can shift with
MMseqs2 version changes without indicating a dnaapler regression.

If an output change is intentional (a new feature, or a dependency bump that
legitimately changes reorientation), regenerate the snapshots with:

    DNAAPLER_UPDATE_GOLDEN=1 pytest tests/test_golden.py

and commit the updated tests/test_data/golden/*.json files.
"""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
from Bio import SeqIO

test_data = Path("tests/test_data")
overall_test_data = test_data / "overall_inputs"
golden_dir = test_data / "golden"

UPDATE = os.environ.get("DNAAPLER_UPDATE_GOLDEN") == "1"


def exec_command(cmnd):
    """Executes a shell command, raising if it exits non-zero."""
    proc = subprocess.Popen(
        cmnd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err.decode('utf8', 'replace')}")
    return out.decode("utf8") if out is not None else None


def _sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


def snapshot_fasta(path):
    """Snapshot a FASTA as an ordered list of records (id, description, length, hash)."""
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq)
        records.append(
            {
                "id": rec.id,
                "description": rec.description,
                "length": len(seq),
                "sha256": _sha256(seq),
                "head": seq[:60],
            }
        )
    return {"type": "fasta", "records": records}


def snapshot_gfa(path):
    """Snapshot a GFA: S-line sequences as hashes, all other lines verbatim."""
    non_s_lines = []
    s_lines = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "S":
                seq = parts[2]
                s_lines[parts[1]] = {
                    "length": len(seq),
                    "sha256": _sha256(seq),
                    "tags": parts[3:],
                }
            else:
                non_s_lines.append(line.rstrip("\n"))
    return {"type": "gfa", "non_s_lines": non_s_lines, "s_lines": s_lines}


def snapshot_text(path):
    """Snapshot a small text file verbatim."""
    with open(path) as fh:
        return {"type": "text", "content": fh.read()}


def snapshot_file(path):
    suffix = Path(path).suffix
    if suffix == ".fasta":
        return snapshot_fasta(path)
    if suffix == ".gfa":
        return snapshot_gfa(path)
    return snapshot_text(path)


def snapshot_outputs(output_dir, filenames):
    return {name: snapshot_file(os.path.join(output_dir, name)) for name in filenames}


def check_golden(name, snapshot):
    """Compare snapshot to the stored golden, or (re)write it when updating."""
    golden_path = golden_dir / f"{name}.json"
    if UPDATE:
        golden_dir.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
        return
    assert golden_path.exists(), (
        f"Missing golden file {golden_path}. "
        "Regenerate with: DNAAPLER_UPDATE_GOLDEN=1 pytest tests/test_golden.py"
    )
    expected = json.loads(golden_path.read_text())
    assert snapshot == expected, (
        f"Output for '{name}' differs from golden {golden_path}. "
        "If this change is intentional, regenerate with: "
        "DNAAPLER_UPDATE_GOLDEN=1 pytest tests/test_golden.py"
    )


# (name, command template, [output files to snapshot])
# {overall}, {test_data} and {out} are filled in below.
GOLDEN_CASES = [
    (
        "chromosome_dnaa",
        "dnaapler chromosome -i {overall}/chromosome.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        "plasmid_repa",
        "dnaapler plasmid -i {overall}/plasmid.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        "phage_terl",
        "dnaapler phage -i {overall}/NC_007458.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        # tophit alignment has no valid start codon -> overlapping-ORF reorientation path
        "phage_terl_no_start_codon",
        "dnaapler phage -i {overall}/SAOMS1.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        "archaea_cog1474",
        "dnaapler archaea -i {overall}/CP001742.1_archaea.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        # largest CDS is on the negative strand - guards issue #102
        "largest_neg_strand",
        "dnaapler largest -i {test_data}/NC_007458_rc.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        "nearest",
        "dnaapler nearest -i {overall}/chromosome.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        # mystery uses the default random seed (13) and is therefore deterministic
        "mystery_seed13",
        "dnaapler mystery -i {overall}/chromosome.fasta -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta"],
    ),
    (
        # GFA input: reorients circular contigs and writes both GFA and complete FASTA
        "all_gfa",
        "dnaapler all -i {overall}/all_test.gfa -o {out} -t 1 -f",
        ["dnaapler_reoriented.fasta", "dnaapler_reoriented.gfa"],
    ),
    (
        # GFA with no circular sequences: passthrough GFA + linear FASTA + summary
        "all_gfa_no_circular",
        "dnaapler all -i {overall}/no_circular.gfa -o {out} -t 1 -f",
        [
            "dnaapler_reoriented.fasta",
            "dnaapler_reoriented.gfa",
            "dnaapler_all_reorientation_summary.tsv",
        ],
    ),
]


@pytest.mark.parametrize(
    "name,cmd_tmpl,outputs", GOLDEN_CASES, ids=[c[0] for c in GOLDEN_CASES]
)
def test_golden(tmp_path, name, cmd_tmpl, outputs):
    out = str(tmp_path / name)
    cmd = cmd_tmpl.format(overall=overall_test_data, test_data=test_data, out=out)
    exec_command(cmd)
    snapshot = snapshot_outputs(out, outputs)
    check_golden(name, snapshot)
