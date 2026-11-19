# make_inputs.py
import json, re, sys
from pathlib import Path

fasta_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "fastas")

pips = {  # same acyl chains for all
    "PI3P": "CCCC(=O)OC[C@@H](COP(=O)(O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)OP(=O)(O)O)O)O)O)OC(=O)CCC",
    "PI4P": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC1[C@@H]([C@H](C([C@H]([C@H]1O)O)OP(=O)(O)O)O)O)OC(=O)CCCCCCCCCCCCCCC",
    "PI45P2": "CCCCCCCC(=O)OC[C@@H](COP(=O)(O)O[C@H]1[C@H]([C@H]([C@@H]([C@H]([C@H]1O)OP(=O)(O)O)OP(=O)(O)O)O)O)OC(=O)CCCCCCC",
    "ctrl_PI": "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC1[C@@H]([C@H](C([C@H]([C@H]1O)O)O)O)O)OC(=O)CCCCCCCCCCCCCCCCC",
}

def read_fasta(path):
    name, seq = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = re.sub(r"[^\w.-]", "_", line[1:].split()[0]), []
        elif line:
            seq.append(line.upper())
    if name:
        yield name, "".join(seq)

domains = {}
for f in sorted(fasta_dir.glob("*")):
    if f.suffix.lower() in {".fa", ".fasta", ".faa"}:
        for name, seq in read_fasta(f):
            if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
                print(f"WARNING: {name} has non-standard characters, skipping")
                continue
            domains[name] = seq

jobs = [
    {
        "name": f"{d}__{p}",
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": 1}},
            {"ligand": {"ligand": smi, "count": 1}},
        ],
    }
    for d, seq in domains.items()
    for p, smi in pips.items()
]

json.dump(jobs, open("batch.json", "w"), indent=2)
print(len(domains), "domains x", len(pips), "ligands =", len(jobs), "jobs")
