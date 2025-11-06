import os
import re
from pathlib import Path
import csv
import math
import matplotlib.pyplot as plt

FASTA_DIR = Path("fasta")
OUT_DIR = Path("output")
RECOGNITION = "GAATTC"  # EcoRI
PATTERN = re.compile(RECOGNITION, re.IGNORECASE)

OUT_DIR.mkdir(exist_ok=True)



def parse_fasta(path):
    name = None
    seq_parts = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip()
            else:
                seq_parts.append(line)
    seq = "".join(seq_parts).upper()
    return name if name else path.name, seq

def ecori_fragments_linear(seq):
    cuts = []
    for m in PATTERN.finditer(seq):
        cutpos = m.start() + 1  # cut between G (index m.start()) and A (m.start()+1)
        cuts.append(cutpos)
    cuts = sorted(set(cuts))
    if not cuts:
        return [len(seq)]
    frags = []
    frags.append(cuts[0])
    for i in range(1, len(cuts)):
        frags.append(cuts[i] - cuts[i-1])
    frags.append(len(seq) - cuts[-1])
    return [int(x) for x in frags if x > 0]

# collect fasta files
fasta_files = sorted(FASTA_DIR.glob("*.fasta"))
if not fasta_files:
    raise SystemExit("No .fasta files found in 'fasta/' directory. Place your 10 genome files there and re-run.")

results = []
all_fragments = {}

for fp in fasta_files:
    rec_name, seq = parse_fasta(fp)
    fragments = ecori_fragments_linear(seq)
    fragments_sorted = sorted(fragments, reverse=True)
    all_fragments[fp.name] = fragments_sorted
    results.append({
        "filename": fp.name,
        "record_name": rec_name,
        "genome_length_bp": len(seq),
        "num_ecori_sites": max(0, len(fragments)-1) if len(fragments)>0 else 0,
        "num_fragments": len(fragments_sorted),
        "largest_fragment_bp": max(fragments_sorted),
        "smallest_fragment_bp": min(fragments_sorted),
        "fragments": fragments_sorted
    })

# write summary CSV
summary_csv = OUT_DIR / "ecori_summary.csv"
with open(summary_csv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["filename","record_name","genome_length_bp","num_ecori_sites","num_fragments","largest_fragment_bp","smallest_fragment_bp"])
    for r in results:
        w.writerow([r["filename"], r["record_name"], r["genome_length_bp"], r["num_ecori_sites"], r["num_fragments"], r["largest_fragment_bp"], r["smallest_fragment_bp"]])

# write fragments CSV (one row per fragment)
frags_csv = OUT_DIR / "ecori_fragments.csv"
with open(frags_csv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["filename","fragment_index","fragment_bp"])
    for r in results:
        for i, f in enumerate(r["fragments"], start=1):
            w.writerow([r["filename"], i, f])

# combined gel plot
filenames = list(all_fragments.keys())
n = len(filenames)
plt.figure(figsize=(10, max(3, 0.4*n + 2)))
y_positions = list(range(n, 0, -1))
all_sizes = [s for fr in all_fragments.values() for s in fr]
min_size = min(all_sizes)
max_size = max(all_sizes)

for idx, fname in enumerate(filenames):
    frags = all_fragments[fname]
    ys = [y_positions[idx]] * len(frags)
    xs = frags
    plt.scatter(xs, ys, marker='|', s=300)
    plt.text(min_size * 0.9, y_positions[idx], fname, va='center', fontsize=8)

plt.xscale('log')
plt.xlabel("Fragment size (bp) [log scale]")
plt.yticks([])
plt.title("Simulated EcoRI digest gel (combined lanes)")
plt.tight_layout()
combined_path = OUT_DIR / "ecori_combined_gel.png"
plt.savefig(combined_path, dpi=200)
plt.close()

# individual plots
for fname, frags in all_fragments.items():
    plt.figure(figsize=(8,2))
    ys = list(range(len(frags), 0, -1))
    plt.scatter(frags, ys, marker='|', s=400)
    plt.xscale('log')
    plt.yticks([])
    plt.xlabel("Fragment size (bp) [log scale]")
    plt.title(f"EcoRI digest — {fname} — {len(frags)} fragments")
    plt.tight_layout()
    outpath = OUT_DIR / f"{fname}_ecori_gel.png"
    plt.savefig(outpath, dpi=200)
    plt.close()

# which genome has the most fragments
most = max(results, key=lambda r: r["num_fragments"])
print(f"{most['filename']} has the most fragments: {most['num_fragments']} fragments (genome length {most['genome_length_bp']} bp).")
print(f"Summary CSV: {summary_csv}")
print(f"Fragments CSV: {frags_csv}")
print(f"Combined gel image: {combined_path}")
print(f"Individual gel images saved to: {OUT_DIR} (one file per genome)")
