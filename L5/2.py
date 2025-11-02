import os
import random
import time
import matplotlib.pyplot as plt

# parameters
FASTA_DIR = "fasta"
OUT_DIR = "outputs"
N_READS = 2000
MIN_READ = 100
MAX_READ = 150
MIN_OVERLAP = 10
PLOT_FILE = "gc_vs_time.png"

random.seed(17)





def read_fasta(path):
    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()

# sample n_reads random substrings from genome
def sample_reads(genome, n_reads):
    L = len(genome)
    reads = []
    for _ in range(n_reads):
        rlen = random.randint(MIN_READ, MAX_READ)
        if rlen >= L:
            start = 0
        else:
            start = random.randint(0, L - rlen)
        reads.append(genome[start:start + rlen])
    return reads

# longest exact suffix(a) == prefix(b) with minimum threshold
def best_overlap(a, b, min_ov=MIN_OVERLAP):
    m = min(len(a), len(b))
    for l in range(m, min_ov - 1, -1):
        if a[-l:] == b[:l]:
            return l
    return 0

# naive greedy assembler (exact matches only)
def greedy_assemble(reads):
    contigs = []
    reads = reads[:]  # work on a copy
    while reads:
        contig = reads.pop()
        extended = True
        while extended:
            extended = False
            best_i = -1
            best_l = 0
            best_dir = 0
            for i, r in enumerate(reads):
                ol = best_overlap(contig, r)
                if ol > best_l:
                    best_l, best_i, best_dir = ol, i, 1
                ol = best_overlap(r, contig)
                if ol > best_l:
                    best_l, best_i, best_dir = ol, i, -1
            if best_l >= MIN_OVERLAP:
                r = reads.pop(best_i)
                if best_dir == 1:
                    contig = contig + r[best_l:]
                else:
                    contig = r + contig[best_l:]
                extended = True
        contigs.append(contig)
    return contigs

# compute gc percent of a sequence
def gc_percent(seq):
    if not seq:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return 100.0 * (g + c) / len(seq)

def main():
    if not os.path.isdir(FASTA_DIR):
        raise SystemExit(f"fasta directory not found: {FASTA_DIR}")

    os.makedirs(OUT_DIR, exist_ok=True)

    results = []  # list of tuples: (label, gc%, time_ms, n_contigs, largest_contig)
    fasta_files = sorted(
        [f for f in os.listdir(FASTA_DIR) if os.path.isfile(os.path.join(FASTA_DIR, f))]
    )

    for fn in fasta_files:
        path = os.path.join(FASTA_DIR, fn)
        label = os.path.splitext(fn)[0]
        seq = read_fasta(path)

        # sample reads
        reads = sample_reads(seq, N_READS)

        # time the naive assembly
        t0 = time.perf_counter()
        contigs = greedy_assemble(reads)
        t1 = time.perf_counter()
        time_ms = (t1 - t0) * 1000.0

        contigs.sort(key=len, reverse=True)
        n_contigs = len(contigs)
        largest = len(contigs[0]) if contigs else 0
        gc = gc_percent(seq)

        # save small per-genome outputs for inspection
        with open(os.path.join(OUT_DIR, f"{label}_reads.fasta"), "w") as fh:
            for i, r in enumerate(reads[:1000]):  # keep reads file small
                fh.write(f">r{i}\n{r}\n")
        with open(os.path.join(OUT_DIR, f"{label}_contigs.fasta"), "w") as fh:
            for i, c in enumerate(contigs[:1000]):
                fh.write(f">c{i}\n{c}\n")

        results.append((label, gc, time_ms, n_contigs, largest))

    # create scatter plot: x = gc%, y = time_ms
    xs = [r[1] for r in results]
    ys = [r[2] for r in results]
    labels = [r[0] for r in results]

    plt.figure(figsize=(9, 6))
    plt.scatter(xs, ys)
    for x, y, lab in zip(xs, ys, labels):
        plt.annotate(lab, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)
    plt.xlabel("C+G percentage")
    plt.ylabel("Assembly time (ms)")
    plt.title("Assembly time vs C+G percentage")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=150)



if __name__ == "__main__":
    main()