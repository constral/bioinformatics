from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt

GENETIC_CODE = {
    'UUU':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu',
    'UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser',
    'UAU':'Tyr','UAC':'Tyr','UAA':'STOP','UAG':'STOP',
    'UGU':'Cys','UGC':'Cys','UGA':'STOP','UGG':'Trp',
    'CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu',
    'CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'CAU':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AUU':'Ile','AUC':'Ile','AUA':'Ile','AUG':'Met',
    'ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
    'AAU':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'AGU':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
    'GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val',
    'GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'GAU':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'
}

AMINO_ACID_INFO = {
    'Leu': "Leucine: essential; high in meat/dairy/soy.",
    'Ser': "Serine: non-essential; common in many foods.",
    'Thr': "Threonine: essential; lower in some grains.",
    'Arg': "Arginine: conditionally-essential; nuts, seeds, meat.",
    'Lys': "Lysine: essential; abundant in legumes, meat, dairy; low in many grains.",
    'Glu': "Glutamate: non-essential; very common in foods.",
    'Ala': "Alanine: non-essential; abundant in animal proteins.",
    'Asn': "Asparagine: non-essential; found in dairy, meat, potatoes.",
    'Gln': "Glutamine: non-essential; common in many proteins.",
    'Val': "Valine: essential; found in high-protein foods and soy.",
    'Phe': "Phenylalanine: essential; abundant in protein-rich foods."
}

LOW_PROTEIN_FOODS = {
    "Fats & Oils": ["Olive oil, butter, margarine, mayonnaise"],
    "Sugars & Starches": ["Table sugar, honey, syrups, sorbet"],
    "Low-protein Fruits/Vegetables": ["Apples, grapes, cucumber, lettuce"],
    "Beverages": ["Water, coffee, tea (no milk)"]
}





# Read FASTA and return concatenated uppercase sequence
def parse_fasta(path):

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {path}")
    lines = []
    with p.open() as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            lines.append(line.strip())
    return "".join(lines).upper()
# DNA -> mRNA (T -> U), case-insensitive
def transcribe(dna):

    return dna.upper().replace('T', 'U')

# Count codons in frame 0
def get_codon_frequencies(rna):
    codons = (rna[i:i+3] for i in range(0, len(rna) - 2, 3))
    return Counter(c for c in codons if len(c) == 3)

# Map codon counts to amino-acid counts, skip STOP
def get_amino_acid_frequencies(codon_counts):

    aa_counts = Counter()
    for codon, cnt in codon_counts.items():
        aa = GENETIC_CODE.get(codon)
        if aa and aa != 'STOP':
            aa_counts[aa] += cnt
    return aa_counts

# Simple bar chart for top N codons
def plot_top_codons(codon_counts, title, top_n=10):

    top = codon_counts.most_common(top_n)
    if not top:
        return
    codons, counts = zip(*top)
    plt.figure(figsize=(10, 5))
    plt.bar(codons, counts)
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Count")
    plt.tight_layout()

# Print top amino acids
def print_top_amino_acids(aa_counts, label, top_n=3):

    print(f"\nTop {top_n} amino acids ({label}):")
    for i, (aa, cnt) in enumerate(aa_counts.most_common(top_n), 1):
        print(f"{i}. {aa}: {cnt}")

# Print short food guidance based on top amino acids
def find_and_print_food_recommendations(top_aa_set):

    top_list = sorted(top_aa_set)
    print("\nTop amino acids across analyses:", ", ".join(top_list) or "None")
    print("\nNotes:")
    for aa in top_list:
        print(f"- {AMINO_ACID_INFO.get(aa, aa + ': no info available')}")
    print("\nIf you want to avoid these, prefer low-protein categories like:")
    for cat, items in LOW_PROTEIN_FOODS.items():
        print(f"- {cat}: {', '.join(items)}")



def main():
    covid_file = "covid.fasta"
    flu_file = "influenza.fasta"

    # COVID
    covid_seq = parse_fasta(covid_file)
    covid_rna = transcribe(covid_seq)
    covid_codons = get_codon_frequencies(covid_rna)
    covid_aa = get_amino_acid_frequencies(covid_codons)
    plot_top_codons(covid_codons, "COVID-19: Top 10 Codons")

    # Influenza
    flu_seq = parse_fasta(flu_file)
    flu_rna = transcribe(flu_seq)
    flu_codons = get_codon_frequencies(flu_rna)
    flu_aa = get_amino_acid_frequencies(flu_codons)
    plot_top_codons(flu_codons, "Influenza: Top 10 Codons")

    print("\n" + "="*40)
    print("GENOME ANALYSIS SUMMARY")
    print("="*40)

    covid_top10 = {c for c, _ in covid_codons.most_common(10)}
    flu_top10 = {c for c, _ in flu_codons.most_common(10)}
    common = covid_top10 & flu_top10

    if common:
        print("\nCommon codons in top 10 (codon: covid_count, flu_count):")
        for codon in sorted(common, key=lambda c: covid_codons[c], reverse=True):
            print(f"{codon}: {covid_codons[codon]}, {flu_codons[codon]}")
    else:
        print("\nNo common codons found in the top 10 lists.")

    print_top_amino_acids(covid_aa, "COVID-19")
    print_top_amino_acids(flu_aa, "Influenza")

    # Food suggestions based on top 3 amino acids from each genome
    top_aa = {aa for aa, _ in covid_aa.most_common(3)} | {aa for aa, _ in flu_aa.most_common(3)}
    find_and_print_food_recommendations(top_aa)

    # Show plots
    print("\nDisplaying charts (close windows to finish).")
    plt.show()

if __name__ == "__main__":
    main()
