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

def transcribe(dna):
    return dna.upper().replace('T', 'U')

def translate(rna):
    start = rna.find('AUG')
    if start == -1:
        return "No start codon (AUG) found."

    aas = []
    for i in range(start, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) < 3:
            break
        aa = GENETIC_CODE.get(codon, '?')  # '?' = unknown codon
        if aa == 'STOP':
            break
        aas.append(aa)
    return "-".join(aas)

def dna_to_protein(dna):
    return translate(transcribe(dna))

if __name__ == "__main__":
    dna = "GGCATGTACCCGGATTGTCGTCAAGCCGCTAGCATAA"
    print("DNA:", dna)
    rna = transcribe(dna)
    print("mRNA:", rna)
    print("Protein:", translate(rna))