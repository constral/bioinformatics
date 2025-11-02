# Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence:
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

# 1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.
# 2. For each combination, find out the percentage inside the S sequence.
# 3. Show the percentage for each combination in the output of your implementation.

# for all dinucleotides
for a in ["A", "C", "G", "T"]:
    for b in ["A", "C", "G", "T"]:
        dinuc = a + b
        print(dinuc)
        cnt = 0
        i = 0
        
        # parse the sequence
        for i in range(0, len(S)-1):
            if S[i:i+2] == dinuc:
                cnt = cnt + 1
        print(str(cnt / len(S) * 100) + "%")
        
        

# for all dinucleotides
for a in ["A", "C", "G", "T"]:
    for b in ["A", "C", "G", "T"]:
        for c in ["A", "C", "G", "T"]:
            trinuc = a + b + c
            print(trinuc)
            cnt = 0
            i = 0
            # parse the sequence
            for i in range(0, len(S)-1):
                if S[i:i+3] == trinuc:
                    cnt = cnt + 1
            print(str(cnt / len(S) * 100) + "%")
        
            



