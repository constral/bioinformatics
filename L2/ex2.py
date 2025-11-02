# find, in sequence S, only the dinucleotides and trinucleotides, w/o the use of the brute force engine.
# for this, one must verify this combination, starting from the beginning, until the end of the sequence.
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

# add up dinucleotides
dinucCounts = {}
for i in range(0, len(S)-1):
    k = S[i:i+2]
    # increment if already, initialize if not
    if k in dinucCounts:
        dinucCounts[k] = dinucCounts[k] + 1
    else:
        dinucCounts[k] = 1

# add up trinucleotides
trinucCounts = {}
for i in range(0, len(S)-1):
    k = S[i:i+3]
    if k in trinucCounts:
        trinucCounts[k] = trinucCounts[k] + 1
    else:
        trinucCounts[k] = 1
    

# after parsing everything, show total counts
print("drinucleotides")
for k in dinucCounts:
    print(k, str(dinucCounts[k] / len(S) * 100) + "%")

print("trinucleotides")
for k in trinucCounts:
    print(k, str(trinucCounts[k] / len(S) * 100) + "%")