# The melting temperature (Tm) is the temperature at which one-half of a particular DNA duplex will dissociate
# and become a single strand of DNA. Primer length and sequence are of critical importance in designing the parameters
# of a successful amplification. The melting temperature of a nucleic acid duplex increases both with its length,
# and with increasing GC content.



import math

seq = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
# ion concentration of the solution
Na = 50.0


length = len(seq)
a = seq.count("A")
t = seq.count("T")
g = seq.count("G")
c = seq.count("C")



def simple_formula():
    #  A simple formula for calculation of the (Tm) is:  Tm = 4(G + C) + 2(A + T) °C
    result =  4*(g+c) + 2*(a+t)
    return result

def actual_formula(Na):
    
    # Na can take a value of 0.001
    na_M = Na / 1000.0

    # %GC
    percent_gc = ((g+c) / length) * 100.0

    # An alternative formula is: Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
    result = 81.5 + 16.6 * math.log10(na_M)
    result += 0.41 * percent_gc - 600.0 / length
    return result



print("simple: melting temp = " + str(simple_formula()) + "C")
print("actual formula: melting temp = " + str(actual_formula(Na)) + "C")