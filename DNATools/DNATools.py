
#Tools
import collections
from structures import *

#Cheking if valid DNA
def validateSeq(dna_seq):
    """Checking the strand for a valid DNA sequence."""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
        return tmpseq

#Coutning nucleotides
def countNucFrequency(seq):
    """Counting the number of nucleotides."""
    #tmpFreqDict = {"A": 0, "T": 0, "C": 0, "G": 0}
    #for nuc in seq:
    #    tmpFreqDict[nuc] += 1
    #return tmpFreqDict
    return dict(collections.Counter(seq))

#Transcription (DNA->RNA)
def transcription(seq):
    """Transcription from DNA to RNA by replacing Thymine with Uracil."""
    return seq.replace("T", "U")

#Compliment DNA strand
def reverse_compliment(seq):
    """Generate a complimentary DNA strand by swapping Adenine and Thymine and Guanine with Cytosine."""
    # return  "".join([DNA_Complement[nuc] for nuc in seq])
    mapping = str.maketrans("ATCG", "TAGC")
    return seq.translate(mapping)

#AT content in a DNA or RNA sequence
def at_content(seq):
    """AT content in a DNA/RNA sequence"""
    return round((seq.count("A") + seq.count("T")) / len(seq) * 100)

#AT content sub sequence
def at_content_subsec(seq, k=20):
    """AT content in a DNA or RNA sub-sequence length k. k=5(default)"""
    res = []
    for i in range(0, len(seq) - k +1, k):
        subseq = seq[i:i +k]
        res.append(at_content(subseq))
    return res

#CG content in a DNA or RNA sequence
def cg_content(seq):
    """GC content in a DNA/RNA sequence"""
    return round((seq.count("C") + seq.count("G")) / len(seq) * 100)

#CG content sub sequence
def cg_content_subsec(seq, k=20):
    """CG content in a DNA or RNA sub-sequence length k. k=5(default)"""
    res = []
    for i in range(0, len(seq) - k +1, k):
        subseq = seq[i:i +k]
        res.append(cg_content(subseq))
    return res

#Translation of DNA into amino acid
def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an amino acid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]



