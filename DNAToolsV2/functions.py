from bio_structs import *
from collections import Counter
import random


class bio_seq:
    """DNA sequence class. Defalt value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization and validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"The provided data does not seem to be correct {self.seq_type} sequence"

    # DNA FUNCTIONS
    def __validate(self):
        """Checking the strand for a valid DNA sequence."""
        return set(Nucleotides).issuperset(self.seq)

    def get_seq_biotype(self):
        """Return the sequence's type"""
        return self.seq_type

    def get_seq_info(self):
        """Return 4 strings with all sequence info"""
        return f"[Label]: {self.label}\n[Type]: {self.seq_type}\n[Sequence  ]: {self.seq}\n[Compliment]: " \
               f"{self.compliment()}\n[Nucloetide Frequency]: " \
               f"{self.nucleotide_frequency()}\n[Length]: " \
               f"{len(self.seq)}\n[AT Content]: {self.at_content()}\n[AT Content Subsec]: " \
               f"{self.at_content_subsec()}\n[GC Content]: {self.cg_content()}\n[GC Content Subsec]: " \
               f"{self.cg_content_subsec()}"

    def get_seq_protein(self):
        """Return the sequence Protein info"""
        rna_seq = self.transcription()
        return f"[Compliment RNA Strand]: {rna_seq}\n[Translation]: {self.translate_seq(rna_seq)}\n[Codon Usage]:" \
               f" {self.codon_usage('L')} "

    def generate_random_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence"""
        seq = ''.join([random.choice(Nucleotides) for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def nucleotide_frequency(self):
        """Count number of nucleotides in a sequence"""
        return dict(Counter(self.seq))

    def compliment(self):
        """Generate a complimentary DNA strand by swapping Adenine and Thymine and Guanine with Cytosine."""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        return self.seq.translate(mapping)

    def at_content(self):
        """AT content in a DNA/RNA sequence"""
        return round((self.seq.count("A") + self.seq.count("T")) / len(self.seq) * 100)

    def at_content_subsec(self, k=20):
        """AT content in a DNA or RNA sub-sequence length k. k=20(default)"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count("A") + subseq.count("T")) / len(subseq) * 100))
        return res

    def cg_content(self):
        """GC content in a DNA/RNA sequence"""
        return round((self.seq.count("C") + self.seq.count("G")) / len(self.seq) * 100)

    def cg_content_subsec(self, k=20):
        """CG content in a DNA or RNA sub-sequence length k. k=20(default)"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count("C") + subseq.count("G")) / len(subseq) * 100))
        return res

    # RNA
    def transcription(self):
        """Generate a compliment RNA strand (replace Thymine with Uracil)."""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "UAGC")
        else:
            return "This is not a DNA strand"
        self.rna_seq = self.seq.translate(mapping)
        return self.rna_seq

    # Protein
    def translate_seq(self, rna_seq, init_pos=0):
        """Translates a DNA sequence into an amino acid sequence"""
        return [RNA_Codons[rna_seq[pos:pos + 3]] for pos in range(init_pos, len(rna_seq) - 2, 3)]

    def codon_usage(self, amino_acid):
        """Provide the frequency of a codon"""
        tmpList = []
        for i in range(0, len(self.rna_seq) - 2, 3):
            if RNA_Codons[self.rna_seq[i:i + 3]] == amino_acid:
                tmpList.append(self.rna_seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for rna_seq in freqDict:
            freqDict[rna_seq] = round(freqDict[rna_seq] / totalWight, 2)
        return freqDict

    # def gen_reading_frames(self):
    #     """Generating the three reading frames of RNA sequences"""
    #     frames = [self.translate_seq(self.rna_seq, 0),
    #               self.translate_seq(self.rna_seq, 1),
    #               self.translate_seq(self.rna_seq, 2),
    #               ]
    #     seq = '\n'.join(map(str, frames))
    #     return seq
    # \n [Reading Frames]:\n{self.gen_reading_frames()}
