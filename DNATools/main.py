# Testing file

from DNATools import *
from utilities import colored
import random

# Generating a random DNA string for testing
randDNAStr = "".join([random.choice(Nucleotides) for nuc in range(50)])
testStrand = ""

# Print Results
DNAstr = validateSeq(randDNAStr)
print(f"1) DNA Sequence: {colored(DNAstr)}\n")
print(f"2) DNA Sequence length: {len(DNAstr)}\n")
print( colored(f"3. Nucleotide Frequency: {countNucFrequency(DNAstr)}\n"))
print(f"4) Transcription(RNA): {colored(transcription(DNAstr))}\n")
print(f"5) Compliment Stand: {colored(compliment(DNAstr))}\n")
print(f"6) DNA Strand and Compliment:\n5' {colored(DNAstr)} 3'")
print(f"   {''.join(['|' for c in range(len(DNAstr))])}")
print(f"3' {colored(compliment(DNAstr))} 5'\n")
print(f"7) AT content: {at_content(DNAstr)}%\n")
print(f"8) AT content in subsection k=5: {at_content_subsec(DNAstr, k=5)}\n")
print(f"9) CG content: {cg_content(DNAstr)}%\n")
print(f"10) CG content in subsection k=5: {cg_content_subsec(DNAstr, k=5)}\n")
print(f"11) Amino acids Sequence of DNA: {translate_seq(DNAstr, 0)}\n")
print("12) Reading Frames: ")
for frame in gen_reading_frames(DNAstr):
    print(frame)