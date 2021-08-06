
from functions import functions
from utilities import read_FASTA, readTextFile, writeTextFile

test_dna = functions()
test_dna.generate_random_seq(50, "DNA")

print(test_dna.get_seq_info())
print(test_dna.get_seq_protein())

writeTextFile("test.txt", test_dna.seq)



