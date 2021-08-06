
from bio_seq import bio_seq
from utilities import read_FASTA, readTextFile, writeTextFile

test_dna = bio_seq()
test_dna.generate_random_seq(50, "DNA")

print(test_dna.get_seq_info())
print(test_dna.get_seq_protein())

writeTextFile("test.txt", test_dna.seq)



