import utils
import argparse
from Bio import SeqIO

# reading wild genome file, variant genome file and set parameters from the terminal
parser = argparse.ArgumentParser()
parser.add_argument('variant', help='directory of the file with the variant genome')
parser.add_argument('wild', help='directory of the file with the wild genome')
parser.add_argument('-k', type=int, default = 50, help="length of the k-mers")
parser.add_argument('-f_min', type=int, default = 11, help="minimum frequency considered error free")
args = parser.parse_args()

# saving the reads of both genomes as two lists
print("Reading FASTA files")
records_v = list(SeqIO.parse(args.variant, "fasta"))
records_w = list(SeqIO.parse(args.wild, "fasta"))
print("Reading completed")

# saving the parameters k and f_min
k=args.k
f_min=args.f_min

# creating two dictionaries for both variant and wild genomes that have the k-mers as keys and the frequencies as values
counts_v=utils.count_all(records_v,k)
counts_w=utils.count_all(records_w,k)

# deleting from the dictionary the pairs which have values less than f_min
cut_counts_v=utils.clear_dict(counts_v,f_min)
cut_counts_w=utils.clear_dict(counts_w,f_min)

# removing the common k-mers that appear in both dictionaries
utils.delete_copies(cut_counts_v, cut_counts_w)

# concatenating k-mers and saving the new sequencies in a list
seq_v=utils.concatenate(cut_counts_v,k)
seq_w=utils.concatenate(cut_counts_w,k)

# printing the sequences containing the mutation
utils.print_results(seq_v,seq_w)