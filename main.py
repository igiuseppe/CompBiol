import utils
import argparse
Import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('wild', help='directory of the file with the wild genome')
parser.add_argument('variant', help='directory of the file with the variant genome')
parser.add_argument("k", type=int, default = 60, help="length of the k-mers")
parser.add_argument("f_min", type=int, default = 11, help="minimum frequency considered error free")
args = parser.parse_args()

records_v = list(SeqIO.parse(args.variant, "fasta"))
records_w = list(SeqIO.parse(args.wild, "fasta"))

k=args.k
f_min=args.f_min


counts_v=utils.counts_all(records_v,k)
counts_w=utils.counts_all(records_w,k)

counts_v=utils.clear_dict(counts_v,f_min)
counts_w=utils.clear_dict(counts_w,f_min)

utils.delete_copies(counts_v, counts_w)

seq_v=utils.concatenate(counts_v,k)
seq_w=utils.concatenate(counts_w,k)

mat=utils.calc_dist(seq_v,seq_w)

ind=np.where(mat=mat.min())

utils.print_results(ind, mat)

