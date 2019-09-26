import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
def rotate_seqs(inSeq, inCut, outSeq, gap_length=100):
	for rc in SeqIO.parse(inSeq,'fasta'):
		cut_index = len(rc.seq)/inCut
		print >>sys.stderr, rc.id,cut_index
		if not re.compile(r'circular=no', re.I).search(rc.description):
			rc.seq = rc.seq[cut_index:] + rc.seq[:cut_index]
		else:
			rc.seq = rc.seq #[cut_index:] + Seq('N'*gap_length) + rc.seq[:cut_index]
#			print >>sys.stderr, rc.id, "not circular, add N"
		SeqIO.write(rc, outSeq, 'fasta')

if __name__ == '__main__':
	rotate_seqs(inSeq=sys.argv[1], inCut=int(sys.argv[2]), gap_length=100, outSeq=sys.stdout)
