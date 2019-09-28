#coding: utf-8
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from xopen import xopen as open

def vcf2dict(inVcf):
	d_vcf = {}
	for line in open(inVcf):
		if line.startswith('#'):
			continue
		temp = line.rstrip().split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = temp[:9]
		start = int(POS) - 1
		end = start + len(REF)
		QUAL = float(QUAL)
		if CHROM in d_vcf and (start, end) in d_vcf[CHROM]:
			if QUAL <= d_vcf[CHROM][(start, end)][2]:
				continue
			else:
				d_vcf[CHROM][(start, end)] = (REF, ALT, QUAL)
		else:
			try: d_vcf[CHROM][(start, end)] = (REF, ALT, QUAL)
			except KeyError: d_vcf[CHROM] = {(start, end): (REF, ALT, QUAL)}
	return rm_dup(d_vcf)
def rm_dup(d_vcf):
	for CHROM, d_var in d_vcf.items():
		rm_dup2(d_vcf, CHROM, d_var)
	return d_vcf
def rm_dup2(d_vcf, CHROM, d_var):
	dup = 0
	left = 0
	for (start, end), (REF, ALT, QUAL) in sorted(d_var.items()):
		if start < left:
			dup += 1
			print >>sys.stderr, 'Neighbor variants have overlap at %s %s, with %s. ' % (CHROM, [(start, end), (REF, ALT, QUAL)], last), 
			if QUAL > last[1][2]:
				try: del d_vcf[CHROM][last[0]]
				except KeyError: pass 	# 有连续的overlap，而上步已del
				print >>sys.stderr, 'Use %s.' % ([(start, end), (REF, ALT, QUAL)],)
			else:
				del d_vcf[CHROM][(start, end)]
				print >>sys.stderr, 'Use %s.' % (last, )
		left = end
		last = [(start, end), (REF, ALT, QUAL)]
	if dup == 0:
		return None
	else:
		d_var = d_vcf[CHROM]	# update
		return rm_dup2(d_vcf, CHROM, d_var)
def vcf2fasta(inVcf, inRef, outRef):
	d_vcf = vcf2dict(inVcf)
	if not isinstance(outRef, file):
		outChanges = outRef + '.changes'
		outRef = open(outRef, 'w')
		close = True
	else:
		outChanges = inRef + '.changes'
		close = False
	f = open(outChanges, 'w')
	for rc in SeqIO.parse(inRef, 'fasta'):
		polish = 1
		try: d_var = d_vcf[rc.id]
		except KeyError: polish = 0
		if not (polish and d_var):
			SeqIO.write(rc, outRef, 'fasta')
			continue
		segments = []
		left = 0
		cleft = 0
		for (start, end), (REF, ALT, _) in sorted(d_var.items()):
			if start < left:
				print >>sys.stderr, 'neighbor variants have overlap at %s %s %s, with %s' % (rc.id, (start, end), (REF, ALT), last)
				continue
			if set(list(REF.upper())) & set(['N']):  # fix bug
				left_base = str(rc.seq[start-1]).upper()
				right_base = str(rc.seq[end]).upper()
				if left_base == 'N' and right_base == "N":
					print >>sys.stderr, 'variant is within gap, skipping', rc.id, (start, end), (REF, ALT)
					continue
				elif left_base == 'N' or right_base == "N":	
					print >>sys.stderr, 'variant neighbor to gap, skipping', rc.id, (start, end), (REF, ALT)
					continue
##				raise ValueError('neighbor variants have overlap at %s %s %s' % (rc.id, (start, end), (REF, ALT)))
			segments += [str(rc.seq[left:start])]
			cstart = cleft + start - left
			segments += [ALT]
			cend = cstart + len(ALT)
			print >>f, '%s:%s-%s\t%s:%s-%s\t%s\t%s' % (rc.id, start, end, rc.id, cstart, cend, REF, ALT)
			left = end
			cleft = cend
			last = [(start, end), (REF, ALT)]
		segments += [str(rc.seq[left:])]
		rc.seq = Seq(''.join(segments))
		SeqIO.write(rc, outRef, 'fasta')
	if close:
		outRef.close()

if __name__ == '__main__':
	vcf2fasta(inVcf=sys.argv[1], inRef=sys.argv[2], outRef=sys.stdout)
			
			
