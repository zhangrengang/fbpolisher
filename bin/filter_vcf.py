import sys
from xopen import xopen as open
def filter_vcf(inVcf, outVcf, minQUAL=20, minDP=10):
	for line in open(inVcf):
		if line.startswith('#'):
			outVcf.write(line)
			continue
		temp = line.rstrip().split('\t')
		if len(temp) <= 9:
			continue
		sample = temp[9]
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = temp[:9]
		gt = sample.split(':')[0]
		gts = gt.split('/')
		if '0' in set(gts): # 0/0 or 0/1 etc.
			continue
		if float(QUAL) < minQUAL: # low QUAL
			continue
		dp = sample.split(':')[FORMAT.split(':').index('DP')]
		if int(dp) < minDP:	# low DP
			continue
		alt = gts[0]
		ALTs = ALT.split(',')
		##INFO = 'ALT=%s;SAMPLE=%s' % (ALT, sample)
		ALT = ALTs[int(alt)-1]
		##FORMAT = 'GT'
		##sample = '1'
		line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample]
		print >> outVcf, '\t'.join(line)
if __name__ == '__main__':
	filter_vcf(inVcf=sys.argv[1], outVcf=sys.stdout, minQUAL=20, minDP=10)
