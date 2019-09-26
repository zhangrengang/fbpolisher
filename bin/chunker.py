import sys
import os
from RunCmdsMP import run_job
bindir = os.path.dirname(os.path.realpath(__file__))

def chunk_reads(reads_list, out_dir, ckpt=None, chunk_size=500000, seqfmt='fastq', **job_args):
	template = '{cat} {file} | python {dir}/split_records.py --chunk-size {size} --prefix {prefix} --format {format} -pfn > {filenames}'
	pre_tplt = '{dir}/{id}-R{end}'
	if ckpt is None:
		ckpt = '{}.list'.format(out_dir)
	reads1, reads2, samples = parse_list(reads_list)
	i = 0
	cmds = []
	for r1, r2, sample in zip(reads1, reads2, samples):
		i += 1
		for read, end in zip([r1, r2], [1, 2]):
			if not exists(read):
				continue
			cat = xcat(read)
			prefix = pre_tplt.format(dir=out_dir, id=i, end=end)
			split_list = prefix + '.list'
			cmd = template.format(
				cat=cat, file=read,
				dir=bindir, size=chunk_size,
				prefix=prefix, format=seqfmt,
				filenames=split_list
				)
			cmds.append(cmd)
	cmd_file = '{}.cmds'.format(out_dir)
	uncompleted = run_job(cmd_file, cmds, **job_args)

	# output list
	i = 0
	fp = open(ckpt, 'w')
	for r1, r2, sample in zip(reads1, reads2, samples):
		i += 1
		r1_list, r2_list = [pre_tplt.format(dir=out_dir, id=i, end=end) + '.list' \
							for read, end in zip([r1, r2], [1, 2])]
		r1_files = [line.strip() for line in open(r1_list)]
		if exists(r2_list):
			r2_files = [line.strip() for line in open(r2_list)]
		else:
			r2_files = [''] * len(r1_files)
		assert len(r1_file) == len(r2_files)
		for sr1, sr2 in zip(r1_files, r2_files):
			print >>fp, '\t'.join([sr1,sr2,sample])
	fp.close()
	return ckpt
def parse_list(reads_list):
	'''input file in mutiple format:
1-column: R1
2-column: R1 R2
2-column: R1 sample
3-column: R1 R2 sample
'''
	reads1, reads2, samples = [], [], []
	for line in open(reads_list):
		cols = line.strip().split()
		if not cols:
			continue
		if len(cols) >= 3:
			r1, r2, sample = cols[:3]
			assert exists(r1) and exists(r2)
			reads1.append(r1)
			reads2.append(r2)
			samples.append(sample)
		elif len(cols) == 1:
			r1, r2, sample = cols[0], None, None
			assert exists(r1)
		elif len(cols) == 2:
			r1, r2, sample = cols[0], cols[1], None
			if not exists(r2):
				r2, sample = None, r2
		sample = os.path.basename(str(sample))
		reads1.append(r1)
		reads2.append(r2)
		samples.append(sample)
	return reads1, reads2, samples
def exists(path):
	return os.path.exists(path)
def xcat(input_file):
	dcat = {
		'gz': 'zcat',
		'bz2': 'bzcat',
		'xz': 'xzcat',
		'bam': 'samtools fastq',
	}
	suffix = os.path.splitext(input_file)[-1].strip('.')
	try:
		return dcat[suffix]
	except KeyError:
		return 'cat'
