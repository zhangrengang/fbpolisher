import os
import sys
import re
import shutil
from RunCmdsMP import run_job
SORT_TEMPLATE = 'samtools sort {sort_opts} > {outbam}'
FILT_TEMPLATE = 'samtools view -b {filt_opts}'
RG = r' -R "@RG\tID:{sample}\tSM:{sample}" '
MAP_TEMPLATE = {
	'bwa mem': 'bwa mem {ref} {r1} {r2} {map_opts} -M ' + RG,
	'minimap2': 'minimap2 -a -x sr {ref} {r1} {r2} {map_opts}' + RG,
	}
MERGE_TEMPLATE = 'samtools merge -pc -f {outbam} {inbams} {merge_opts}'
RMDUP_TEMPLATE = 'sambamba markdup {rmdup_opts} {inbams} {outbam}'
INDEX_TEMPLATE = 'samtools index {inbam} {index_opts}'
def map_reads(ref, reads_list, out_dir, prefix=None, 
			mapper='bwa mem', map_opts='', filt_opts='', sort_opts='',
			merge_opts='', markdup=True, rmdup_opts='',
			index_opts='',
			**job_args):
	if prefix is None:
		prefix = out_dir
	if job_args['mem'] is None:
		job_args['mem'] = pre_mem(ref)

	uncompleted = 0	
	# build index
	index_prefix = '{}/ref'.format(out_dir)
	cmds = build_index(ref, index_prefix, mapper)
	if cmds:
		cmd_file = '{}.index.cmds'.format(out_dir)
		uncompleted += run_job(cmd_file, cmds, **job_args)
		ref = index_prefix
	# map
	cmds = []
	d_bams = {}
	i = 0
	for line in open(reads_list):
		if line.startswith('#'):
			continue
		i += 1
		cols = line.strip().split('\t')
		if not cols:
			continue
		r1, r2, sample = cols[:3]
		outbam = '{}/{}.sorted.bam'.format(out_dir, i)
		cmd = MAP_TEMPLATE[mapper].format(
					ref=ref, r1=r1, r2=r2, map_opts=map_opts, 
					sample=sample, )
		if filt_opts:
			cmd += ' | ' + FILT_TEMPLATE.format(filt_opts=filt_opts)
		cmd += ' | ' + SORT_TEMPLATE.format(sort_opts=sort_opts,
						outbam=outbam)
		try: d_bams[sample] += [outbam]
		except KeyError: d_bams[sample] = [outbam]
		cmds += [cmd]
	cmd_file = '{}.map.cmds'.format(out_dir)
	
	job_args['cpu'] = max(job_args['cpu'], pre_cpu(cmd, func=sum))
	uncompleted += run_job(cmd_file, cmds, **job_args)
	
	# merge bams
	cmds = []
	bams = []
	for sample, inbams in d_bams.items():
		outbam = '{}.{}.tmp.bam'.format(prefix, sample)
		inbams = ' '.join(inbams)
		cmd = MERGE_TEMPLATE.format(outbam=outbam, inbams=inbams, merge_opts=merge_opts)
		finalbam = '{}.{}.bam'.format(prefix, sample)
		if markdup:
			cmd += ' && ' + RMDUP_TEMPLATE.format(outbam=finalbam, inbams=outbam, rmdup_opts=rmdup_opts) + \
				   ' && rm {tmpbam}'.format(tmpbam=outbam)
		else:
			cmd += ' && mv {} {}'.format(outbam, finalbam) + \
				   ' &&' +  INDEX_TEMPLATE.format(inbam=finalbam, index_opts=index_opts)
		cmds += [cmd]
		bams += [finalbam]
	cmd_file = '{}.merge.cmds'.format(out_dir)
	job_args['cpu'] = max(job_args['cpu'], pre_cpu(cmd, func=max))
	uncompleted += run_job(cmd_file, cmds, **job_args)
	if uncompleted == 0:
		shutil.rmtree(out_dir)
	else:
		raise ValueError('mapping failed')
	return bams

def pre_cpu(cmd, func=sum):
	ncpus = []
	try:
		cpus = re.compile(r'(-t|-@|--threads)\s(\d+)').findall(cmd)
		for par, cpu in cpus:
			ncpus += [int(cpu)]
	except AttributeError:
		ncpu += [1]
	return max(1, func(ncpus))
def pre_mem(ref):
	size = os.path.getsize(ref)
	return '{}G'.format(int(max(size/1e9,1)*2))
def build_index(ref, prefix, mapper):
	if mapper.startswith('bwa'):
#		indice = [ '{}.{}'.format(ref, suffix) for suffix in ['amb', 'ann', 'bwt', 'pac', 'sa']]
#		if all_noempty(indice):
#			pass
		cmd = 'bwa index {} -p {}'.format(ref, prefix)
		return [cmd]
	else:
		return []

def _test():
	ref, reads_list, out_dir = sys.argv[1:4]
	map_reads(ref, reads_list, out_dir)
if __name__ == '__main__':
	_test()
