import sys
MERGE_TEMPLATE = 'samtools merge {outbam} {inbams} {merge_opts}'
INDEX_TEMPLATE = 'samtools index {inbam} {index_opts}'
CALL_TEMPLATE = {
	'freebayes': 'freebayes -f {genome} -r {region} {bam} {call_opts} > {outvcf}'
	}
def call_variants(ref, bams, out_dir, prefix=None, merge_opts='', index_opts='',
		bin_size=1000000, caller='freebayes', call_opts='', filt_opts='-f QUAL > 20'):
	if prefix is None:
		prefix = out_dir
	uncompleted = 0
	cmds = []
	faidx = ref + '.fai'
	if not os.path.exists(faidx):
		cmds += ['samtools faidx {}'.format(ref)]
	if len(bams) == 1:
		bam = bams[0]
		if not exist_bai(bam):
			cmds += [INDEX_TEMPLATE.format(bam=bam, index_opts=index_opts)]
	elif len(bams) > 1:
		bam = '{}/merged.bam'.format(out_dir)
		inbams = ' '.join(bams)
		cmd = MERGE_TEMPLATE.format(outbam=bam, inbams=inbams, merge_opts=merge_opts)
		cmd += ' && ' + INDEX_TEMPLATE.format(bam=bam, index_opts=index_opts)
		cmds += [cmd]
	cmd_file = '{}.precall.cmds'.format(out_dir)
	if cmds:
		uncompleted += run_job(cmd_file, cmds, **job_args)
	
	# call
	cmds = []
	vcfs = []
	i = 0
	for line in open(faidx):
		fields = line.strip().split("\t")
		chrom_name = fields[0]
		chrom_length = int(fields[1])
		for start in range(0, chrom_length, bin_size):
			i += 1
			end = start + bin_size
			if end > chrom_length:
				end = chrom_length
			region = '{chrom}:{start}-{end}'.format(chrom=chrom_name, start=start, end=end)
			outvcf = '{}/{}.{}.vcf'.format(out_dir, i, region)
			cmd = CALL_TEMPLATE[caller].format(genome=ref, region=region,
						bam=bam, call_opts=call_opts, outvcf=outvcf)
			cmds += [cmd]
			vcfs += [outvcf]
	cmd_file = '{}.call.cmds'.format(out_dir)
	job_args['cpu'] = 1
	job_args['mem'] = '2G'
	uncompleted += run_job(cmd_file, cmds, **job_args)
	
	# merge
	invcfs = ' '.join(vcfs)
	template = 'cat {invcfs} | vcffirstheader | vcffilter {filt} | bcftools sort - -O v -m 2G | vcfuniq | bgzip -c > {outvcf}'
	outvcf = '{}.vcf.gz'.format(prefix, )
	cmd = template.format(invcfs=invcfs, filt=filt_opts, outvcf=outvcf)
	cmd_file = '{}.vcf.cmds'.format(out_dir)
	cmds = [cmd]
	job_args['cpu'] = 4
	uncompleted += run_job(cmd_file, cmds, **job_args)
	if uncompleted == 0:
		shutil.rmtree(out_dir)
	return outvcf
def exist_bai(bam):
	if os.path.exists(bam+'.bai') or os.path.exists(bam[:-1]+'i'):
		return True
	else:
		return False
