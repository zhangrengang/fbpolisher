#~/bin/env python
import os
import sys
import argparse
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s -%(levelname)s- %(message)s')
logger = logging.getLogger(__name__)
bindir = os.path.dirname(os.path.realpath(__file__))
sys.path = [bindir + '/bin'] + sys.path
os.environ['PATH'] = '{}:{}'.format(bindir, os.environ['PATH'])
from chunker import chunk_reads
from rotate_seqs import rotate_seqs
from mapper import map_reads
from caller import call_variants
from filter_vcf import filter_vcf
from vcf2fasta import vcf2fasta

__version__ = '1.0'
def Args():
	parser = argparse.ArgumentParser(version=__version__,
				formatter_class=argparse.RawDescriptionHelpFormatter,
				description='Polishing noisy genome assembly with short reads of high-quality.',
				epilog='''Examples of the reads list file:
>> for paired-end data:
lib_1.R1.fq.gz lib_1.R2.fq.gz
lib_2.R1.fq.gz lib_2.R2.fq.gz
>> for single-end data:
lib_1.R1.fq.gz
lib_2.R1.fq.gz'''
					)
	parser.add_argument("genome", action="store",type=str,
					help="genome to be polished in fasta format [required]")
	parser.add_argument("reads", action="store",type=str,
					help="list of Illumina reads data [required]")
	parser.add_argument("-iter", dest="iterations", action="store",
					default=3, type=int, metavar='INT',
					help="polish iterations [default=%(default)s]")
	parser.add_argument("-pre", dest="prefix", action="store",
					default=None, type=str, metavar='STR',
					help="output prefix [default=%(default)s]")
	parser.add_argument("-tmp", dest="tmp_dir", action="store",
					default='./tmp', type=str, metavar='DIR',
					help="directory for temporary files [default=%(default)s]")
	parser.add_argument("-circular", action="store_true",
					default=False,
					help="genome is circular [default=%(default)s]")

	group_mapp = parser.add_argument_group('map reads', 'options to map reads')
	group_mapp.add_argument("-mapper", action="store", metavar='CMD',
                    default='bwa mem', choices=['bwa mem', 'minimap2'],
                    help="mapper to map reads against genome [default=%(default)s]")
	group_mapp.add_argument("-map-opts", action="store",
                    default='-t 4', type=str, metavar='OPTS',
                    help="options for the mapper [default='%(default)s']")
	group_mapp.add_argument("-filt-opts", action="store",
                    default='-@ 1 -F 2304', type=str, metavar='OPTS',
                    help="options for `samtool view` [default='%(default)s']")
	group_mapp.add_argument("-sort-opts", action="store",
                    default='-@ 2', type=str, metavar='OPTS',
                    help="options for `samtools sort` [default='%(default)s']")
	group_mapp.add_argument("-merge-opts", action="store",
                    default='-@ 10', type=str, metavar='OPTS',
                    help="options for `samtools merge` [default='%(default)s']")
	group_mapp.add_argument("-markdup-opts", action="store",
                    default='-t 10 -r', type=str, metavar='OPTS',
                    help="options for `sambamba markdup` [default='%(default)s']")
	group_mapp.add_argument("-index-opts", action="store",
                    default='-@ 2', type=str, metavar='OPTS',
                    help="options for `samtools index` [default='%(default)s']")
	group_mapp.add_argument("-no-markdup", action="store_true",
                    default=False,
                    help="do not mark duplicates, which should be true for PCR-free library [default=%(default)s]")
	group_mapp.add_argument("-chunk", dest="chunk_size", action="store",
                    default=1000000, type=int, metavar='INT',
                    help="number of reads per chunk [default=%(default)s]")

	group_call = parser.add_argument_group('call cariants', 'options to call variants')
	group_call.add_argument("-caller", action="store", metavar='CMD',
                    default='freebayes', choices=['freebayes', ],
                    help="caller to call variants [default=%(default)s]")
	group_call.add_argument("-call-opts", action="store", metavar='OPTS',
                    default='--ploidy 2 --min-repeat-entropy 0 --min-alternate-count 2 --min-alternate-fraction 0.2 --use-best-n-alleles 4 --genotype-qualities', type=str,
                    help="options for the caller [default='%(default)s']")
	group_call.add_argument("-vcfilt-opts", action="store",
                    default='-f "QUAL > 20"', type=str, metavar='OPTS',
                    help="options for `vcffilter` [default='%(default)s']")
	group_call.add_argument("-bin", dest="bin_size",  action="store",
                    default=500000, type=int, metavar='INT',
                    help="bin size of genome for caller [default=%(default)s]")

	group_appl = parser.add_argument_group('apply variants', 'options to apply polish')
	group_appl.add_argument("-qual", action="store",
                    default=20, type=int, metavar='INT',
                    help="min QUAL of variants [default=%(default)s]")
	group_appl.add_argument("-dp", action="store",
                    default=10, type=int, metavar='INT',
                    help="min DP of variants [default=%(default)s]")
	group_task = parser.add_argument_group('task control', 'options to run tasks')
	group_task.add_argument('-m', "-mode", dest='mode', action="store",
                    default='local', type=str, choices=['local', 'grid'],
                    help="run mode [default=%(default)s]")
	group_task.add_argument("-tc", dest='tasks', metavar='INT',
                    action="store", default=50, type=int,
                    help="maximum number of simultaneously running tasks [default=%(default)s]")
	group_task.add_argument("-retry", action="store",
                    default=3, type=int, metavar='INT',
                    help="try times for tasks [default=%(default)s]")
	group_task.add_argument("-redo", action="store_true",
                    default=False,
                    help="re-run compeleted tasks [default=%(default)s]")
	group_task.add_argument("-grid-opts", action="store", type=str,
					default="-l h_vmem={mem} -pe mpi {cpu}", metavar='OPTS',
					help='grid options [default="%(default)s"]')
	group_task.add_argument("-cpu", action="store",
                    default=None, type=int, metavar='INT',
                    help="CPU required for each task [default=auto]")
	group_task.add_argument("-mem", action="store",
                    default=None, type=str, metavar='MEM',
                    help="Memory required for each task [default=auto]")

#	parser.print_help()
	args = parser.parse_args()
	if args.prefix is None:
		args.prefix = os.path.splitext(os.path.basename(args.genome))[0]
	makedirs(args.tmp_dir)
	return args

def pipeline(args):
	# glocal
	job_args = {
		'tc_tasks': args.tasks,
		'mode': args.mode,
		'grid_opts': args.grid_opts,
		'retry': args.retry,
		'cont': (not args.redo),
		'cpu': args.cpu,
		'mem': args.mem,
		}
	map_args = dict(
		mapper=args.mapper,
		map_opts=args.map_opts,
		filt_opts=args.filt_opts,
		sort_opts=args.sort_opts,
		merge_opts=args.merge_opts,
		rmdup_opts=args.markdup_opts,
		index_opts=args.index_opts,
		markdup=(not args.no_markdup),
		)
	call_args = dict(
		caller=args.caller,
		call_opts=args.call_opts,
		merge_opts=args.merge_opts,
		index_opts=args.index_opts,
		bin_size=args.bin_size,
		filt_opts=args.vcfilt_opts,
		)
	refbase = args.prefix
	# split reads
	out_dir = '{}/reads'.format(args.tmp_dir)
	makedirs(out_dir)
	out_dir = realpath(out_dir)
	logger.info('chunking reads into {}'.format(out_dir))
	reads_list = chunk_reads(reads_list=args.reads, out_dir=out_dir, chunk_size=args.chunk_size, **job_args)

	ref = args.genome #realpath(args.genome)
	for itr in range(1, args.iterations+1):
		logger.info('Polishing iteration-{}'.format(itr))
		# rotate
		if args.circular:
			rotated_genome = '{}.rotated.{}.fa'.format(refbase, itr)
			rotated_genome = realpath(rotated_genome)
			logger.info('rotating genome to {}'.format(rotated_genome))
			with open(rotated_genome, 'w') as fp:
				rotate_seqs(ref, args.iterations, fp)
			ref = rotated_genome
		# mapping
		out_dir = '{}/{}.mapping.{}'.format(args.tmp_dir, refbase, itr)
		makedirs(out_dir)
		out_dir = realpath(out_dir)
		logger.info('mapping to genome {}'.format(ref))
#		logger.info(map_args)
#		logger.info(job_args)
		map_args.update(job_args)
#		logger.info(map_args)
		bams = map_reads(ref, reads_list, out_dir=out_dir, **map_args)
		# calling
		out_dir = '{}/{}.calling.{}'.format(args.tmp_dir, refbase, itr)
		makedirs(out_dir)
		out_dir = realpath(out_dir)
		logger.info('calling from {}'.format(bams))
		call_args.update(job_args)
		vcf = call_variants(ref, bams, out_dir=out_dir, **call_args)

		# filter errors
		filt_vcf = '{}/{}.calling.{}.filt.vcf'.format(args.tmp_dir, refbase, itr)
		logger.info('filter out errors from {}'.format(vcf))
		with open(filt_vcf, 'w') as fp:
			filter_vcf(vcf, fp, minQUAL=args.qual, minDP=args.dp)
		# apply
		out_genome = '{}.polished.{}.fa'.format(refbase, itr)
		out_genome = realpath(out_genome)
		logger.info('applying variants {}'.format(filt_vcf))
		with open(out_genome, 'w') as fp:
			vcf2fasta(filt_vcf, ref, fp)
		ref = out_genome
	logger.info('Finished')
def exists(path):
	return os.path.exists(path)
def nonempty(path):
	return os.path.exists(path) and os.path.getsize(path) > 0
def realpath(path):
	return os.path.realpath(path)
def makedirs(path):
	if not exists(path):
		os.makedirs(path)
if __name__ == '__main__':
	pipeline(Args())
