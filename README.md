# fbpolisher
It is developed to polish noisy genome assembly with high-quality short reads. 

### Installation ###
#### Dependencies:
+ [python 2.7](https://www.python.org/)  
	+ [dramaa](https://github.com/pygridtools/drmaa-python)
	+ [parallel python](https://www.parallelpython.com/): quickly install by `pip install pp`  

The binaries of following program are also provided in `bin/`. If them do not work, please re-install.
+ [bwa](https://github.com/lh3/bwa)
+ [samtools](https://github.com/samtools/samtools)
+ [sambamba](https://github.com/biod/sambamba)
+ [freebayes](https://github.com/ekg/freebayes)
+ [bcftools](https://github.com/samtools/bcftools)
+ [vcflib](https://github.com/vcflib/vcflib)

```
https://github.com/zhangrengang/fbpolisher
```

### Quick Start ###
```
git clone https://github.com/zhangrengang/fbpolisher
cd fbpolisher
cd test

# run
python ../fbpolisher.py ref.fa reads.list
```
#### Grid support ####

For SGE:
```
python ../fbpolisher.py ref.fa reads.list -m grid -grid-opts '-tc {tc} -l h_vmem={mem} -pe mpi {cpu}' -grid-array 'if [$SGE_TASK_ID -eq {id} ]; then\n{cmd}\nfi'
```
For Slurm:
```
python ../fbpolisher.py ref.fa reads.list -m grid -grid-opts '-c {cpu} -a %{tc} ' -grid-array 'if [$SLURM_ARRAY_TASK_ID -eq {id} ]; then\n{cmd}\nfi'
```

### Usage ###
```
$ python fbpolisher.py -h
usage: fbpolisher.py [-h] [-v] [-iter INT] [-pre STR] [-tmp DIR] [-circular]
                     [-mapper CMD] [-map-opts OPTS] [-filt-opts OPTS]
                     [-sort-opts OPTS] [-merge-opts OPTS] [-markdup-opts OPTS]
                     [-index-opts OPTS] [-no-markdup] [-chunk INT]
                     [-caller CMD] [-call-opts OPTS] [-vcfilt-opts OPTS]
                     [-bin INT] [-qual INT] [-dp INT] [-m {local,grid}]
                     [-tc INT] [-retry INT] [-redo] [-grid-opts OPTS]
                     [-grid-array OPTS] [-cpu INT] [-mem MEM]
                     genome reads

Polishing noisy genome assembly with short reads of high-quality.

positional arguments:
  genome                genome to be polished in fasta format [required]
  reads                 list of Illumina reads data [required]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -iter INT             polish iterations [default=3]
  -pre STR              output prefix [default=None]
  -tmp DIR              directory for temporary files [default=./tmp]
  -circular             genome is circular [default=False]

map reads:
  options to map reads

  -mapper CMD           mapper to map reads against genome [default=bwa mem]
  -map-opts OPTS        options for the mapper [default='-t 4']
  -filt-opts OPTS       options for `samtool view` [default='-@ 1 -F 2304']
  -sort-opts OPTS       options for `samtools sort` [default='-@ 2']
  -merge-opts OPTS      options for `samtools merge` [default='-@ 10']
  -markdup-opts OPTS    options for `sambamba markdup` [default='-t 10 -r']
  -index-opts OPTS      options for `samtools index` [default='-@ 2']
  -no-markdup           do not mark duplicates, which should be true for PCR-
                        free library [default=False]
  -chunk INT            number of reads per chunk [default=1000000]

call cariants:
  options to call variants

  -caller CMD           caller to call variants [default=freebayes]
  -call-opts OPTS       options for the caller [default='--ploidy 2 --min-
                        repeat-entropy 0 --min-alternate-count 2 --min-
                        alternate-fraction 0.2 --use-best-n-alleles 4
                        --genotype-qualities']
  -vcfilt-opts OPTS     options for `vcffilter` [default='-f "QUAL > 20"']
  -bin INT              bin size of genome for caller [default=500000]

apply variants:
  options to apply polish

  -qual INT             min QUAL of variants [default=20]
  -dp INT               min DP of variants [default=10]

task control:
  options to run tasks

  -m {local,grid}, -mode {local,grid}
                        run mode [default=local]
  -tc INT               maximum number of simultaneously running tasks
                        [default=100]
  -retry INT            try times for tasks [default=3]
  -redo                 re-run compeleted tasks [default=False]
  -grid-opts OPTS       grid options [default='-tc {tc} -l h_vmem={mem} -pe
                        mpi {cpu}']
  -grid-array OPTS      template to run job-array tasks [default='if [
                        $SGE_TASK_ID -eq {id} ]; then\n{cmd}\nfi']
  -cpu INT              CPU required for each task [default=auto]
  -mem MEM              Memory required for each task [default=auto]

Examples of the reads list file:
>> for paired-end data:
lib_1.R1.fq.gz lib_1.R2.fq.gz
lib_2.R1.fq.gz lib_2.R2.fq.gz
>> for single-end data:
lib_1.R1.fq.gz
lib_2.R1.fq.gz
```

