#!/bin/sh
echo "=== PhiX =========="
zcat testout/PhiX.sqi.gz | ./sqi2fq 1 51 | fastq_to_fasta -Q33 | \
bowtie2 -p 16 -x phix_control/phix_control --end-to-end --reorder -k 1 --no-unal -f --no-head /dev/stdin > /dev/null

echo "=== Unknown =========="
zcat testout/Unknown.sqi.gz | ./sqi2fq 1 51 | fastq_to_fasta -Q33 | \
bowtie2 -p 16 -x phix_control/phix_control --end-to-end --reorder -k 1 --no-unal -f --no-head /dev/stdin > /dev/null

