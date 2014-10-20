ALLINDICES=phix.1.bt2 phix.2.bt2 phix.3.bt2 phix.4.bt2 phix.rev.1.bt2 phix.rev.2.bt2

all: ${ALLINDICES}

phix.1.bt2:
phix.2.bt2:
phix.3.bt2:
phix.4.bt2:
phix.rev.1.bt2:
phix.rev.2.bt2:
	bowtie2-build -o 4 --seed 0 phix.fasta phix

clean:
	rm -f ${ALLINDICES}
