CLASSPATH=${CURDIR}/picard/picard-1.91.jar:${CURDIR}/picard/sam-1.91.jar:${CURDIR}

all: net/sf/picard/illumina/parser/TailSeqCIFDumper.class sqi2fq
	python setup.py install --user
	python test.py

	echo "#!/bin/sh" > cifdumper.sh
	echo "java -cp ${CLASSPATH} net.sf.picard.illumina.parser.TailSeqCIFDumper "'$$*' >> cifdumper.sh
	chmod 755 cifdumper.sh

net/sf/picard/illumina/parser/TailSeqCIFDumper.class: TailSeqCIFDumper.java
	mkdir -p net/sf/picard/illumina/parser
	javac -cp ${CLASSPATH} $< && mv -f TailSeqCIFDumper.class net/sf/picard/illumina/parser

sqi2fq: sqi2fq.c
	cc -O3 -Wall -o $@ $<
