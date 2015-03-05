all: sqi2fq tailseq-retrieve-signals
#	python setup.py install --user
#	python test.py

sqi2fq: sqi2fq.c
	cc -O3 -Wall -o $@ $<

tailseq-retrieve-signals: tailseq-retrieve-signals.c ssw.c ssw.h
	cc -O3 -Wall -o $@ tailseq-retrieve-signals.c ssw.c ssw.h -lz
