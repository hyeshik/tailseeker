CC=	gcc

all: bin/sqi2fq bin/tailseq-retrieve-signals
#	python setup.py install --user
#	python test.py

bin/sqi2fq: src/sqi2fq/sqi2fq.c
	mkdir -p bin
	${CC} -O3 -Wall -o $@ $<

bin/tailseq-retrieve-signals: src/retrieve-signals/tailseq-retrieve-signals.c src/retrieve-signals/ssw.c src/retrieve-signals/ssw.h
	${CC} -O3 -Wall -o $@ src/retrieve-signals/tailseq-retrieve-signals.c src/retrieve-signals/ssw.c src/retrieve-signals/ssw.h -lz
