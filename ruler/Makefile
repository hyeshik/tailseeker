CC=	gcc

RETRIEVE_SIGNALS_SOURCES=	\
	src/retrieve-signals/tailseq-retrieve-signals.c \
	src/retrieve-signals/cifreader.c \
	src/retrieve-signals/bclreader.c \
	src/retrieve-signals/phix_control.c \
	src/retrieve-signals/ssw.c

all: bin/sqi2fq bin/tailseq-retrieve-signals
#	python setup.py install --user
#	python test.py

bin/sqi2fq: src/sqi2fq/sqi2fq.c
	mkdir -p bin
	${CC} -O3 -Wall -o $@ $<

bin/tailseq-retrieve-signals: ${RETRIEVE_SIGNALS_SOURCES}
	${CC} -O3 -Wall -o $@ ${RETRIEVE_SIGNALS_SOURCES} -lz
