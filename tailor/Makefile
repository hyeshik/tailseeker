CC=	gcc
#CFLAGS=	-O3 -Wall -Werror -pg
CFLAGS=	-O3 -Wall -Werror

PROG=	bin/sqi2fq bin/tailseq-retrieve-signals

RETRIEVE_SIGNALS_OBJECTS=	\
	src/retrieve-signals/tailseq-retrieve-signals.o \
	src/retrieve-signals/cifreader.o \
	src/retrieve-signals/bclreader.o \
	src/retrieve-signals/altcalls.o \
	src/retrieve-signals/phix_control.o \
	src/retrieve-signals/controlaligner.o \
	src/retrieve-signals/my_strstr.o \
	src/retrieve-signals/ssw.o

SQI2FQ_OBJECTS=	\
	src/sqi2fq/sqi2fq.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${DFLAGS} ${INCLUDES} $< -o $@

all: bin/sqi2fq bin/tailseq-retrieve-signals
#	python setup.py install --user
#	python test.py

clean:
	rm -f ${PROG} ${RETRIEVE_SIGNALS_OBJECTS} ${SQI2FQ_OBJECTS}

bin/sqi2fq: ${SQI2FQ_OBJECTS}
	mkdir -p bin
	${CC} ${CFLAGS} -o $@ $<

bin/tailseq-retrieve-signals: ${RETRIEVE_SIGNALS_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${RETRIEVE_SIGNALS_OBJECTS} -lz

