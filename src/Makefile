CC=	gcc
#CFLAGS=	-O3 -Wall -Werror -pg
CFLAGS=	-O3 -Wall -Werror
bindir=	../bin

PROG=	${bindir}/sqi2fq ${bindir}/tailseq-retrieve-signals

RETRIEVE_SIGNALS_OBJECTS=	\
	retrieve-signals/tailseq-retrieve-signals.o \
	retrieve-signals/cifreader.o \
	retrieve-signals/bclreader.o \
	retrieve-signals/altcalls.o \
	retrieve-signals/phix_control.o \
	retrieve-signals/controlaligner.o \
	retrieve-signals/my_strstr.o \
	retrieve-signals/ssw.o

SQI2FQ_OBJECTS=	\
	sqi2fq/sqi2fq.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${DFLAGS} ${INCLUDES} $< -o $@

all: ${bindir}/sqi2fq ${bindir}/tailseq-retrieve-signals
#	python setup.py install --user
#	python test.py

clean:
	rm -f ${RETRIEVE_SIGNALS_OBJECTS} ${SQI2FQ_OBJECTS}

distclean: clean
	rm -f ${PROG}

${bindir}/sqi2fq: ${SQI2FQ_OBJECTS}
	mkdir -p ${bindir}
	${CC} ${CFLAGS} -o $@ $<

${bindir}/tailseq-retrieve-signals: ${RETRIEVE_SIGNALS_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${RETRIEVE_SIGNALS_OBJECTS} -lz

