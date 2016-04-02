CC=	gcc
#CFLAGS=	-O3 -Wall -Werror -pg
CFLAGS=	-O3 -Wall -Werror
LIBS=	-lz -lm
bindir=	../bin

PROG=	${bindir}/sqi2fq ${bindir}/tailseq-sigproc

SIGPROC_OBJECTS=	\
	sigproc/tailseq-sigproc.o \
	sigproc/cifreader.o \
	sigproc/bclreader.o \
	sigproc/altcalls.o \
	sigproc/phix_control.o \
	sigproc/controlaligner.o \
	sigproc/my_strstr.o \
	sigproc/findpolya.o \
	sigproc/signalproc.o \
	sigproc/misc.o \
	sigproc/ssw.o

SQI2FQ_OBJECTS=	\
	sqi2fq/sqi2fq.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${DFLAGS} ${INCLUDES} $< -o $@

all: ${bindir}/sqi2fq ${bindir}/tailseq-sigproc
#	python setup.py install --user
#	python test.py

clean:
	rm -f ${SIGPROC_OBJECTS} ${SQI2FQ_OBJECTS}

distclean: clean
	rm -f ${PROG}

${bindir}/sqi2fq: ${SQI2FQ_OBJECTS}
	mkdir -p ${bindir}
	${CC} ${CFLAGS} -o $@ $<

${bindir}/tailseq-sigproc: ${SIGPROC_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${SIGPROC_OBJECTS} ${LIBS}

