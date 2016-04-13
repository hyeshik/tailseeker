CC?=		gcc
#CFLAGS=	-O0 -Wall -Werror -pg -ggdb ${HTSLIB_CFLAGS}
CFLAGS=		-O3 -Wall -Werror -DNDEBUG ${HTSLIB_CFLAGS}
SIGPROC_LIBS=	-lz -lm -lpthread ${HTSLIB_LIBS}
DEDUP_LIBS=	-lm
bindir=		../bin

HTSLIB_CFLAGS:=	$(shell pkg-config --cflags htslib)
HTSLIB_LIBS:=	$(shell pkg-config --libs htslib)

PROG=	${bindir}/tailseq-sigproc

SIGPROC_OBJECTS= \
	altcalls.o \
	bclreader.o \
	cifreader.o \
	controlaligner.o \
	findpolya.o \
	parseconfig.o \
	phix_control.o \
	signalproc.o \
	spotanalyzer.o \
	tailseq-sigproc.o \
	contrib/ini.o \
	contrib/misc.o \
	contrib/my_strstr.o \
	contrib/ssw.o

DEDUP_OBJECTS= \
	tailseq-dedup.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${INCLUDES} $< -o $@

all: ${bindir}/tailseq-sigproc ${bindir}/tailseq-dedup
#	python setup.py install --user
#	python test.py

clean:
	rm -f ${SIGPROC_OBJECTS} ${DEDUP_OBJECTS}

distclean: clean
	rm -f ${PROG}

${bindir}/tailseq-sigproc: ${SIGPROC_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${SIGPROC_OBJECTS} ${SIGPROC_LIBS}

${bindir}/tailseq-dedup: ${DEDUP_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${DEDUP_OBJECTS} ${DEDUP_LIBS}

