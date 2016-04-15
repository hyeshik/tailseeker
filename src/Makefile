CC?=		gcc
#CFLAGS=	-O0 -Wall -Werror -pg -ggdb ${HTSLIB_CFLAGS}
CFLAGS=		-O3 -Wall -Werror -DNDEBUG ${HTSLIB_CFLAGS}
SIGPROC_LIBS=	-lz -lm -lpthread ${HTSLIB_LIBS}
DEDUP_LIBS=	-lm
WRITEFASTQ_LIBS=-lm -lz ${HTSLIB_LIBS}
bindir=		../bin

HTSLIB_CFLAGS:=	$(shell pkg-config --cflags htslib)
HTSLIB_LIBS:=	$(shell pkg-config --libs htslib)

PROG=	${bindir}/tailseq-sigproc ${bindir}/tailseq-dedup \
	${bindir}/tailseq-writefastq

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
	utils.o \
	contrib/ini.o \
	contrib/misc.o \
	contrib/my_strstr.o \
	contrib/ssw.o

DEDUP_OBJECTS= \
	tailseq-dedup.o

WRITEFASTQ_OBJECTS= \
	utils.o \
	tailseq-writefastq.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${INCLUDES} $< -o $@

all: ${PROG}

clean:
	rm -f ${SIGPROC_OBJECTS} ${DEDUP_OBJECTS} ${WRITEFASTQ_OBJECTS}

distclean: clean
	rm -f ${PROG}

${bindir}/tailseq-sigproc: ${SIGPROC_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${SIGPROC_OBJECTS} ${SIGPROC_LIBS}

${bindir}/tailseq-dedup: ${DEDUP_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${DEDUP_OBJECTS} ${DEDUP_LIBS}

${bindir}/tailseq-writefastq: ${WRITEFASTQ_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${WRITEFASTQ_OBJECTS} ${WRITEFASTQ_LIBS}
