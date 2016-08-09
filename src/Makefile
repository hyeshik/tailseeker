CC?=		gcc
#CFLAGS=	-O0 -Wall -Werror -pg -ggdb -DINI_MAX_LINE=1024 ${HTSLIB_CFLAGS}
CFLAGS=		-O3 -Wall -DNDEBUG -DINI_MAX_LINE=1024 ${HTSLIB_CFLAGS}

IMPORT_LIBS=	-lz -lm -lpthread ${HTSLIB_LIBS}
SIGPROC_LIBS=	-lz -lm -lpthread ${HTSLIB_LIBS}
DEDUP_PERFECT_LIBS=	-lm
DEDUP_APPROX_LIBS=	-lm -lpthread ${HTSLIB_LIBS}
WRITEFASTQ_LIBS=	-lm -lz ${HTSLIB_LIBS}
ARCH_FLAGS=	-msse2 -DUSE_SSE2
bindir=		../bin

HTSLIB_CFLAGS:=	$(shell pkg-config --cflags htslib)
HTSLIB_LIBS:=	$(shell pkg-config --libs htslib)


PROG=	${bindir}/tailseq-import ${bindir}/tailseq-sigproc \
	${bindir}/tailseq-dedup-perfect ${bindir}/tailseq-writefastq \
	${bindir}/tailseq-dedup-approx

IMPORT_OBJECTS= \
	importer/altcalls.o \
	importer/bclreader.o \
	importer/cifreader.o \
	importer/controlaligner.o \
	importer/findpolya.o \
	importer/phix_control.o \
	importer/signalproc.o \
	importer/spotanalyzer.o \
	importer/parseconfig.o \
	importer/tailseq-import.o \
	utils.o \
	contrib/ini.o \
	contrib/misc.o \
	contrib/my_strstr.o \
	contrib/ssw.o

SIGPROC_OBJECTS= \
	sigproc/altcalls.o \
	sigproc/bclreader.o \
	sigproc/cifreader.o \
	sigproc/controlaligner.o \
	sigproc/findpolya.o \
	sigproc/parseconfig.o \
	sigproc/phix_control.o \
	sigproc/signalproc.o \
	sigproc/spotanalyzer.o \
	sigproc/tailseq-sigproc.o \
	utils.o \
	contrib/ini.o \
	contrib/misc.o \
	contrib/my_strstr.o \
	contrib/ssw.o

DEDUP_PERFECT_OBJECTS= \
	deduplicator/tailseq-dedup-perfect.o

DEDUP_APPROX_OBJECTS= \
	deduplicator/tailseq-dedup-approx.o

WRITEFASTQ_OBJECTS= \
	utils.o \
	exporter/tailseq-writefastq.o

.SUFFIXES:.c .o

.c.o:
	${CC} -c ${CFLAGS} ${ARCH_FLAGS} ${INCLUDES} $< -o $@

all: ${PROG}

clean:
	rm -f ${IMPORT_OBJECTS} ${SIGPROC_OBJECTS} ${DEDUP_PERFECT_OBJECTS} \
		${WRITEFASTQ_OBJECTS} ${DEDUP_APPROX_OBJECTS}
	rm -rf cdhit

distclean: clean
	rm -f ${PROG}

${bindir}/tailseq-import: ${IMPORT_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${IMPORT_OBJECTS} ${IMPORT_LIBS}

${bindir}/tailseq-sigproc: ${SIGPROC_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${SIGPROC_OBJECTS} ${SIGPROC_LIBS}

${bindir}/tailseq-dedup-perfect: ${DEDUP_PERFECT_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${DEDUP_PERFECT_OBJECTS} ${DEDUP_PERFECT_LIBS}

${bindir}/tailseq-dedup-approx: ${DEDUP_APPROX_OBJECTS}
	${CXX} ${CFLAGS} -fopenmp -o $@ ${DEDUP_APPROX_OBJECTS} ${DEDUP_APPROX_LIBS}

${bindir}/tailseq-writefastq: ${WRITEFASTQ_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${WRITEFASTQ_OBJECTS} ${WRITEFASTQ_LIBS}

