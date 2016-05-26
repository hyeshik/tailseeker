CC?=		gcc
CXX?=		g++
#CFLAGS=	-O0 -Wall -Werror -pg -ggdb ${HTSLIB_CFLAGS}
#CXXFLAGS=	-O0 -DNO_OPENMP -Wall -pg -ggdb
CFLAGS=		-O3 -Wall -DNDEBUG ${HTSLIB_CFLAGS}
CXXFLAGS=	-O3 -DNO_OPENMP -DNDEBUG

SIGPROC_LIBS=	-lz -lm -lpthread ${HTSLIB_LIBS}
DEDUP_PERFECT_LIBS=	-lm
DEDUP_APPROX_LIBS=	-lm -lpthread ${HTSLIB_LIBS}
WRITEFASTQ_LIBS=	-lm -lz ${HTSLIB_LIBS}
ARCH_FLAGS=	-msse2 -DUSE_SSE2
bindir=		../bin

HTSLIB_CFLAGS:=	$(shell pkg-config --cflags htslib)
HTSLIB_LIBS:=	$(shell pkg-config --libs htslib)


PROG=	${bindir}/tailseq-sigproc ${bindir}/tailseq-dedup-perfect \
	${bindir}/tailseq-writefastq ${bindir}/tailseq-dedup-approx

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

DEDUP_PERFECT_OBJECTS= \
	tailseq-dedup-perfect.o

DEDUP_APPROX_OBJECTS= \
	tailseq-dedup-approx.o \
	cdhit-proxy.o \
	cdhit/cdhit-common-patched.o

WRITEFASTQ_OBJECTS= \
	utils.o \
	tailseq-writefastq.o

.SUFFIXES:.c .cc .o

.c.o:
	${CC} -c ${CFLAGS} ${ARCH_FLAGS} ${INCLUDES} $< -o $@

.cc.o:
	${CXX} -c ${CXXFLAGS} ${INCLUDES} $< -o $@

all: ${PROG}

clean:
	rm -f ${SIGPROC_OBJECTS} ${DEDUP_PERFECT_OBJECTS} ${WRITEFASTQ_OBJECTS} \
		${DEDUP_APPROX_OBJECTS}
	rm -rf cdhit

distclean: clean
	rm -f ${PROG}

cdhit/cdhit-common-patched.cc: cdhit/cdhit-common.c++ cdhit/cdhit-common.h
	patch -f -o $@ $< cdhit-patches.diff

cdhit-proxy.o: cdhit/cdhit-common.h

cdhit/cdhit-common.c++:
	mkdir -p cdhit
	wget -O $@ 'https://raw.githubusercontent.com/weizhongli/cdhit/master/cdhit-common.c%2B%2B'

cdhit/cdhit-common.h:
	mkdir -p cdhit
	wget -O $@ 'https://raw.githubusercontent.com/weizhongli/cdhit/master/cdhit-common.h'

${bindir}/tailseq-sigproc: ${SIGPROC_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${SIGPROC_OBJECTS} ${SIGPROC_LIBS}

${bindir}/tailseq-dedup-perfect: ${DEDUP_PERFECT_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${DEDUP_PERFECT_OBJECTS} ${DEDUP_PERFECT_LIBS}

${bindir}/tailseq-dedup-approx: ${DEDUP_APPROX_OBJECTS}
	${CXX} ${CFLAGS} -o $@ ${DEDUP_APPROX_OBJECTS} ${DEDUP_APPROX_LIBS}

${bindir}/tailseq-writefastq: ${WRITEFASTQ_OBJECTS}
	${CC} ${CFLAGS} -o $@ ${WRITEFASTQ_OBJECTS} ${WRITEFASTQ_LIBS}

