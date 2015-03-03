#!/bin/sh
java -cp /atp/hyeshik/p/tailseeker/ruler/src/picard/picard-1.91.jar:/atp/hyeshik/p/tailseeker/ruler/src/picard/sam-1.91.jar:/atp/hyeshik/p/tailseeker/ruler/src net.sf.picard.illumina.parser.TailSeqCIFDumper $*
