# Tailseeker 3.1

Tailseeker is the official pipeline for TAIL-seq, which measures poly(A) tail
lengths and 3′-end modifications with Illumina SBS sequencers.

It is not yet fully stable for generic uses. Please feel free to ask anything
to Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt; whenever you are stuck on a problem while
using it. 


## Analysis levels

Users can choose the extent of analysis by Tailseeker to let Tailseeker do
almost everything, or just minimal tail length measurement. The options and
the list of supported genomes are as followed:

| Level | Genomes | Features |
| ----- | ------- | -------- |
| 1 | Any | Poly(A) length measurement (≥ 5nt)<br>Non-A additions to poly(A) tails<br>PCR duplicate removal<br>Quality check for poly(A) length measurement |
| 2 | BDGP6 *(D. melanogaster)*<br>JGIxl91 *(Xenopus laevis)* | All features from level 1<br>Poly(A) length refinement based on genome sequence<br>Non-templated 3′-end tails<br>Alignments to genome (BAM) |
| 3 | GRCh38 *(Homo sapiens)*<br>GRCm38 *(Mus musculus)*<br>GRCz10 *(Danio rerio)*<br>WBcel235 *(C. elegans)*<br>Rnor\_6.0 *(Rattus norvegicus)* | All features from level 2<br>Gene-level statistics for poly(A) length and non-templated additions<br>Gene-level quantifications |


## Running with Docker

If you have a host running [Docker](https://www.docker.com), you can run
the tailseeker pipeline without installing any. The current image is not
ready for running it on multi-node clusters. For those environments, you're
encouraged to install the software in conventional way as described later.

Download the image and a wrapper script:

    docker pull hyeshik/tailseeker:latest

    wget http://bit.ly/tseek-docker
    chmod 755 tseek-docker

Set the environment variables up:

    # Point the directory holding the raw data from an Illumina sequencer
    export TAILSEEKER_DATADIR=/storage/150922_M01178_0123_00000000-ACB72

    # Create an empty directory for new temporary and output files
    mkdir myproject    # replace myproject with your favorite name
    cd myproject

    # Download a configuration template from the repository or
    # copy from one of your previous projects
    wget -O tailseeker.yaml http://bit.ly/tseek-settings-miseq

    # Modify the settings to adapt to the current task
    vi tailseeker.yaml    # or with your preferred text editor

    # -> Change the value of "dir" under the "sources" section to
    #    "/data" where the raw data directory is mounted in our Docker
    #    container.
    # -> Edit the list of samples, barcodes and spike-ins as you need.

Run the pipeline:

    ../tseek-docker run -j

Then, the results will be located in the current directory.


## Prerequisite tools for the conventional installation

  * Essential dependencies
    * Python 3.3 or higher
    * pkg-config
    * bash
    * wget
    * make and a C compilation toolchain
    * whiptail
    * [htslib](http://www.htslib.org) – `htslib` depends on `zlib` 1.2.4 or later. If you are using an old system released before 2010, you may need to upgrade `zlib` first.
    * Python packages that can be easily installed using `pip` (see below)
      * [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) - 3.5 or higher
      * [colormath](https://pypi.python.org/pypi/colormath/)
      * [matplotlib](http://matplotlib.org)
      * [NumPy](http://numpy.org)
      * [SciPy](http://www.scipy.org)
      * [pandas](http://pandas.pydata.org)
      * [PyYAML](http://pyyaml.org)
  * Required only for optional gene-level statistics
    * [STAR](https://github.com/alexdobin/STAR)
    * [samtools](https://github.com/samtools/samtools)
    * [bedtools](https://github.com/arq5x/bedtools2)
    * [seqtk](https://github.com/lh3/seqtk)
    * [GNU parallel](http://www.gnu.org/software/parallel/)
    * [feather](https://pypi.python.org/pypi/feather-format)
    * [Python lzma module](https://docs.python.org/3/library/lzma.html)
    * [XlsxWriter](https://pypi.python.org/pypi/XlsxWriter)
  * Optional for more sensitive analysis
    * [All Your Bases](http://www.ebi.ac.uk/goldman-srv/AYB/) - requires
      [my patch](https://github.com/hyeshik/AYB2) to work with the recent Illumina
      sequencers.
    * [GSNAP](http://research-pub.gene.com/gmap/)

The toolchains and generic command line utilities can be installed if
you're an administrator on a Debian or Ubuntu system:

    sudo apt install whiptail pkg-config gcc wget make

You can install the Python modules in the list with this command from the
top source directory:

    pip3 install --user -r install/requirements.txt


## Installing

A script in the top directory will check the paths of prerequisite tools and 
guide you to set configurations correctly. Please run:

    ./setup.sh

Then, build reference databases unless you're going to run `tailseeker` in
genome-independent mode, or the level 1 analysis.

    cd refdb/level2 && snakemake -j -- {genome}
    cd refdb/level3 && snakemake -j -- {genome}

Type the identifier of the genome to be used in place of `{genome}`. List of
the available genomes are shown in the first section of this tutorial.


## Running the pipeline

  1. Copy the full output hierarchy from MiSeq or HiSeq to somewhere in
     your machine.
  2. Create an empty work directory. This is used for storing the final result
     files and the intermediate files which you may want to look into when
     something went wrong.
  3. Copy `templates/miseq-v2.yaml` to the work directory as a new name `tailseeker.yaml`.
  4. Edit the copied setting file, `tailseeker.yaml`, to change the options and
     paths to the directories. More options are specified in `conf/default-miseq.conf` and
     `conf/defaults.conf`. They are overridden if you put lines in your `tailseeker.yaml`.
  5. Run the pipeline with one of these commands:
 
     ```sh
     # In case you have an access to a job queuing system of a cluster. Change 150 to the
     # maximum number of jobs that you can put into the queue at a time.
     tseek run -c qsub -j 150

     # In case you have a single multi-core machine,
     tseek run -j
     ```

     All [Snakemake options](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-all-options)
     can be used in `tseek run`, too.

  7. Take a look at the `qcplots/` on the work directory. The plots there show how
     poly(A) length calling was accurate.

  8. Perform the downstream analyses using the output files.


## Read name format (analysis level 1 only)

In an analysis level 1 output, FASTQ files are fulfilled with nucleotide
sequences, quality scores as well as poly(A) tail information in read name.
For the higher analysis levels, read names only include minimal identifiers.
Use `refined-taginfo/*.txt.gz` for tailing status of each read this case.

FASTQ files are located in `fastq/` in the level 1 analysis. It will contain
`_R5.fastq.gz` and `_R3.fastq.gz` files for each sample. `_R5` includes
the sequences from the 5′-end of the RNA fragments, which is generally
sequenced by read 1. `_R3` is from the other end.
     
Each sequence entry has the identifier names in the following structure:

    ```
    (1)   (2)      (3) (4) (5)(6)
    a1101:00003863:0012:17:10:TT
    
    (1) Tile number with an internal lane identifier.
    (2) Serial number of the sequence, which is unique in the tile.
    (3) Flags in hexadecimal representing data processing procedure of the read.
    (4) Length of poly(A) tail.
    (5) Length of additions modifications to poly(A).
    (6) Post-poly(A) nucleotide additions.
    ```

Flags on the third field are encoded by combinations of the following bits:

| Bit (decimal) | Bit (hexadecimal) | Description |
| ------------- | ----------------- | ----------- |
|       1       |      0x0001       | A poly(A) tail is detected |
|       2       |      0x0002       | Delimiter sequence is matched with one or more mismatch |
|       4       |      0x0004       | Have a post-poly(A) modification |
|       8       |      0x0008       | Poly(A) length is measured using fluorescence signal |
|      16       |      0x0010       | Index sequence is matched to a sample with one or more mismatches |
|      32       |      0x0020       | Delimiter sequence is found at a shifted position |
|      64       |      0x0040       | One or more cycle in 3′-read are dark (no fluorescence signal) |
|     128       |      0x0080       | Delimiter sequence is not found |
|     256       |      0x0100       | Basecalling quality of balancer region is bad |
|     512       |      0x0200       | Nucleotide composition of balancer region is biased |
|    1024       |      0x0400       | Fluorescence signal in balancer region is irregular or too dark |
|    2048       |      0x0800       | Number of dark cycle in read 2 exceeds the threshold |
|    4096       |      0x1000       | (level 2) 5′-read and 3′-read are aligned to two very distant positions in the genome |
|    8192       |      0x2000       | (level 2) 3′-read is aligned to a position adjacent to an expected polyadenylation site |


# Software licenses

## The tailseeker suite

Copyright (c) 2013-2016 Hyeshik Chang

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.


## `strstr` implementation from the FreeBSD libc (`src/contrib/my_strstr.c`)

Copyright (c) 1990, 1993
     The Regents of the University of California.  All rights reserved.

This code is derived from software contributed to Berkeley by
Chris Torek.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the University nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.


## SIMD Smith-Waterman alignment library (`src/contrib/ssw.c` and `src/contrib/ssw.h`)

Copyright (c) 2012-2015 Boston College

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.



## INIH configuration file parser (`src/contrib/ini.c` and `src/contrib/ini.h`)

Copyright (c) 2009, Ben Hoyt
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of Ben Hoyt nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY BEN HOYT ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL BEN HOYT BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
