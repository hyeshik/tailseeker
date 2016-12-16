# Tailseeker 3.1

Tailseeker is the official pipeline for TAIL-seq, which measures poly(A) tail
lengths and 3′-end modifications with Illumina SBS sequencers.

It is not yet fully stable for generic uses. Please feel free to ask anything
to Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt; whenever you are stuck on a problem while
using it. 


# Table of Contents

   * [Analysis levels](#analysis-levels)
   * [Running with Docker](#running-with-docker)
   * [Non-Docker installations](#non-docker-installations)
      * [Installing tailseeker](#installing-tailseeker)
         * [Installing from a full binary bundle](#installing-from-a-full-binary-bundle)
         * [Installing from a source package](#installing-from-a-source-package)
            * [Prerequisite software](#prerequisite-software)
            * [Installation](#installation)
      * [Running the pipeline](#running-the-pipeline)
         * [Generating genome reference databases](#generating-genome-reference-databases)
         * [Running the pipeline](#running-the-pipeline-1)
   * [Pre-built genome resource packages](#pre-built-genome-resource-packages)
   * [Data outputs](#data-outputs)
      * [Read name format (analysis level 1 only)](#read-name-format-analysis-level-1-only)
   * [Software licenses](#software-licenses)
      * [The tailseeker suite](#the-tailseeker-suite)
      * [strstr implementation from the FreeBSD libc (<code>src/contrib/my_strstr.c</code>)](#strstr-implementation-from-the-freebsd-libc-srccontribmy_strstrc)
      * [SIMD Smith-Waterman alignment library (src/contrib/ssw.c and <code>src/contrib/ssw.h</code>)](#simd-smith-waterman-alignment-library-srccontribsswc-and-srccontribsswh)
      * [INIH configuration file parser (src/contrib/ini.c and <code>src/contrib/ini.h</code>)](#inih-configuration-file-parser-srccontribinic-and-srccontribinih)


# Analysis levels

Users can choose the extent of analysis by Tailseeker to let Tailseeker do
almost everything, or just minimal tail length measurement. The options and
the list of supported genomes are as followed:

| Level | Genomes | Features |
| ----- | ------- | -------- |
| 1 | Any | Poly(A) length measurement (≥ 5nt)<br>Non-A additions to poly(A) tails<br>PCR duplicate removal<br>Quality check for poly(A) length measurement |
| 2 | BDGP6 *(D. melanogaster)*<br>JGIxl91 *(Xenopus laevis)* | All features from level 1<br>Poly(A) length refinement based on genome sequence<br>Non-templated 3′-end tails<br>Alignments to genome (BAM) |
| 3 | GRCh38 *(Homo sapiens)*<br>GRCm38 *(Mus musculus)*<br>GRCz10 *(Danio rerio)*<br>WBcel235 *(C. elegans)*<br>Rnor\_6.0 *(Rattus norvegicus)* | All features from level 2<br>Gene-level statistics for poly(A) length and non-templated additions<br>Gene-level quantifications |


# Running with Docker

If you have a host running [Docker](https://www.docker.com), you can run
the tailseeker pipeline without installing any. For Apple macOS or
Microsoft Windows users, this is the only easy way to run tailseeker
without extensive effort. The current image is not ready for running
it on multi-node HPC clusters. For those environments, you're encouraged
to install the software in
[conventional way as described later](#non-docker-installations).

Download the image and a wrapper script:

    docker pull hyeshik/tailseeker:latest

    curl -L http://bit.ly/tseek-docker > tseek
    chmod 755 tseek

Prepare a project configuration on
[this page](http://hyeshik.github.io/tailseeker/generate-settings.html).
Fill `/data` in “Data dir.” instead of the original paths.

Set the environment variables up:

    # Point the directory holding the raw data from an Illumina sequencer
    export TAILSEEKER_DATADIR=/storage/150922_M01178_0123_00000000-ACB72

    # Create an empty directory for new temporary and output files
    mkdir myproject    # replace myproject with your favorite name
    cd myproject
    cat > tailseeker.yaml
    # and paste the content generated from the settings web page.
    # Press Ctrl-D.

Run the pipeline:

    ../tseek -j

Then, the results will be located in the current directory.

When you run an analysis with references to the genome (level 2 and 3),
you need to extend the Docker image to supplement a genome reference
database. Build the Docker image from an empty directory like this:

    curl -L http://bit.ly/Dockerfile-withref > Dockerfile
    docker build -t tailseeker:GRCz10 --build-arg genome=GRCz10 .

Then, you'll need to define an environment variable before running
`tseek` to use your own Docker image.

    export TAILSEEKER_IMAGE=tailseeker:GRCz10


# Non-Docker installations

## Installing tailseeker

You can install *tailseeker* from either a source distribution or a binary package.
The binary package includes many of pre-compiled external programs that were built on
a x64 Linux box with Ubuntu 16.04. For the other environments, it is recommended to
use the source package to install it.

### Installing from a full binary bundle

Download a tarball from the [download section](https://github.com/hyeshik/tailseeker/releases).
Extract the files into an appropriate place inside your filesystem.

    wget {the download URL}
    tar -xzf tailseeker-3.x.x-bundle-ubuntu_xenial.tar.gz
    cd tailseeker-3.x.x-bundle-ubuntu_xenial

Install Python modules that are used in tailseeker using `pip`.

    pip3 install --user --upgrade --requirement install/requirements.txt

Add the `bin/` subdirectory of the tailseeker top directory to your `PATH`. To continue
using tailseeker later, you will need to add this to a shell startup script such as
`.bashrc` or `.zshrc` according to your login shell.

    export PATH="{PATH_TO}/tailseeker-3.x.x-bundle-ubuntu_xenial/bin:$PATH"

Now, you can invoke the tailseeker pipeline with `tseek` command from anywhere. Proceed to
[generate the genome reference database](#generating-genome-reference-databases).

### Installing from a source package

#### Prerequisite software

Here're the list of software that must be installed before using tailseeker.

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


#### Installation

A script in the top directory will check the paths of prerequisite tools and 
guide you to set configurations correctly. Please run:

    ./setup.sh

Proceed to [generate the genome reference database](#generating-genome-reference-databases).


## Running the pipeline

### Generating genome reference databases

First of all, build reference databases unless you're going to run `tailseeker` in
genome-independent mode, or the level 1 analysis.

    cd {tailseeker home}/refdb/level2 && snakemake -j -- {genome}
    cd {tailseeker home}/refdb/level3 && snakemake -j -- {genome}

Type the identifier of the genome to be used in place of `{genome}`. List of
the available genomes are shown in the first section of this tutorial.

### Running the pipeline

  1. Copy the full output hierarchy from MiSeq or HiSeq to somewhere in
     your machine.
  2. Create an empty work directory. This is used for storing the final result
     files and the intermediate files which you may want to look into when
     something went wrong.
  3. Prepare a settings file on
     [this page](http://hyeshik.github.io/tailseeker/generate-settings.html). Paste the
     content into a new file `tailseeker.yaml` inside the work directory.
  4. Run the pipeline with one of these commands:
 
     ```sh
     # In case you have an access to a job queuing system of a cluster. Change 150 to the
     # maximum number of jobs that you can put into the queue at a time.
     tseek -c qsub -j 150

     # In case you have a single multi-core machine,
     tseek -j
     ```

     All [Snakemake options](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-all-options)
     can be used in `tseek`, too.

  5. Take a look at the `qcplots/` on the work directory. The plots there show how
     poly(A) length calling was accurate.

  6. Perform the downstream analyses using the output files.


# Pre-built genome resource packages

Instead of building a resource database by yourself, you can download one of the
pre-built packages that are updated from time to time. Here're are the pointers
for those files.

| Date            | Species                  | Genome                                                               | Download |
| --------------- | ------------------------ | -------------------------------------------------------------------- | -------- |
| Dec 5, 2016     | *Danio rerio*            | [GRCz10](http://www.ensembl.org/Danio_rerio/Info/Index)              | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.157175.svg)](https://doi.org/10.5281/zenodo.157175) |
| Dec 5, 2016     | *Drosophila melanogater* | [BDGP6](http://www.ensembl.org/Drosophila_melanogaster/Info/Index)   | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192528.svg)](https://doi.org/10.5281/zenodo.192528) |
| Dec 5, 2016     | *Caenorhabditis elegans* | [WBcel235](http://www.ensembl.org/Caenorhabditis_elegans/Info/Index) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192527.svg)](https://doi.org/10.5281/zenodo.192527) |
| Dec 15, 2016    | *Mus musculus*           | [GRCm38](http://www.ensembl.org/Mus_musculus/Info/Index) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.203939.svg)](https://zenodo.org/record/203939) |
| Dec 15, 2016    | *Homo sapiens*           | [GRCh38](http://www.ensembl.org/Homo_sapiens/Info/Index) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.204805.svg)](https://doi.org/10.5281/zenodo.204805)<br/> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.203943.svg)](https://doi.org/10.5281/zenodo.203943) |
| Dec 16, 2016    | *Xenopus laevis*         | [JGIxl91](http://www.xenbase.org/common/displayGBrowse.do?source=xl9_1) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.205747.svg)](https://doi.org/10.5281/zenodo.205747) |


# Data outputs

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
