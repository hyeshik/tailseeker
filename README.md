# Tailseeker 2

Tailseeker is the official pipeline for TAIL-seq, which measures poly(A) tail
lengths and 3′-end modifications with Illumina SBS sequencers.

It is not yet fully stable for generic uses. Please feel free to ask anything
to Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt; whenever you are stuck on a problem while
using it. 

## Prerequisite tools

  * Python 3.3 or higher
  * Generic build tools
    * pkg-config
    * bash
    * wget
    * make and a C compiler toolchain
    * whiptail
  * Command line tools for bioinformatics
    * [htslib](http://www.htslib.org)
    * [All Your Bases](http://www.ebi.ac.uk/goldman-srv/AYB/) - requires
      [my patch](https://github.com/hyeshik/AYB2) to work with the recent Illumina
      sequencers.
  * Python modules
    * [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) - 3.5 or higher
    * [GHMM](http://ghmm.org) – use
      [my helper script](https://github.com/hyeshik/tailseeker/blob/master/support/install-ghmm.sh) to
      install it for Python 3.
    * [NumPy](http://numpy.org)
    * [SciPy](http://www.scipy.org)
    * [pandas](http://pandas.pydata.org)
    * [BioPython](http://biopython.org/wiki/Main_Page)
    * [scikit-learn](http://scikit-learn.org/stable/)
    * [matplotlib](http://matplotlib.org)
    * [PyYAML](http://pyyaml.org)

The toolchains and generic command line utilities can be installed if
you're managing a Ubuntu box:

    sudo apt-get install whiptail pkg-config gcc wget make

You can install the Python modules in the list with this command from the
top source directory:

    pip3 install --user -r src/requirements.txt


## Installing

A script in the top directory will check the paths of prerequisite tools and 
guide you to set configurations correctly. Please run:

    ./setup.sh


## Running the pipeline

  1. Copy the full output hierarchy from MiSeq or HiSeq to somewhere in
     your machine.
  2. Create an empty scratch directory. This is used for storing temporary files
     which is accessed by usually heavily I/O-bound tasks. It is recommended
     to locate this in a SSD or a mirrored RAID volume.
  3. Create an empty work directory. This is used for storing the final result
     files and the intermediate files which you may want to look into when
     something went wrong.
  4. Copy `templates/miseq-v2.yaml` to the work directory as a new name `tailseeker.yaml`.
  5. Edit the copied setting file, `tailseeker.yaml`, to change the options and
     paths to the directories. More options are specified in `conf/default-miseq.conf` and
     `conf/defaults.conf`. They are overridden if you put lines in your `tailseeker.yaml`.
  6. Run the pipeline with one of these commands:
 
     ```sh
     # In case you have an access to a job queuing system of a cluster. Change 150 to the
     # maximum number of jobs that you can put into the queue at a time.
     tailseeker run -c qsub -j 150

     # In case you have a single multi-core machine,
     tailseeker run -j
     ```

     All [Snakemake options](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-all-options)
     can be used in `tailseeker run`, too.

  7. Take a look at the `qcplots/` on the work directory. The plots there show how
     poly(A) length calling was accurate.

  8. Use FASTQ files in `fastq/` for the subsequent analyses. It will contain
     `_R5.fastq.gz` and `_R3.fastq.gz` files for each sample. `_R5` includes
     the sequences from the 5′-end of the RNA fragments, which is generally
     sequenced by read 1. `_R3` is from the other end. Each sequence entry
     has the identifier names in the following structure:

    ```
      +------------------------- Tile number with an internal lane identifier.
      |       +----------------- Serial number of the sequence, which is unique
      |       |                  in the tile.
    a1101:00003863:010:002
                    |   |
                    |   +------- Length of additions modifications to poly(A).
                    +----------- Length of poly(A) tail.
                    
    When a sequence AATTTTTTTTTTGTACGGAT is found in the _R3 file with the name
    above, it can be understood that:
     - Poly(A) length is 10 nt.
     - There is a two nt-long U tail to the poly(A) tail.
       (_R3 is in its reverse complement form.)
     - The 3′-end sequence of 3′ UTR or a transcript body is ATCCGTAC. This is
       usually not reliable when poly(A) tail is longer than ~12 nt due to
       the untrackable phasing issues in long homopolymers.
    ```


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


## `strstr` implementation from the FreeBSD libc (`src/sigproc/contrib/my_strstr.c`)

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


## SIMD Smith-Waterman alignment library (`src/sigproc/contrib/ssw.c` and `src/sigproc/contrib/ssw.h`)

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



## INIH configuration file parser (`src/sigproc/contrib/ini.c` and `src/sigproc/contrib/ini.h`)

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
