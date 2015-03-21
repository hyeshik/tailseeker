#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

shell.prefix('set -e; set -o pipefail; ')
shell.executable(os.popen('which bash').read().strip()) # pipefail is supported by bash only.

RFAM_FASTA_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/{accession}.fa.gz'

CONTAMINANTS_RFAM_ACCESSIONS = [
    'RF00001', # 5S ribosomal RNA
    'RF00002', # 5.8S ribosomal RNA
    'RF01960', # Eukaryotic small subunit ribosomal RNA
    'RF02543', # Eukaryotic large subunit ribosomal RNA
    'RF00012', # Small nucleolar RNA U3
    'RF00017', # Metazoan signal recognition particle RNA
]

rule download_rfam_fasta:
    output: 'tmp/Rfam/{accession}.fa.gz'
    run:
        url = RFAM_FASTA_URL.format(accession=wildcards.accession)
        shell('wget -O {output} {url}')

rule rfam_contaminant_seq_names:
    input: 'tmp/Rfam/{accession}.fa.gz'
    output: temp('tmp/{species}/contaminants-ids-{accession}')
    run:
        species_w_space = wildcards.species.replace('_', ' ')
        shell("zgrep '^>.*{species_w_space}' {input} | \
                sed -e 's,^>\([^ ]*\).*,\\1,g' > {output}")

rule make_contaminants_fasta:
    input:
        rfam_ids='tmp/{species}/contaminants-ids-{accession}',
        rfam_fasta='tmp/Rfam/{accession}.fa.gz'
    output: temp('tmp/{species}/contaminants-{accession}.fa')
    shell: 'faSomeRecords {input.rfam_fasta} {input.rfam_ids} {output}'

rule merge_contaminants_fasta:
    input:
        expand('tmp/{{species}}/contaminants-{accession}.fa',
               accession=CONTAMINANTS_RFAM_ACCESSIONS)
    output: 'contaminants/{species}.fa'
    shell: 'cat {input} | awk \'/^>/ {{ print $0; }} \
                                /^[^>]/ {{ gsub(/U/, "T"); print $0; }}\' > {output}'

rule build_gsnap_index:
    input: 'contaminants/{species}.fa'
    output: 'contaminants/{species}.gmap/{species}.gmap.genomecomp'
    shell: 'gmap_build -D contaminants -d {wildcards.species}.gmap -k 12 -q 1 {input}'

rule build_contaminants_star_index:
    input: 'contaminants/{species}.fa'
    output: 'contaminants/{species}.star/Genome'
    threads: 32
    params: outputdir='contaminants/{species}.star'
    shell: 'STAR --runMode genomeGenerate --genomeDir {params.outputdir} \
                --genomeSAindexNbases 4 --genomeFastaFiles {input} --runThreadN {threads}'

# ex: syntax=snakemake
