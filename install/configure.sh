#!/bin/sh
#
# Copyright (c) 2016 Institute for Basic Science
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

BACKTITLE="Tailseeker 3"
PATHCONF="${TOPDIR}/conf/paths.conf"
PYTHON=python3

required_executables_level1="\
${PYTHON}:Python_3
pkg-config:pkg-config
bash:bash
wget:wget
make:make
snakemake:Snakemake
bgzip:htslib"
required_pkgconfig_level1="\
htslib:htslib"
required_python3mod_level1="\
colormath:colormath
matplotlib:matplotlib
numpy:numpy
pandas:pandas
scipy:scipy
snakemake:Snakemake
yaml:PyYAML"

required_executables_level2="\
tabix:htslib
STAR:STAR_RNA-seq_aligner
samtools:samtools
bedtools:bedtools
parallel:GNU_parallel"
required_pkgconfig_level2=""
required_python3mod_level2=""

required_executables_level3=""
required_pkgconfig_level3=""
required_python3mod_level3="\
lzma:lzma
feather:feather-format"


RED="[01;31m"
WHITE="[1m"
RESET="[00m"

requirements_not_found=""

# =======================================================
# Choose the analysis level
# =======================================================

analysis_level=\
$($WHIPTAIL --title "Configure: Analysis level" \
--backtitle "$BACKTITLE" --default-item "Level 3" \
--menu "
You can run Tailseeker for different ranges of provided
functions. It can be used for essential tail length
measurement only (level 1) or for getting the full gene-
level statistics (level 3). Higher level will include
all result data from lower level, but it takes much more
space, time, and external dependencies. Please choose
one. You can re-configure it later." \
18 60 3 \
"Level 1" " Tail length measurement and QC." \
"Level 2" " Refinement of calls based on the genome." \
"Level 3" " Gene-level statistics." 3>&1 1>&2 2>&3)

if [ $? -ne 0 ]; then
  exit 0
fi

analysis_level=$(echo $analysis_level | sed -e 's,Level ,,')

echo "tailseeker: $TOPDIR" > $PATHCONF
echo "whiptail: $WHIPTAIL" >> $PATHCONF
echo "analysis_level: $analysis_level" >> $PATHCONF


# =======================================================
# Choose whether to use the AYB base-caller
# =======================================================

if [ $analysis_level -ge 2 ]; then
  $WHIPTAIL --title "Configure: AYB basecalling" \
--backtitle "$BACKTITLE" \
--yesno "\
The standard basecaller of the Illumina sequencer is
sometimes fail to produce high quality basecalls for
5'-side read (usually read 1) when long poly(A) tail is
attached on its 3'-side. AYB is much stronger at reliable
basecalling in this case, gives 10-50% more mappable reads
for tags adjacent to poly(A) tails > 100 nt.

Do you want to use AYB for basecalling?" 15 62

  if [ "$?" -eq 0 ]; then
    required_executables_level1="$required_executables_level1
AYB:AYB_(All_Your_Bases)"
  fi
fi


# =======================================================
# Choose whether to use the GSNAP aligner
# =======================================================

use_gsnap=no
if [ $analysis_level -ge 2 ]; then
  $WHIPTAIL --title "Configure: GSNAP alignment" \
--backtitle "$BACKTITLE" --defaultno \
--yesno "\
Tailseeker 3 uses STAR as the primary sequence alignment
software. Optionally, GSNAP can be used to rescue unmapped
sequences by STAR. Tags with long poly(A) tails sometimes
suffer from low sequence read quality, and this leads
to failure in identification of tags. Usually, GSNAP
rescues ~0.5% of total reads and ~10% among long poly(A)
tails. However, it lengthens total running time of the
pipeline by three-fold at least.

Do you want to use GSNAP together with STAR?" 17 62
  
  if [ "$?" -eq 0 ]; then
    required_executables_level2="$required_executables_level2
gsnap:GSNAP_(GMAP)
gmap_build:GSNAP_(GMAP)
gtf_splicesites:GSNAP_(GMAP)
gtf_introns:GSNAP_(GMAP)
iit_store:GSNAP_(GMAP)"
    use_gsnap=yes
  fi
fi


# =======================================================
# Choose GSNAP indexing sensitivity
# =======================================================

if [ "$use_gsnap" = yes -a $analysis_level -ge 2 ]; then
  $WHIPTAIL --title "Configure: GSNAP sensitivity" \
--backtitle "$BACKTITLE" --defaultno \
--yesno "\
Building GSNAP genome indices with shorter k-mer and more
frequent sampling interval boosts the alignment sensitivity
by 2-3% for the difficult alignments. Again, it makes the
GSNAP tasks even slower by ~2x.

Do you want to tune GSNAP index building for better
sensitivity?" 13 64
  
  if [ "$?" -eq 0 ]; then
    echo "gsnap_sensitive_index: yes" >> $PATHCONF
  else
    echo "gsnap_sensitive_index: no" >> $PATHCONF
  fi
fi


# =======================================================
# Set up the lists of external dependencies
# =======================================================

required_executables=""
required_pkgconfig=""
required_python3mod=""

if [ $analysis_level -ge 1 ]; then
  required_executables="$required_executables_level1"
  required_pkgconfig="$required_pkgconfig_level1"
  required_python3mod="$required_python3mod_level1"
fi
if [ $analysis_level -ge 2 ]; then
  required_executables="$required_executables $required_executables_level2"
  required_pkgconfig="$required_pkgconfig $required_pkgconfig_level2"
  required_python3mod="$required_python3mod $required_python3mod_level2"
fi
if [ $analysis_level -ge 3 ]; then
  required_executables="$required_executables $required_executables_level3"
  required_pkgconfig="$required_pkgconfig $required_pkgconfig_level3"
  required_python3mod="$required_python3mod $required_python3mod_level3"
fi


# =======================================================
# Check paths of executables
# =======================================================

echo "==> Checking paths of required executables"

for req in $required_executables; do
  command=$(echo "$req" | cut -d: -f1)
  command_quoted=$(echo "$command" | sed -e 's,-,_,g')
  friendly_name=$(echo "$req" | cut -d: -f2 | sed -e 's,_, ,g')
  fullpath=$(which $command)

  if [ -z "$fullpath" ]; then
    echo "$RED * $command -> NOT FOUND$RESET"
    requirements_not_found="$requirements_not_found 
 * $friendly_name"
  else
    echo " * $command -> $fullpath"
    echo "$command_quoted: $fullpath" >> $PATHCONF
  fi
done

echo ""


# =======================================================
# Check installations of pkg-config packages
# =======================================================

echo "==> Checking installed versions of packages in pkg-config"

for req in $required_pkgconfig; do
  pkgname=$(echo "$req" | cut -d: -f1)
  friendly_name=$(echo "$req" | cut -d: -f2 | sed -e 's,_, ,g')
  modversion=$(pkg-config --modversion $pkgname)

  if [ -z "$modversion" ]; then
    echo "$RED * $pkgname -> NOT FOUND$RESET"
    requirements_not_found="$requirements_not_found 
 * $friendly_name"
  else
    echo " * $pkgname -> $modversion found"
  fi
done

echo ""


# =======================================================
# Check installations of Python 3 modules/packages
# =======================================================

echo "==> Checking Python 3 modules and packages"

for req in $required_python3mod; do
  modname=$(echo "$req" | cut -d: -f1)
  friendly_name=$(echo "$req" | cut -d: -f2 | sed -e 's,_, ,g')

  if ! "$PYTHON" -c "import $modname" >/dev/null 2>&1; then
    echo "$RED * $modname -> NOT FOUND$RESET"
    requirements_not_found="$requirements_not_found 
 * $friendly_name"
  else
    modulepath="$($PYTHON -c 'import '$modname'; print('$modname'.__file__)')"
    echo " * $modname -> $(dirname $modulepath)"
  fi
done

echo ""


# =======================================================
# Show error message and instruction if requirements were not meet
# =======================================================

if [ "$requirements_not_found" ]; then
  cat - <<END

#### ERROR ####

One or more required components to run Tailseeker is not detected
in an appropriate location. Please retry the configuration after
installing the following software. If it does not find an installed
program or module, please correct PATH, PKG_CONFIG_PATH or
PYTHONPATH environment variables.

Unmet requirements:$requirements_not_found

END
  exit 1
fi


# =======================================================
# Install a command line entry script
# =======================================================

writable_dirs=$(echo "$PATH" | awk -F: '{
  NODESCR="\x7f"
  for (i = 1; i <= NF; i++)
    if (system("test ! \\( -w \""$i"\" -a -d \""$i"\" \\)"))
      print $i " " NODESCR
}')
num_writable_dirs=$(echo "$writable_dirs" | wc -l)

if [ -z "$writable_dirs" ]; then
  cat - <<END
No directory in your PATH is writable. Please create a new
directory that will hold your own binaries, and add it to
your PATH environment variable.
END
  exit 3
fi

install_dir=$($WHIPTAIL --title "Configure: install the tailseeker command" \
--backtitle "$BACKTITLE" \
--menu "
Choose a directory where to install the \"tailseeker\" command.
It is recommended to keep the directory in your PATH." \
$(($num_writable_dirs + 10)) 65 $num_writable_dirs $writable_dirs \
3>&1 1>&2 2>&3)

rm -f $install_dir/tailseeker
cat $TOPDIR/install/tailseeker.in | \
sed -e "s,%%TAILSEEKER_DIR%%,$TOPDIR,g" > $install_dir/tailseeker
script_mode=$(printf "%04o\n" $((0777 & ~$(umask))))
chmod $script_mode $install_dir/tailseeker


# ============================
# Show the final instructions
# ============================

if [ "$analysis_level" -eq 1 ]; then
  cat - <<END

Congratulations!

We are ready to run Tailseeker. Please refer the documentation for what
to do next.

Good Luck!

Hyeshik Chang

END
else
  if [ "$analysis_level" -ge 3 ]; then
    refdbdirs="${TOPDIR}/refdb/level2 and ${TOPDIR}/refdb/level3"
  else
    refdbdirs="${TOPDIR}/refdb/level2"
  fi

  getintomsg="$(echo "Get into $refdbdirs, and run this command:" | \
                fmt -w 65)"

  cat - <<END

Congratulations!

We are almost ready to run Tailseeker.

Before running your first analysis, you need to prepare indices for the
reference genome sequences and gene annotations.

$getintomsg

  snakemake -j -- {genome}

You need to specify the name of genome to prepare in place of
{genome}. The list of available genomes are followed below.

END

  echo "In refdb/level2 (reference database for level 2 analysis):"
  (cd ${TOPDIR}/refdb/level2; env QUIET=yes snakemake -q)

  echo ""

  if [ "$analysis_level" -ge 3 ]; then
    echo "In refdb/level3 (reference database for level 3 analysis):"
    (cd ${TOPDIR}/refdb/level3; env QUIET=yes snakemake -q)
  fi

  cat - <<END

You'll need to carefully monitor the memory consumption with tools like
top(1) or htop(1) during the build in case your system does not have
plenty of memory. Building a STAR index takes up more than 30 GiB of
RAM for many genomes.

Please refer the documentation for what to do next.

Good Luck!

Hyeshik Chang

END
fi

# ex: ts=8 sw=2 sts=2 et
