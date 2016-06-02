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

TOPDIR_REL=$(dirname $0)
TOPDIR=$(cd $TOPDIR_REL; pwd)

RED="[01;31m"
WHITE="[1m"
RESET="[00m"

cat - <<END
  ______     ${RED}_${RESET} __            __              |  ${WHITE}from${RESET}
 /_  __/__ _${RED}(_)${RESET} /__ ___ ___ / /_____ ____    |  ${WHITE}Center for RNA Research${RESET}
  / / / _ \`/ / (_-</ -_) -_)  '_/ -_) __/    |  ${WHITE}Institute for Basic Science${RESET}
 /_/  \\_,_/_/_/___/\\__/\\__/_/\\_\\\\__/_/       |  ${WHITE}Seoul, South Korea${RESET}

 The most recent version of the software is available from
   https://github.com/hyeshik/tailseeker.

 Send questions, comments, bug-reports, and etc. to
   Hyeshik Chang <hyeshik@snu.ac.kr>.


END

WHIPTAIL=$(which whiptail)
if [ -z "$WHIPTAIL" ]; then
  cat - <<END
## ERROR ##

\`whiptail' is required for the installation process. Please follow one
of the following instructions:

Debian/Ubuntu Linux:
  apt-get install whiptail

RedHat Linux or CentOS:
  yum install -y whiptail

Conda:
  conda install -c bioconda newt

END
  exit 1
fi

BASH=$(which bash)
if [ -z "$BASH" ]; then
  cat - <<END
## ERROR ##

\`bash' is required for the installation process. Please follow one
of the following instructions:

Debian/Ubuntu Linux:
  apt-get install bash

RedHat Linux or CentOS:
  yum install -y bash

Conda:
  conda install -c bioconda bash

END
  exit 1
fi

echo "Press enter to start to setup the package."
read __

env "TOPDIR=$TOPDIR" "WHIPTAIL=$WHIPTAIL" "SHELL=$BASH" "$BASH" "$TOPDIR/install/configure.sh"

# ex: ts=8 sw=2 sts=2 et
