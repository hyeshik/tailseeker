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

TOPDIR="$(pwd)"
CONFDIR="$(pwd)/conf"
PYTHON=python3

executable() {
  found="$(which $2 2>/dev/null)"

  if [ ! -x "$found" ]; then
    if [ -z "$3" ]; then
      PROGNAME="$2"
    else
      PROGNAME="$3"
    fi

    if [ -z "$4" ]; then
      MESSAGE="Please install $PROGNAME."
    else
      MESSAGE="$4"
    fi

    echo "> An executable \`$2' is not found. $MESSAGE"
    return 2
  fi

  echo "* $2 => ${found}"
  echo "$2: ${found}" >> "$1"
}

python3mod() {
  if [ -z "$2" ]; then
    MODNAME="$1"
  else
    MODNAME="$2"
  fi

  if ! "$PYTHON" -c "import $1" >/dev/null 2>&1; then
    echo "> A Python 3 module $MODNAME could not be found. $3"
    return 2;
  fi

  modulepath="$($PYTHON -c 'import '$1'; print('$1'.__file__)')"
  echo "* $1 => ${modulepath}"
}

check_requirements() {
  echo "==> Checking required external programs ...\n"

  PATHSCONFTMP="$CONFDIR/.paths.conf"

  echo "tailor: $(pwd)" > "$PATHSCONFTMP"

  if executable "$PATHSCONFTMP" "python3" "Python 3" && \
     executable "$PATHSCONFTMP" "wget" && \
     executable "$PATHSCONFTMP" "make" && \
     executable "$PATHSCONFTMP" "whiptail" && \
     executable "$PATHSCONFTMP" "snakemake" && \
     executable "$PATHSCONFTMP" "bedtools" && \
     executable "$PATHSCONFTMP" "samtools" && \
     executable "$PATHSCONFTMP" "bgzip" "htslib" && \
     executable "$PATHSCONFTMP" "tabix" "htslib" && \
     executable "$PATHSCONFTMP" "AYB" "AYB (All Your Bases)"; then
    EXTERNAL_PROGRAMS_READY=yes
  else
    return 2
  fi

  echo "\n==> Checking required Python 3 modules ...\n"

  if python3mod "ghmm" "GHMM" "Try running support/install-ghmm.sh."; then
    GHMM_READY=yes
  else
    return 2
  fi

  if python3mod "numpy" "NumPy" && \
     python3mod "scipy" "SciPy" && \
     python3mod "Bio" "BioPython" && \
     python3mod "sklearn" "scikit-learn" && \
     python3mod "matplotlib" && \
     python3mod "yaml" "PyYAML" && \
     python3mod "snakemake"; then
    PYTHON_MODULES_READY=yes
  else
    cat - <<END
One or more required Python modules are not available. Consider running
this command:

  pip3 install -r src/requirements.txt

If you are not a superuser in this machine, try this:

  pip3 install --user -r src/requirements.txt

END
    return 2
  fi

  if [ "$EXTERNAL_PROGRAMS_READY" = "yes" -a \
       "$GHMM_READY" = "yes" -a \
       "$PYTHON_MODULES_READY" = "yes" ]; then
    echo "\nAll prerequisites are ready.\n"
    mv -f ${PATHSCONFTMP} ${CONFDIR}/paths.conf
  else
    return 2
  fi
}

build_internal_programs()
{
  cd "$TOPDIR/src"

  echo ""
  rm -rf "$TOPDIR/bin"
  make || return 1
  $PYTHON setup.py build $1 || return 1
}

install_internal_programs()
{
  cd "$TOPDIR/src"

  echo ""
  make clean || return 1
  $PYTHON setup.py install $1 || return 1
}


# ----


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


Press enter to start to setup the package.
END

read __

if ! check_requirements; then
  exit 2
fi


cat - <<END
====================================================================
We're going to build and install few small C programs and Python
extension modules written in C. This process will require a compiler
and development files like header files. Please specify additional
arguments to Python distutils if you need. Generally, none is
required for a superuser, but you'll need to add --user if you don't
have write permissions to the system. You may need to add
$HOME/.local/bin to your PATH when you are installing these
in the "user" mode.

END

echo -n "command> $PYTHON setup.py install "

read distutils_options

if ! build_internal_programs "$distutils_options"; then
  echo "\nFailed while building internal programs. Please find the error"
  echo "messages above.\n"
  return 2
fi

if ! install_internal_programs "$distutils_options"; then
  echo "\nFailed while installing internal programs. Please find the error"
  echo "messages above.\n"
  return 2
fi


# ex: ts=8 sw=2 sts=2 et
