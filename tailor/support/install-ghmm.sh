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

TMPDIR="tmp-ghmm-install.$$"
GHMM_SRCDIST_DIR="ghmm-0.9-rc3"
GHMM_SRCDIST_URL="http://downloads.sourceforge.net/project/ghmm/ghmm/ghmm%200.9-rc3/ghmm-0.9-rc3.tar.gz"

ORIGDIR="$(pwd)"

if [ ! "$PYTHON" ]; then
  PYTHON=python3
fi

if [ "$($PYTHON -c 'import sys; print(sys.version_info[0])')" = "3" ]; then
  echo "This script is going to install the GHMM python package into $PYTHON."
  sleep 5
else
  echo "A python 3 executable cannot be found at \`$PYTHON'."
  exit 1
fi

rm -rf "$TMPDIR"
mkdir -p "$TMPDIR"

cd "$TMPDIR"

# Download and extract the source distribution
wget -O - "$GHMM_SRCDIST_URL" | tar -xzf -

cd "$GHMM_SRCDIST_DIR"
patch -p1 <<END
diff -ruN ghmm-0.9-rc3.orig/ghmmwrapper/ghmm.py ghmm-0.9-rc3/ghmmwrapper/ghmm.py
--- ghmm-0.9-rc3.orig/ghmmwrapper/ghmm.py	2013-08-03 02:30:46.000000000 +0900
+++ ghmm-0.9-rc3/ghmmwrapper/ghmm.py	2016-01-14 12:46:15.910745884 +0900
@@ -113,13 +113,12 @@
 import ghmmhelper
 import modhmmer
 import re
-import StringIO
+import io
 import copy
 import math
 import sys
 import os
 import logging
-from string import join
 from textwrap import fill
 
 # Initialize logging to stderr
@@ -153,6 +152,9 @@
 ghmmwrapper.ghmm_rng_init()
 ghmmwrapper.time_seed()
 
+def join(l, s):
+    return s.join(l)
+
 
 #-------------------------------------------------------------------------------
 #- Exceptions ------------------------------------------------------------------
@@ -330,7 +332,7 @@
             self._lengthOfCharacters = None
         else:
             if len(lens) == 1:
-                self._lengthOfCharacters = lens.keys()[0]
+                self._lengthOfCharacters = list(lens.keys())[0]
             else:
                 self._lengthOfCharacters = None
 
@@ -347,7 +349,7 @@
         strout = ["GHMM Alphabet:\n"]
         strout.append("Number of symbols: " + str(len(self)) + "\n")
         strout.append("External: " + str(self.listOfCharacters) + "\n")
-        strout.append("Internal: " + str(range(len(self))) + "\n")
+        strout.append("Internal: " + str(list(range(len(self)))) + "\n")
         return join(strout,'')
 
 
@@ -390,7 +392,7 @@
         """
         result = copy.deepcopy(emissionSequence)
         try:
-            result = map(lambda i: self.index[i], result)
+            result = [self.index[i] for i in result]
         except IndexError:
             raise KeyError
         return result
@@ -417,7 +419,7 @@
         """
         result = copy.deepcopy(internalSequence)
         try:
-            result = map(lambda i: self.listOfCharacters[i], result)
+            result = [self.listOfCharacters[i] for i in result]
         except IndexError:
             raise KeyError
         return result
@@ -451,7 +453,7 @@
     Creates an Alphabet with internal and external representation of range(a,b)
     @return Alphabet
     """
-    return Alphabet(range(a,b))
+    return Alphabet(list(range(a,b)))
 
 
 # To be used for labelled HMMs. We could use an Alphabet directly but this way it is more explicit.
@@ -676,7 +678,7 @@
 
 
         # check if ghmm is build with asci sequence file support
-        if isinstance(sequenceInput, str) or isinstance(sequenceInput, unicode):
+        if isinstance(sequenceInput, str) or isinstance(sequenceInput, str):
             if ghmmwrapper.ASCI_SEQ_FILE:
                 if  not os.path.exists(sequenceInput):
                     raise IOError('File ' + str(sequenceInput) + ' not found.')
@@ -921,7 +923,7 @@
 
 
         # reads in the first sequence struct in the input file
-        if isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, unicode):
+        if isinstance(sequenceSetInput, str) or isinstance(sequenceSetInput, str):
             if sequenceSetInput[-3:] == ".fa" or sequenceSetInput[-6:] == ".fasta":
                 # assuming FastA file:
                 alfa = emissionDomain.toCstruct()
@@ -936,7 +938,7 @@
                 with the conditional \"GHMM_OBSOLETE\".")
             else:
                 if not os.path.exists(sequenceSetInput):
-                    raise IOError, 'File ' + str(sequenceSetInput) + ' not found.'
+                    raise IOError('File ' + str(sequenceSetInput) + ' not found.')
                 else:
                     tmp = self.seq_read(sequenceSetInput)
                     if len(tmp) > 0:
@@ -983,7 +985,7 @@
 
 
         if seq.seq_number <= 6:
-            iter_list = range(seq.seq_number)
+            iter_list = list(range(seq.seq_number))
         else:
             iter_list = [0,1,'X',seq.seq_number-2,seq.seq_number-1]
 
@@ -1282,7 +1284,7 @@
 
     def __call__(self, fileName, modelIndex=None, filetype=None):
 
-        if not isinstance(fileName,StringIO.StringIO):
+        if not isinstance(fileName,io.StringIO):
             if not os.path.exists(fileName):
                 raise IOError('File ' + str(fileName) + ' not found.')
 
@@ -1396,7 +1398,7 @@
             if background_dist != {}:
                 # transformation to list for input into BackgroundDistribution,
                 # ensure the rigth order
-                for i in range(len(code2name.keys())-1):
+                for i in range(len(list(code2name.keys()))-1):
                     bg_list.append(background_dist[code2name[i]])
 
                 bg = BackgroundDistribution(emission_domain, bg_list)
@@ -1416,7 +1418,7 @@
 
             if background_dist != {}:
                 ids = [-1]*m.N
-                for s in hmm_dom.state.values():
+                for s in list(hmm_dom.state.values()):
                     ids[s.index-1] = s.background # s.index ranges from [1, m.N]
 
                 m.setBackground(bg, ids)
@@ -1517,7 +1519,7 @@
 
             match = res.match(line)
             if match:
-                fileLike = StringIO.StringIO(string)
+                fileLike = io.StringIO(string)
                 modelList.append(self.openSingleHMMER(fileLike))
                 string = ""
                 match = None
@@ -1677,7 +1679,7 @@
                     state.out_states, state.out_id, state.out_a = ghmmhelper.extract_out(A[i])
 
                     #set "in" probabilities
-                    A_col_i = map(lambda x: x[i], A)
+                    A_col_i = [x[i] for x in A]
                     # Numarray use A[,:i]
                     state.in_states, state.in_id, state.in_a = ghmmhelper.extract_out(A_col_i)
                     #fix probabilities in reestimation, else 0
@@ -2009,7 +2011,7 @@
             else:
                 outstr += '  '+str(i+1)
             outstr += " :(order= " + str(self.cbackground.getOrder(i))+ "): "
-            outstr += " "+join(map(f,d[i]),', ')+"\n"
+            outstr += " "+join(list(map(f,d[i])),', ')+"\n"
         return outstr
 
 
@@ -2032,11 +2034,11 @@
         distNum = self.cbackground.n
         orders = ghmmwrapper.int_array2list(self.cbackground.order, distNum)
         B = []
-        for i in xrange(distNum):
+        for i in range(distNum):
             order = orders[i]
             size = int(pow(self.m,(order+1)))
             b = [0.0]*size
-            for j in xrange(size):
+            for j in range(size):
                 b[j] = ghmmwrapper.double_matrix_getitem(self.cbackground.b,i,j)
             B.append(b)
         return (distNum,orders,B)
@@ -2054,7 +2056,7 @@
 
     def updateName2id(self):
         """adds all background names to the dictionary name2id"""
-        for i in xrange(self.cbackground.n):
+        for i in range(self.cbackground.n):
             tmp = self.cbackground.getName(i)
             if tmp is not None:
                 self.name2id[tmp] = i
@@ -2255,7 +2257,7 @@
         (alpha,scale)  = self.forward(sequence)
         beta = self.backward(sequence,scale)
 
-        return map(lambda v,w : map(lambda x,y : x*y, v, w), alpha, beta)
+        return list(map(lambda v,w : list(map(lambda x,y : x*y, v, w)), alpha, beta))
 
 
     def joined(self, emissionSequence, stateSequence):
@@ -2562,14 +2564,14 @@
         strout = []
         if model_type == kNotSpecified:
             return 'kNotSpecified'
-        for k in types.keys():
+        for k in list(types.keys()):
             if model_type & k:
                 strout.append(types[k])
         return ' '.join(strout)
 
     def updateName2id(self):
         """adds all state names to the dictionary name2id"""
-        for i in xrange(self.cmodel.N):
+        for i in range(self.cmodel.N):
             self.name2id[i] = i
             if(self.cmodel.getStateName(i) != None):
                 self.name2id[self.cmodel.getStateName(i)] = i
@@ -2632,7 +2634,7 @@
             order = [0]*hmm.N
 
         if hmm.N <= 4:
-            iter_list = range(self.N)
+            iter_list = list(range(self.N))
         else:
             iter_list = [0,1,'X',hmm.N-2,hmm.N-1]
 
@@ -2945,7 +2947,7 @@
         # check for valid background id
         for d in stateBackground:
             if type(d) == str:
-                assert self.background.name2id.has_key(d), "Error:  Invalid background distribution name."
+                assert d in self.background.name2id, "Error:  Invalid background distribution name."
                 d = self.background.name2id[d]
             assert d in range(self.background.cbackground.n), "Error: Invalid background distribution id."
 
@@ -3092,7 +3094,7 @@
         label = ghmmwrapper.int_array2list(hmm.label, self.N)
 
         if hmm.N <= 4:
-            iter_list = range(self.N)
+            iter_list = list(range(self.N))
         else:
             iter_list = [0,1,'X',hmm.N-2,hmm.N-1]
 
@@ -3250,9 +3252,9 @@
 
         f = lambda i: self.labelDomain.external(self.getLabel(i))
         if seqNumber == 1:
-            labels = map(f, vPath)
+            labels = list(map(f, vPath))
         else:
-            labels = [map(f, vp) for vp in vPath]
+            labels = [list(map(f, vp)) for vp in vPath]
 
         return (labels, log_p)
 
@@ -3290,7 +3292,7 @@
             ghmmwrapper.free(labeling)
 
         if emissionSequences.cseq.seq_number > 1:
-            return (map(self.externalLabel, allLabels), allLogs)
+            return (list(map(self.externalLabel, allLabels)), allLogs)
         else:
             return (self.externalLabel(allLabels[0]), allLogs[0])
 
@@ -3541,7 +3543,7 @@
         f = lambda x: "%.2f" % (x,)  # float rounding function
 
         if hmm.N <= 4:
-            iter_list = range(self.N)
+            iter_list = list(range(self.N))
         else:
             iter_list = [0,1,'X',hmm.N-2,hmm.N-1]
 
@@ -3757,7 +3759,7 @@
                 viterbiPath, log_p = (None, float("-infinity"))
 
             if viterbiPath != None:
-                onePath = ghmmwrapper.int_array2list(viterbiPath, seq_len/self.cmodel.dim)
+                onePath = ghmmwrapper.int_array2list(viterbiPath, int(seq_len/self.cmodel.dim))
             else:
                 onePath = []
 
@@ -3919,7 +3921,7 @@
         f = lambda x: "%.2f" % (x,)  # float rounding function
 
         if hmm.N <= 4:
-            iter_list = range(self.N)
+            iter_list = list(range(self.N))
         else:
             iter_list = [0,1,'X',hmm.N-2,hmm.N-1]
 
@@ -4934,7 +4936,7 @@
         import xml.dom.minidom
         from ghmm_gato import xmlutil
 
-        if not (isinstance(fileName_file_or_dom, StringIO.StringIO) or
+        if not (isinstance(fileName_file_or_dom, io.StringIO) or
                 isinstance(fileName_file_or_dom, xml.dom.minidom.Document)):
             if not os.path.exists(fileName_file_or_dom):
                 raise IOError('File ' + str(fileName_file_or_dom) + ' not found.')
@@ -4965,8 +4967,8 @@
             cmodel.name = 'Unused'
 
         alphabets = hmm_dom.getAlphabets()
-        cmodel.number_of_alphabets = len(alphabets.keys())
-        sizes = [len(alphabets[k]) for k in alphabets.keys()]
+        cmodel.number_of_alphabets = len(list(alphabets.keys()))
+        sizes = [len(alphabets[k]) for k in list(alphabets.keys())]
         cmodel.size_of_alphabet = ghmmwrapper.list2int_array(sizes)
 
         # set number of d_seqs to zero. If you want to use them you have to
@@ -4976,7 +4978,7 @@
         # c array of states allocated
         cstates = ghmmwrapper.dpstate_array_alloc(cmodel.N)
         # python list of states from xml
-        pystates = hmm_dom.state.values()
+        pystates = list(hmm_dom.state.values())
 
         silent_flag = 0
         silent_states = []
@@ -4985,7 +4987,7 @@
         maxOffsetY = 0
 
         transitionClassFlag = 0
-        maxTransitionIndexDiscrete = len(alphabets.keys())
+        maxTransitionIndexDiscrete = len(list(alphabets.keys()))
         maxTransitionIndexContinuous = 0
 
         # from build matrices in xmlutil:
@@ -5136,7 +5138,7 @@
         cmodel.model_type += silent_flag
         cmodel.silent = ghmmwrapper.list2int_array(silent_states)
         distribution = DiscreteDistribution(DNA)
-        emissionDomains = [Alphabet(hmm_dom.hmmAlphabets[alphabet].name.values()) for alphabet in alphabets]
+        emissionDomains = [Alphabet(list(hmm_dom.hmmAlphabets[alphabet].name.values())) for alphabet in alphabets]
         model = PairHMM(emissionDomains, distribution, cmodel)
         model.states = pystates
         model.transitionFunctions = hmm_dom.transitionFunctions
diff -ruN ghmm-0.9-rc3.orig/ghmmwrapper/ghmmhelper.py ghmm-0.9-rc3/ghmmwrapper/ghmmhelper.py
--- ghmm-0.9-rc3.orig/ghmmwrapper/ghmmhelper.py	2013-06-07 11:40:28.000000000 +0900
+++ ghmm-0.9-rc3/ghmmwrapper/ghmmhelper.py	2016-01-14 12:46:15.910745884 +0900
@@ -188,7 +188,7 @@
 
     # parsing indixes belonging to postive probabilites
     for j in range(cos):
-        transmat_col_state = map( lambda x: x[state], transmat[j])
+        transmat_col_state = [x[state] for x in transmat[j]]
         for i in range(len(transmat_col_state)):
 
             if transmat_col_state[i] != 0.0 and i not in lis:
diff -ruN ghmm-0.9-rc3.orig/ghmmwrapper/modhmmer.py ghmm-0.9-rc3/ghmmwrapper/modhmmer.py
--- ghmm-0.9-rc3.orig/ghmmwrapper/modhmmer.py	2013-06-07 11:40:30.000000000 +0900
+++ ghmm-0.9-rc3/ghmmwrapper/modhmmer.py	2016-01-14 12:46:15.910745884 +0900
@@ -33,7 +33,7 @@
 #             last change by $Author: cic99 $.
 #
 ################################################################################
-import sys,re,string,StringIO
+import sys,re,string,io
 from xml.dom import minidom
 
 def gotoLine(f,res):
@@ -51,16 +51,16 @@
 
 def build_matrix(n,m):
     "builds n x m matrix with lists"
-    matrix = range(n)
+    matrix = list(range(n))
     for i in range(n):
-        matrix[i] = range(m)
+        matrix[i] = list(range(m))
         for j in range(m):
             matrix[i][j] = 0
     return matrix
 
 def sum_mrows(mat):
     "sums the rows of a matrix"
-    mout = range(len(mat))
+    mout = list(range(len(mat)))
     for i in range(len(mat)):
         s = 0
         for j in range(len(mat[i])):
@@ -102,7 +102,7 @@
 def map_entries(dic,lis):
     "translates the letters to the number of the columns"
     dicout = {}
-    for k in dic.keys():
+    for k in list(dic.keys()):
         dicout[k] = []
         for i in range(len(dic[k])):
             dicout[k].append((lis.index(dic[k][i][0]),lis.index(dic[k][i][1])))
@@ -119,7 +119,7 @@
     "writes <strcontent> in file <strf>"
     try:
         f = open(strf,"w")
-    except IOError,info:
+    except IOError as info:
         sys.stderr.write(str(info) + "\n")
         sys.exit(1)
     try:
@@ -175,7 +175,7 @@
                 f = open(strf,"r")
             else:
                 f = strf    
-        except IOError,info:
+        except IOError as info:
                  sys.stderr.write(str(info) + "\n")
             #     sys.exit(1)
         
@@ -257,7 +257,7 @@
         #print "parsen fertig"
 		
     def __str__(self):
-        print "oben"
+        print("oben")
         hmm_str = "N= " + str(self.n)  +", M= " + str(self.m) + "\n"
         hmm_str += "Transitions: \n"
         for row in self.matrans:
@@ -265,7 +265,7 @@
         hmm_str += "Emissions: \n"		
         for row in self.maems:
             hmm_str += str(row) + "\n" 
-        print hmm_str	
+        print(hmm_str)	
         return hmm_str
 
 
@@ -296,7 +296,7 @@
         #build matrix for transition: N B E J C T M1 I1 D1 M2 I2 D2 ... Mn In Dn
         p1 = self.manull[0]
 
-        print "globalConfig:  self_p = ", p1
+        print("globalConfig:  self_p = ", p1)
         
         self.matrans[0][0] = 1.0 - p1  #  N->B
         self.matrans[0][1] = p1       # N->N
@@ -441,9 +441,9 @@
                 nodnode.appendChild(noddata)
                 #node:data:emissions
                 if k in ["M","I"]:
-                    strem = string.join(map(str,self.maems[self.lisMID.index(k)][i]),", ")
+                    strem = string.join(list(map(str,self.maems[self.lisMID.index(k)][i])),", ")
                 elif k in ["N","C"]:
-                    strem = string.join(map(str,self.maems[2][0]),", ")
+                    strem = string.join(list(map(str,self.maems[2][0])),", ")
                 else:
                     strem = (self.m-1)*"0.0, " + "0.0"
                 nodnode.appendChild(xml_newdatanode(doc,"data","key","emissions",strem))
diff -ruN ghmm-0.9-rc3.orig/ghmmwrapper/sclass_change.c ghmm-0.9-rc3/ghmmwrapper/sclass_change.c
--- ghmm-0.9-rc3.orig/ghmmwrapper/sclass_change.c	2013-06-07 11:40:30.000000000 +0900
+++ ghmm-0.9-rc3/ghmmwrapper/sclass_change.c	2016-01-14 12:46:15.910745884 +0900
@@ -45,6 +45,24 @@
 #include <ghmm/sequence.h>
 #include <ghmm/smodel.h>
 
+#if PY_VERSION_HEX >= 0x03000000
+
+#define PyClass_Check(obj) PyObject_IsInstance(obj, (PyObject *)&PyType_Type)
+#define PyInt_Check(x) PyLong_Check(x)
+#define PyInt_AsLong(x) PyLong_AsLong(x)
+#define PyInt_FromLong(x) PyLong_FromLong(x)
+#define PyInt_FromSize_t(x) PyLong_FromSize_t(x)
+#define PyString_Check(name) PyBytes_Check(name)
+#define PyString_FromString(x) PyUnicode_FromString(x)
+#define PyString_Format(fmt, args)  PyUnicode_Format(fmt, args)
+#define PyString_AsString(str) PyBytes_AsString(str)
+#define PyString_Size(str) PyBytes_Size(str)    
+#define PyString_InternFromString(key) PyUnicode_InternFromString(key)
+#define Py_TPFLAGS_HAVE_CLASS Py_TPFLAGS_BASETYPE
+#define PyString_AS_STRING(x) PyUnicode_AS_STRING(x)
+#define _PyLong_FromSsize_t(x) PyLong_FromSsize_t(x)
+
+#endif
 
 
 /* smo is a ghmm_cmodel struct
diff -ruN ghmm-0.9-rc3.orig/ghmmwrapper/setup.py ghmm-0.9-rc3/ghmmwrapper/setup.py
--- ghmm-0.9-rc3.orig/ghmmwrapper/setup.py	2013-08-06 00:55:08.000000000 +0900
+++ ghmm-0.9-rc3/ghmmwrapper/setup.py	2016-01-14 12:47:19.406419410 +0900
@@ -46,7 +46,7 @@
       ext_modules = [Extension('_ghmmwrapper',
                                ['sclass_change.c', 'pclasschange.c', 'gql.c', 'ghmmwrapper.i'],
                                include_dirs = ['..'],
-                               library_dirs = ['../ghmm/.libs'],
+                               library_dirs = ['.'],
                                libraries = ['ghmm', 'm', 'pthread', 'xml2', 'z'],
                                extra_compile_args = ["-O2", "-pipe", "-Wall"], # -g might help debugging
                                depends = ['wrapper_alphabet.i', 'wrapper_cmodel.i', 'wrapper_cseq.i',
END


./configure --without-python || (echo "<*> Failed to configure GHMM."; exit 1)
make -j4 || (echo "<*> Failed to build GHMM."; exit 1)

# Make an object archive for static linking to the python extension module.
# Don't use object files in ghmm/*.o, they are compiled with position dependency.
ar cru ghmmwrapper/libghmm.a ghmm/.libs/*.o
ranlib ghmmwrapper/libghmm.a

cd ghmmwrapper

"$PYTHON" setup.py install || (echo "<*> Failed to build and install the extension."; exit 1)

echo ""
echo "<*> Finished installing the GHMM package to your $PYTHON."

cd "$ORIGDIR"
rm -rf "$TMPDIR"
