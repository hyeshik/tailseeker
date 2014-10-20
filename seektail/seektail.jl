#!/usr/bin/env julia
#
# Copyright (c) 2013-2014 Institute for Basic Science
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

importall Base

include("bclfile.jl")
include("ciffile.jl")
include("controlfile.jl")

# ---------------

#cif = CIFFile("sample-miseq/Data/Intensities/L001/C1.1/s_1_1101.cif")
#println(length(read(cif, 100)))
#println(length(read(cif, 200)))
#println(length(read(cif)))
#println(length(read(cif)))
#seek(cif, 100000)
#println(length(read(cif)))
#seek(cif, 0)
#println(read(cif, 2))
#println("---------------")
#
#cif = CIFCollection("sample-miseq/Data/Intensities", 1, 1101)
#println(cif.ncluster)
#println(read(cif, 2))

#bcl = BCLFile("sample-miseq/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.bcl")
#println(length(read(bcl, 100)))
#println(length(read(bcl, 200)))
#println(length(read(bcl)))
#println(length(read(bcl)))
#seek(bcl, 100000)
#println(length(read(bcl)))
#seek(bcl, 0)
#println(read(bcl, 20))
#println("---------------")

#bcl = BCLCollection("sample-miseq/Data/Intensities", 1, 1101)
#println(bcl.ncluster)
#println(read(bcl, 2))
#
#bcl = BCLCollection("sample-miseq/Data/Intensities", 1, 1101, 1, 30)
#println(read(bcl, 2))
#
#cif = CIFCollection("sample-miseq/Data/Intensities", 1, 1101, 1, 30)
#println(read(cif, 2))

ctl = ControlFile("sample-miseq/Data/Intensities/BaseCalls/L001/s_1_1101.control")
records = read(ctl)
println(sum(records))

println("Done.")
