#
# Copyright (c) 2014 Institute for Basic Science, Korea
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

#
# Refer "cntl2fastq Conversion User Guide" by Illumina for the descriptions
# of the file formats.
#
# Actually, MiSeq and HiSeq don't mark flags in the .control files when we
# run "GenerateFASTQ" workflows. All records are filled with just zeros.
#

const CONTROL_HEADER_SIZE = 12
const CONTROL_FILE_VERSION = 2


type ControlInfoPerCluster
    encoded::Uint16
end


type ControlFile
    path::String

    handle::IOStream
    ncluster::Uint32
    position::Uint32

    function ControlFile(filename::String)
        handle = open(filename)

        zeropad = read(handle, Uint32)
        if zeropad != 0
            error("The first four bytes of a control file should be padded with zero.")
        end 

        version = ltoh(read(handle, Uint32))
        if version != CONTROL_FILE_VERSION
            error("Control file version $version is not supported.")
        end

        ncluster = ltoh(read(handle, Uint32))
        return new(filename, handle, ncluster, 0)
    end
end


function read(cntl::ControlFile, count::Int = -1)
    if count < 0 
        nrecords = cntl.ncluster - cntl.position
    else
        nrecords = minimum((cntl.ncluster - cntl.position, count))
    end 

    records = ltoh(read(cntl.handle, Uint16, int(nrecords)))
    cntl.position += nrecords

    return map(ControlInfoPerCluster, records)
end


function seek(cntl::ControlFile, newposition::Int)
    seek(cntl.handle, CONTROL_HEADER_SIZE + sizeof(Uint16) * newposition)
    cntl.position = newposition
end

