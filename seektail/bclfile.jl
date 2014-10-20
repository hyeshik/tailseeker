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
# Refer "bcl2fastq Conversion User Guide" by Illumina for the descriptions
# of the file formats.
#

const BCL_HEADER_SIZE = 4
const NOCALL_BASE = b"N"[1]
const CALL_BASES = b"ACGT"
const NOCALL_QUALITY = 2
const PHRED_BASE = 33

type BCLFile
    path::String

    handle::IOStream
    ncluster::Uint32
    position::Uint32

    function BCLFile(filename::String)
        handle = open(filename)
        ncluster = ltoh(read(handle, Uint32))
        return new(filename, handle, ncluster, 0)
    end
end


function read(bcl::BCLFile, count::Int = -1)
    if count < 0 
        nrecords = bcl.ncluster - bcl.position
    else
        nrecords = minimum((bcl.ncluster - bcl.position, count))
    end 

    records = read(bcl.handle, Uint8, int(nrecords))
    sequence = Array(Uint8, nrecords)
    quality = Array(Uint8, nrecords)

    for i in 1:nrecords
        if records[i] == 0
            sequence[i] = NOCALL_BASE
            quality[i] = NOCALL_QUALITY + PHRED_BASE
        else
            sequence[i] = CALL_BASES[1 + records[i] & 3]
            quality[i] = records[i] >> 2 + PHRED_BASE
        end
    end

    bcl.position += nrecords

    return (sequence, quality)
end

function seek(bcl::BCLFile, newposition::Int)
    seek(bcl.handle, BCL_HEADER_SIZE + sizeof(Uint8) * newposition)
    bcl.position = newposition
end


type BCLCollection
    datadir::String
    lane::Int
    tile::Int

    firstcycle::Int
    ncycle::Int
    ncluster::Int
    position::Int
    bclfiles::Array{BCLFile}

    function BCLCollection(datadir, lane, tile, firstcycle=-1, ncycle=-1)
        lanedir = "$datadir/BaseCalls/L00$lane"

        cycledirs = filter(x -> ismatch(r"^C[0-9]+\.1$", x), readdir(lanedir))
        cyclenums = map(x -> int(match(r"^C([0-9]+)\.1$", x).captures[1]), cycledirs)
        if firstcycle == -1
            firstcycle = minimum(cyclenums)
        end
        if ncycle == -1
            ncycle = maximum(cyclenums) - firstcycle + 1
        end

        bclfiles = [BCLFile("$lanedir/C$cycle.1/s_$(lane)_$(tile).bcl")
                    for cycle in firstcycle:(firstcycle+ncycle-1)]
        ncluster = bclfiles[1].ncluster

        if any(x -> x.ncluster != ncluster, bclfiles)
            error("Number of cluster is not consistent over the cycles.")
        end 

        return new(datadir, lane, tile, firstcycle, ncycle, ncluster, 0, bclfiles)
    end 
end

function read(bcl::BCLCollection, count::Int = -1) 
    if count < 0 
        nrecords = bcl.ncluster - bcl.position
    else
        nrecords = minimum((bcl.ncluster - bcl.position, count))
    end 

    sequences = Array(Uint8, (nrecords, bcl.ncycle))
    qualities = Array(Uint8, (nrecords, bcl.ncycle))

    for cycle in 1:bcl.ncycle
        seq, qual = read(bcl.bclfiles[cycle], nrecords)
        sequences[:, cycle] = seq
        qualities[:, cycle] = qual
    end

    sequences_t = permutedims(sequences, [2, 1])
    qualities_t = permutedims(qualities, [2, 1])

    sequences_ascii = [ascii(sequences_t[:, i]) for i in 1:nrecords]
    qualities_ascii = [ascii(qualities_t[:, i]) for i in 1:nrecords]

    return (sequences_ascii, qualities_ascii)
end

