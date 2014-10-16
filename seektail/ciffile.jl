#
# Copyright (c) Institute for Basic Science, Korea
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

const NUM_CHANNELS = 4
const CIF_HEADER_SIZE = 13


type CIFFile
    path::String

    handle::IOStream
    firstcycle::Uint16
    ncluster::Uint32
    position::Uint32

    function CIFFile(filename::String)
        handle = open(filename)

        # Check the magic number in the beginning of the file
        filemagic = readbytes(handle, 3)
        if filemagic != b"CIF"
            error("Not a CIF file.")
        end

        # Read values in header
        version = read(handle, Uint8)
        datasize = read(handle, Uint8)
        firstcycle = ltoh(read(handle, Uint16))
        ncycle = ltoh(read(handle, Uint16))
        ncluster = ltoh(read(handle, Uint32))

        if version != 1
            error("Unsupport CIF version: $version.")
        end

        if datasize != sizeof(Int16)
            error("Unsupported data size ($datasize) in CIF intensities.")
        end

        if ncycle != 1
            error("Expected a single cycle CIF.")
        end

        return new(filename, handle, firstcycle, ncluster, 0)
    end
end


function read(cif::CIFFile, count::Int = -1)
    if count < 0
        nrecords = cif.ncluster - cif.position
    else
        nrecords = minimum((cif.ncluster - cif.position, count))
    end

    records = ltoh(read(cif.handle, Int16, (NUM_CHANNELS, int(nrecords))))
    cif.position += nrecords
    return records
end


function seek(cif::CIFFile, newposition::Int)
    seek(cif.handle, CIF_HEADER_SIZE + sizeof(Int16) * NUM_CHANNELS * newposition)
    cif.position = newposition
end


type CIFCollection
    datadir::String
    lane::Int
    tile::Int

    ncycle::Int
    ncluster::Int
    position::Int
    ciffiles::Array{CIFFile}

    function CIFCollection(datadir, lane, tile)
        lanedir = "$datadir/L00$lane"

        cycledirs = filter(x -> ismatch(r"^C[0-9]+\.1$", x), readdir(lanedir))
        ncycle = maximum(map(x -> int(match(r"^C([0-9]+)\.1$", x).captures[1]), cycledirs))
        ciffiles = [CIFFile("$lanedir/C$cycle.1/s_$(lane)_$(tile).cif") for cycle in 1:ncycle]
        ncluster = ciffiles[1].ncluster

        if any(x -> x.ncluster != ncluster, ciffiles)
            error("Number of cluster is not consistent over the cycles.")
        end

        return new(datadir, lane, tile, ncycle, ncluster, 0, ciffiles)
    end
end

function read(cif::CIFCollection, count::Int = -1)
    if count < 0
        nrecords = cif.ncluster - cif.position
    else
        nrecords = minimum((cif.ncluster - cif.position, count))
    end

    intensities = Array(Int16, (NUM_CHANNELS, cif.ncycle, nrecords))
    for cycle in 1:cif.ncycle
        intensities[:, cycle, :] = read(cif.ciffiles[cycle], nrecords)
    end

    return intensities
end
