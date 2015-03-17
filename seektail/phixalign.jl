importall Base

include("bclfile.jl")
include("ciffile.jl")


type PhiXSpotsReader
    datadir::String
    lane::Int
    tile::Int

    blocksize::Int
    bowtie2command::String
    bowtiethreads::Int
    phixindexpath::String

    read5start::Int
    read5ncycles::Int
    read3start::Int
    read3ncycles::Int

    function PhiXSpotsReader(datadir::String, lane::Int, tile::Int, blocksize::Int=500000,
                                bowtie2command::String="bowtie2",
                                bowtiethreads::Int=8, phixindexpath::String="phix/phix",
                                read5start::Int=1, read5ncycles::Int=51,
                                read3start::Int=58, read3ncycles::Int=251)

        return new(datadir, lane, tile, blocksize, bowtie2command, bowtiethreads, phixindexpath,
                   read5start, read5ncycles, read3start, read3ncycles)
    end
end


type PhiXSpotsReaderState
    bclreader_read5::BCLCollection
    bclreader_read3::BCLCollection
    cifreader_read3::CIFCollection

    function PhiXSpotsReaderState(reader::PhiXSpotsReader)
        bclreader_read5 = BCLCollection(reader.datadir, reader.lane, reader.tile,
                                        reader.read5start, reader.read5ncycles)
        bclreader_read3 = BCLCollection(reader.datadir, reader.lane, reader.tile,
                                        reader.read3start, reader.read3ncycles)
        cifreader_read3 = CIFCollection(reader.datadir, reader.lane, reader.tile,
                                        reader.read3start, reader.read3ncycles)

        if cifreader_read3.ncluster != bclreader_read5.ncluster ||
                cifreader_read3.ncluster != bclreader_read3.ncluster
            error("Some of CIF and BCL files does not have consistent number of clusters in lane $(reader.lane), tile $(reader.tile).")
        end

        return new(bclreader_read5, bclreader_read3, cifreader_read3)
    end
end

function Base.start(reader::PhiXSpotsReader)
    return PhiXSpotsReaderState(reader)
end

function Base.next(reader::PhiXSpotsReader, state::PhiXSpotsReaderState)
    println("next")
    return ("Yay! $(state.cifreader_read3.ncluster)", state)
end

function Base.done(reader::PhiXSpotsReader, state::PhiXSpotsReaderState)
    return is_eof(state.bclreader_read5)
end


reader = PhiXSpotsReader("sample-miseq/Data/Intensities", 1, 2101)
for i in reader
    println(i)
end


#const bowtie2command = "bowtie2"
#const phixindexpath = "phix/phix"
#const bclblocksize = 500000
#const bowtiethreads = 8
#
#bcl = BCLCollection("sample-hiseq/Data/Intensities", 1, 2101, 1, 51)
#is_control = zeros(Bool, bcl.ncluster)
#
#btout, btin, btproc = readandwrite(`$bowtie2command -p $bowtiethreads -f --quiet --end-to-end -k 1 --fast --no-head --no-sq --reorder -x $phixindexpath -U -`)
#
#@async begin
#    n = 1::Int64
#
#    while true
#        sequences, qualities = read(bcl, bclblocksize)
#        length(sequences) > 0 || break
#
#        println(">> Feeding $(length(sequences)) sequences.")
#        for seq in sequences
#            write(btin, ">$n\n$seq\n")
#            n += 1
#        end
#    end
#
#    close(btin)
#end
#
#
#rn = 1::Int64
#ncontrols = 0::Int64
#while true
#    line = readline(btout)
#    line != "" || break
#
#    seqno, flag = map(parseint, split(chomp(line), "\t")[1:2])
#    seqno == rn || error("Aligner output does not align to input one by one.")
#
#    if flag & 4 == 0 # aligned?
#        is_control[rn] = true
#        ncontrols += 1
#    end
#
#    rn += 1
#end
#
#println("PhiX: $ncontrols / $(bcl.ncluster)")
#
