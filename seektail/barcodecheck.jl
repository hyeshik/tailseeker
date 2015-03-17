

include("bclfile.jl")

type BarcodeChecker
    datadir::String
    lane::Int
    tile::Int

    indexrange::UnitRange{Int}
    barcodes::Array{String}
    bclreader::BCLCollection

    function BarcodeChecker(datadir::String, lane::Int, tile::Int,
                            barcodes::Array{String},
                            indexrange::UnitRange{Int}=52:57)

        bclreader = BCLCollection(datadir, lane, tile)
        barcodes = Array(String, 0)

        return new(datadir, lane, tile, indexrange, barcodes, bclreader)
    end
end



barcodes = convert(Array{String}, split("GCCAAT GCCAAT CTTGTA GATCAG TAGCTT ATCACG CGATGT TTAGGC TGACCA ACAGTG CAGATC ACTTGA"))

bchecker = BarcodeChecker("sample-miseq/Data/Intensities", 1, 1101, barcodes, 52:57)
println(bchecker)



#
#        bclreader = BCLCollection(datadir, lane, tile)
##bcl = BCLCollection("sample-hiseq/Data/Intensities", 1, 2101, 1, 51)
#        return new(datadir, lane, tile, blocksize, bowtie2command, bowtiethreads, phixindexpath)
#    end
#
#    function Base.start(reader::TAILseqPass1Reader)
#    end
#
#    function Base.next(reader::TAILseqPass1Reader, state)
#    end
#
#    function Base.done(reader::TAILseqPass1Reader, state)
#    end
#end
#
#
#reader = TAILseqPass1Reader("sample-miseq/Data/Intensities", 1, 2101)
#println(reader)
#
#
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
