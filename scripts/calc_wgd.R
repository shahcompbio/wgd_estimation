library(argparse)
library(GenomicRanges)

# function from https://github.com/gerstung-lab/PCAWG-11/blob/master/code/PCAWG-functions.R
averageHom <- function(bb){
    sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

# function from https://github.com/gerstung-lab/PCAWG-11/blob/master/code/PCAWG-functions.R
# heuristic developed in PCAWG-11
#.classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy

# heuristic from Andrew
.classWgd <- function(ploidy, hom) 2.8 - (1.1/0.85)*hom <= ploidy


readCnTable <- function(cn_path, clonal_freq) {
    cn = read.table(cn_path, header=1)
    bb <- GRanges( # bb file formated for MutationTimeR
      seqnames = Rle( as.character(cn$chromosome) ),
      ranges = IRanges(
        as.numeric(cn$start),
        end = as.numeric(cn$end)
      ),
      major_cn = as.integer(cn$major_1),
      minor_cn = as.integer(cn$minor_1),
      clonal_frequency = rep(clonal_freq, dim(cn)[1]), 
    )
    return(bb)
}

get_args <- function() {
    p <- ArgumentParser(description = "Calculate WGD status from copy number, purity and ploidy")

    p$add_argument("cn", help = "CN tsv file with columns: chromosome, start, end, major_1, minor_1")
    p$add_argument("cf", help = "Clonal frequency (purity)")
    p$add_argument("ploidy", help = "Ploidy")

    p$add_argument("wgd", help = "output wgd tsv")

    return(p$parse_args())
}


main <- function() {
    argv <- get_args()

    clonal_freq <- as.numeric(argv$cf) # purity
    ploidy <- as.numeric(argv$ploidy) # ploidy
    print("[LOG] ploidy") 
    print(ploidy)
    bb <- readCnTable(argv$cn, clonal_freq) # cn_path
    out_path <- argv$wgd
    
    # get WGD status
    hom <- averageHom(bb)
    print("[LOG] hom")
    print(hom) ##@##
    isWgd <- .classWgd(ploidy, hom)
    print("[LOG] isWgd")
    print(isWgd)

    # save output
    LOH <- c(hom)
    Ploidy <- c(ploidy)
    WGD <- c(isWgd)

    df <- data.frame(LOH, Ploidy, WGD)
    write.table(df, out_path, row.names=F, quote=F, sep="\t")
}

main()
