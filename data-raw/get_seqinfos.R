require(GenomeInfoDb)
require(purrr)

genomes <- c("hg38", "hg19", "hg18", "panTro4", "panTro3", "panTro2",
             "bosTau8", "bosTau7", "bosTau6", "canFam3", "canFam2",
             "canFam1", "musFur1", "mm10", "mm9", "mm8", "susScr3",
             "susScr2", "rn6", "rheMac3", "rheMac2", "galGal4", "galGal3",
             "gasAcu1", "danRer7", "apiMel2", "dm6", "dm3", "ce10", "ce6",
             "ce4", "ce2", "sacCer3", "sacCer2")
names(genomes) <- genomes

si <- map(genomes, ~ Seqinfo(genome = .x))
save(si, file = "data/si.RData", compress = "xz")
