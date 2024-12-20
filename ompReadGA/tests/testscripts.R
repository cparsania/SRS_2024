library(IRanges)
library(ompBAM)

# Updating the modifications
pkg_path <- "~/Desktop/Uni /Research/Centenary /Documents/ompReadGA"
setwd(pkg_path)
devtools::clean_dll()
devtools::document(pkg = pkg_path)
devtools::load_all(recompile = TRUE)
devtools::check()

library(ompReadGA)
setwd("~/Desktop/Uni /Research/Centenary /Documents")

result_paired <- ompReadGAlignmentsPaired("subset.bam", 2)
result_normal <- ompReadGAlignments(bam_file = ompBAM::example_BAM("Unsorted"),n_threads_to_use =  3)

?ompReadGAlignments
devtools::check()
