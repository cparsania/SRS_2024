# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

extract_alignment_details <- function(bam_file, n_threads_to_use = 1L) {
    .Call(`_ompReadGA_extract_alignment_details`, bam_file, n_threads_to_use)
}

extract_alignment_details_paired <- function(bam_file, n_threads_to_use = 1L) {
    .Call(`_ompReadGA_extract_alignment_details_paired`, bam_file, n_threads_to_use)
}

