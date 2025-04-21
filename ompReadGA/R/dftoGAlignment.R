#' Convert BAM Alignment Data to GAlignments
#'
#' This function extracts alignment details from a BAM file and returns them as a
#' GAlignments object. Metadata columns include sequences, quality and mapping quality. 
#'
#' @param bam_file Path to the BAM file.
#' @param n_threads_to_use Number of threads to use for processing.
#' @return A GAlignments object containing the alignment details.
#' @examples 
#' # Example Usage 
#' ompReadGAlignments("path/to/file.bam", n_threads_to_use = 4)
#' @importFrom GenomicAlignments GAlignments
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings DNAStringSet
#' @export

ompReadGA <- function(bam_file, n_threads_to_use = 1L, fields_to_return = c("seq", "qual", "mapq")) {
  alignment_data <- extract_alignment_details(bam_file, n_threads_to_use)
  
  alignments <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(alignment_data$seqnames),
    pos = alignment_data$start,
    cigar = alignment_data$cigar,
    strand = S4Vectors::Rle(alignment_data$strand)
  )
  
  if ("seq" %in% fields_to_return) {
    sequences <- Biostrings::DNAStringSet(alignment_data$sequence)
    S4Vectors::mcols(alignments)$seq <- sequences
  }
  
  if ("qual" %in% fields_to_return) {
    S4Vectors::mcols(alignments)$qual <- alignment_data$quality_scores
  }
  
  if ("mapq" %in% fields_to_return) {
    S4Vectors::mcols(alignments)$mapq <- alignment_data$mapq
  }
  
  return(alignments)
}
