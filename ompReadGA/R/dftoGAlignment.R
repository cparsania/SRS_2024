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

ompReadGAlignments <- function(bam_file, n_threads_to_use = 1L) {
  alignment_data <- extract_alignment_details(bam_file, n_threads_to_use)
  sequences <- Biostrings::DNAStringSet(alignment_data$sequence)
  quality <- alignment_data$quality_scores
  alignments <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(alignment_data$seqnames),
    pos = alignment_data$start,
    cigar = alignment_data$cigar,
    strand = S4Vectors::Rle(alignment_data$strand)
  )

  S4Vectors::mcols(alignments)$seq <- sequences
  S4Vectors::mcols(alignments)$qual <- quality
  S4Vectors::mcols(alignments)$mapq <- alignment_data$mapq

  return(alignments)
}
