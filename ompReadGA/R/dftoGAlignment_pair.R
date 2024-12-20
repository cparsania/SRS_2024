#' Extract Paired-End Alignment Details from BAM File
#'
#' Extracts paired-end alignment details from a BAM file and returns them as a
#' GAlignments object. Metadata columns include sequences, quality scores, mapping quality, mate sequence names,
#' mate positions, ReadID, mapping quality and template length. 
#'
#'
#' @param bam_file Path to the BAM file.
#' @param n_threads_to_use Number of threads to use for processing.
#' @return A GAlignments object containing the paired-end alignment details.
#' @examples
#' # Example with a BAM file
#' bam_path <- "example_paired.bam"
#' ompReadGAlignmentsPaired(bam_path, n_threads_to_use = 4)
#' @importFrom GenomicAlignments GAlignments
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings DNAStringSet
#' @export

ompReadGAlignmentsPaired <- function(bam_file, n_threads_to_use = 1L) {
  alignment_data <- extract_alignment_details_paired(bam_file, n_threads_to_use)
  galignments <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(alignment_data$seqnames),
    pos = alignment_data$start,
    cigar = alignment_data$cigar,
    strand = S4Vectors::Rle(alignment_data$strand),
    qwidth = nchar(alignment_data$sequences),
    seqinfo = NULL
  )

  # Add mate information as metadata
  mcols(galignments)$sequences <- Biostrings::DNAStringSet(alignment_data$sequences)
  mcols(galignments)$quality_scores <- alignment_data$quality_scores
  mcols(galignments)$mapq <- alignment_data$mapq
  mcols(galignments)$mate_seqnames <- alignment_data$mate_seqnames
  mcols(galignments)$mate_positions <- alignment_data$mate_positions
  mcols(galignments)$read_id <- alignment_data$read_id
  mcols(galignments)$mate_mapq <- alignment_data$mate_mapq
  mcols(galignments)$template_length <- alignment_data$template_length

  return(galignments)
}
