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

ompReadGAPaired <- function(bam_file, n_threads_to_use = 1L, 
                            fields_to_return = c("sequences", "quality_scores", "mapq", 
                                                 "mate_seqnames", "mate_positions", 
                                                 "read_id", "mate_mapq", "template_length")) {
  alignment_data <- extract_alignment_details_paired(bam_file, n_threads_to_use)
  
  galignments <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(alignment_data$seqnames),
    pos = alignment_data$start,
    cigar = alignment_data$cigar,
    strand = S4Vectors::Rle(alignment_data$strand),
    qwidth = nchar(alignment_data$sequences),
    seqinfo = NULL
  )
  
  # Add metadata based on requested fields
  if ("sequences" %in% fields_to_return) {
    mcols(galignments)$sequences <- Biostrings::DNAStringSet(alignment_data$sequences)
  }
  if ("quality_scores" %in% fields_to_return) {
    mcols(galignments)$quality_scores <- alignment_data$quality_scores
  }
  if ("mapq" %in% fields_to_return) {
    mcols(galignments)$mapq <- alignment_data$mapq
  }
  if ("mate_seqnames" %in% fields_to_return) {
    mcols(galignments)$mate_seqnames <- alignment_data$mate_seqnames
  }
  if ("mate_positions" %in% fields_to_return) {
    mcols(galignments)$mate_positions <- alignment_data$mate_positions
  }
  if ("read_id" %in% fields_to_return) {
    mcols(galignments)$read_id <- alignment_data$read_id
  }
  if ("mate_mapq" %in% fields_to_return) {
    mcols(galignments)$mate_mapq <- alignment_data$mate_mapq
  }
  if ("template_length" %in% fields_to_return) {
    mcols(galignments)$template_length <- alignment_data$template_length
  }
  
  return(galignments)
}
