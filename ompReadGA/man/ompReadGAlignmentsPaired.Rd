% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dftoGAlignment_pair.R
\name{ompReadGAlignmentsPaired}
\alias{ompReadGAlignmentsPaired}
\title{Extract Paired-End Alignment Details from BAM File}
\usage{
ompReadGAlignmentsPaired(bam_file, n_threads_to_use = 1L)
}
\arguments{
\item{bam_file}{Path to the BAM file.}

\item{n_threads_to_use}{Number of threads to use for processing.}
}
\value{
A GAlignments object containing the paired-end alignment details.
}
\description{
Extracts paired-end alignment details from a BAM file and returns them as a
GAlignments object. Metadata columns include sequences, quality scores, mapping quality, mate sequence names,
mate positions, ReadID, mapping quality and template length.
}
\examples{
# Example with a BAM file
bam_path <- "example_paired.bam"
ompReadGAlignmentsPaired(bam_path, n_threads_to_use = 4)
}
