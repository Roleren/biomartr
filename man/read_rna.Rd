% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_rna.R
\name{read_rna}
\alias{read_rna}
\title{Import RNA as Biostrings or data.table object}
\usage{
read_rna(file, format = "fasta", obj.type = "Biostrings", ...)
}
\arguments{
\item{file}{a character string specifying the path to the file
storing the RNA.}

\item{format}{a character string specifying the file format used to store the
genome, e.g. \code{format = "fasta"} (default) or \code{format = "gbk"}.}

\item{obj.type}{a character string specifying the object stype in which the
genomic sequence shall be represented.
Either as \code{obj.type = "Biostrings"} (default) or
as \code{obj.type = "data.table"}.}

\item{...}{additional arguments that are used by
\code{\link[seqinr]{read.fasta}}.}
}
\value{
A data.table storing the gene id in the first column and the
corresponding sequence as string in the second column.
}
\description{
This function reads an organism specific RNA stored in a
defined file format.
}
\details{
This function takes a string specifying the path to the RNA file
of interest as first argument. It is possible to read in different proteome
file standards such as \emph{fasta} or \emph{genebank}.
}
\seealso{
Other rna: 
\code{\link{getRNA}()},
\code{\link{getRNASet}()}

Other readers: 
\code{\link{read_cds}()},
\code{\link{read_genome}()},
\code{\link{read_gff}()},
\code{\link{read_proteome}()}
}
\author{
Hajk-Georg Drost
}
\concept{readers}
\concept{rna}
