% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMetaGenomeAnnotations.R
\name{getMetaGenomeAnnotations}
\alias{getMetaGenomeAnnotations}
\title{Retrieve annotation *.gff files for metagenomes from NCBI Genbank}
\usage{
getMetaGenomeAnnotations(
  name,
  path = file.path("_ncbi_downloads", "metagenome", "annotations"),
  metagenomes.members = dplyr::filter(getMetaGenomeSummary(), organism_name == name &
    total_gene_count > 0)
)
}
\arguments{
\item{name}{metagenome name retrieved by \code{\link{listMetaGenomes}}.}

\item{path}{a character string specifying the location (a folder)
in which the corresponding metagenome annotations shall be stored.
Default is
\code{path} = \code{file.path("_ncbi_downloads","metagenome","annotations")}.}

\item{metagenomes.members}{a tibble with selected assemblies,
default: dplyr::filter(getMetaGenomeSummary(), organism_name == name & total_gene_count > 0)).
This is different to getMetaGenome since it requires a gff to exist, most genbank
assemblies are .gbff files only, which are usually not useful.}
}
\description{
Retrieve available annotation *.gff files for metagenomes
from NCBI Genbank. NCBI Genbank allows users
to download entire metagenomes and their annotations of several metagenome
projects. This function downloads available metagenomes that can then be
downloaded via \code{\link{getMetaGenomes}}.
}
\examples{
\dontrun{
# Frist, retrieve a list of available metagenomes
listMetaGenomes()

# Now, retrieve the 'human gut metagenome'
getMetaGenomeAnnotations(name = "human gut metagenome")
}
}
\seealso{
\code{\link{getMetaGenomes}}, \code{\link{listMetaGenomes}},
\code{\link{getGFF}}
}
\author{
Hajk-Georg Drost
}
