% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCollection.R
\name{getCollection}
\alias{getCollection}
\title{Retrieve a Collection: Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats}
\usage{
getCollection(
  db = "refseq",
  organism,
  reference = TRUE,
  skip_bacteria = TRUE,
  release = NULL,
  assembly_type = "toplevel",
  analyse_genome = FALSE,
  remove_annotation_outliers = FALSE,
  gunzip = FALSE,
  path = file.path("_db_downloads", "collections"),
  mute_citation = FALSE
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{Organism selector id,
there are three options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{skip_bacteria}{Due to its enormous dataset size (> 700MB as of July 2023),
the bacterial summary file will not be loaded by default anymore. If users
wish to gain insights for the bacterial kingdom they needs to actively specify \code{skip_bacteria = FALSE}. When \code{skip_bacteria = FALSE} is set then the
bacterial summary file will be downloaded.}

\item{release}{a numeric, the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used. \code{release = 75} would for human would give the stable
GRCh37 release in ensembl. Value must be > 46, since ensembl did not structure their data
if the standard format before that.}

\item{assembly_type}{character, default c("primary_assembly", "toplevel"). Used for ensembl only,
specifies the genome assembly type. Searches for both primary and toplevel, and if both are found, uses the
first by order (so primary is prioritized by default).
The Primary assembly should usually be used if it exists.
The "primary assembly" contains all the top-level sequence regions,
excluding alternative haplotypes and patches.
If the primary assembly file is not present for a species
(only defined for standard model organisms),
that indicates that there were no haplotype/patch regions,
and in such cases, the 'toplevel file is used.
For more details see:
\href{https://grch37.ensembl.org/info/genome/genebuild/assembly.html}{ensembl tutorial}}

\item{analyse_genome}{logical, default FALSE. If TRUE, get general genome statistics like
gc content etc. For more details, see ?summary_genome}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input \code{annotation_file}?
If yes, then the initial \code{annotation_file} will be overwritten and the removed outlier lines will be stored at \code{\link{tempdir}}
for further exploration.}

\item{gunzip}{a logical, indicating whether or not files should be unzipped.}

\item{path}{a character string specifying the location (a folder) in which
the corresponding collection shall be stored. Default is
\code{path} = \code{file.path("_db_downloads","collections")}.}

\item{mute_citation}{logical, default FALSE, indicating whether citation message should be muted.}
}
\value{
File path to downloaded genome.
}
\description{
Main collection retrieval function for an organism of interest.
By specifying the scientific name of an organism of interest a collection consisting of
the genome file, proteome file, CDS file, RNA file, GFF file, Repeat Masker file, AssemblyStats
file of the organism of interest
can be downloaded and stored locally. Collections can be retrieved from
several databases. For full set of collection elements, see:
biomartr:::supported_biotypes(db)
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory relative to file type, if you get fasta genomes it will be
_ncbi_downloads/genomes'.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/genomes' folder and is accessible within the workspace,
no download process will be performed.
For other file types the same rule applies.
}
\examples{
\dontrun{
# download the collection of Homo sapiens from refseq
# and store the corresponding genome file in '_ncbi_downloads/collection'
 Hsap_collection <- getCollection( db       = "refseq",
             organism = "Homo sapiens",
             path = file.path("_db_downloads","collections"))
# download the collection of Homo sapiens from genbank
# and store the corresponding genome file in '_ncbi_downloads/collection'
 Hsap_collection <- getCollection( db       = "genbank",
             organism = "Homo sapiens",
             path = file.path("_db_downloads","collections"))
# download the collection of Homo sapiens from ensembl
# and store the corresponding genome file in '_ncbi_downloads/collection'
 Hsap_collection <- getCollection( db       = "ensembl",
             organism = "Homo sapiens",
             path = file.path("_db_downloads","collections"))
}

}
\seealso{
Other getBio: 
\code{\link{getBio}()},
\code{\link{getCDS}()},
\code{\link{getGFF}()},
\code{\link{getGenome}()},
\code{\link{getProteome}()},
\code{\link{getRNA}()}

Other collection: 
\code{\link{getCollectionSet}()}
}
\author{
Hajk-Georg Drost
}
\concept{collection}
\concept{getBio}
