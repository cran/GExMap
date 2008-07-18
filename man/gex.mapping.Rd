\name{gex.mapping}
\alias{gex.mapping}
\title{ Function from the GExMap package which format the genome to be used by GExMap }
\description{
 Format the genome accordinf to the scale and the genomic distribution of the user gene list.
}
\usage{
gex.mapping(data, list.chr, scale, I, list.ENS)
}
\arguments{
  \item{data}{ Matrix of the data.data file }
  \item{list.chr}{ Data about the analized chromosome }
  \item{scale}{ Scale to format the genome. }
  \item{I}{ This function is performed chromosome by chromosome, I is the name of the analized chromosome }
  \item{list.ENS}{ List of ENSEMBL genes produced by the \code{load.corr()} function. }
}
\details{
  This function is included in the GExMap package and could not be used independently.
}
\value{
}
\references{  }
\author{ Nicolas Cagnard <nicolas.cagnard@inserm.fr> }
\note{
}
\seealso{GExMap}
\keyword{ misc }

