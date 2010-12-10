\name{gex.mapping}
\alias{gex.mapping}
\title{ Function from the GExMap package which format the genome to be used by GExMap }
\description{
 Format the genome according to the scale and the genomic distribution of the tested gene list.
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
  This function is included in the GExMap package and should not be used independently.
}
\references{
http://gexmap.voila.net/index.html
Research Article JPB/Vol.2.1 051- 059(2009)/January 2009 doi:10.4172/jpb.1000060
GExMap: An Intuitive Visual Tool to Detect and Analyze Genomic Distribution in Microarray-generated Lists of Differentially Expressed Genes
Nicolas Cagnard, Carlo Lucchesi, Gilles Chiocchia 
}
\author{
Nicolas Cagnard <nicolas.cagnard@gmail.com>
}
