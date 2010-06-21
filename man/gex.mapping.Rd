\name{gex.mapping}
\Rdversion{1.1}
\alias{gex.mapping}
\title{
Function from the GExMap package which format the genome to be used by GExMap
}
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
This function is included its format the genomic annotaions into a genomic map data matrix.
}
\value{
The columns of the genomic map data matrix produced by this function are: 
	\item{chr}{ Name of the current chromosome. }
	\item{scale}{ Number for each unit (1, 2, 3, ...) in the genomic order }
	\item{ensembl}{ Sum of the ensemble genes for this genomic unit. }
	\item{list}{ Sum of the genes of the tested list for this genomic unit. }
}
\references{
http://gexmap.voila.net/index.html
Research Article JPB/Vol.2.1 051- 059(2009)/January 2009 doi:10.4172/jpb.1000060
GExMap: An Intuitive Visual Tool to Detect and Analyze Genomic Distribution in Microarray-generated Lists of Differentially Expressed Genes
Nicolas Cagnard, Carlo Lucchesi, Gilles Chiocchia 
}
\author{
Nicolas Cagnard <nicolas.cagnard@yahoo.fr>
}

