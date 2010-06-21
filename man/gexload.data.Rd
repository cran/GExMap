\name{gexload.data}
\Rdversion{1.1}
\alias{gexload.data}
\title{
Function from the GExMap package which upload data files to be used by GExMap
}
\description{
Automatically uploads (or ask to the user) the data.Rdata file and the *.Rdata file specific to the user identifier type.
}
\usage{
gexload.data(source, genome)
}
\arguments{
  \item{source}{ The source folder where would be placed all the Rdata files. }
  \item{genome}{ ENSEMBL genome whith all genomic informations }
}
\details{
Upload the genome Rdata file.
}
\value{
	\item{ data }{ The uploaded genomic data file. }
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
