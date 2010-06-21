\name{gexload.corr}
\Rdversion{1.1}
\alias{gexload.corr}
\title{
Function from the GExMap package which annotates identifier list
}
\description{
Creates a list of ENSEMBL id from the user gene list which contains microarray id
}
\usage{
gexload.corr(data, list, scale, source, res)
}
\arguments{
  \item{data}{ The ENSEMBL genome reference data. }
  \item{list}{ The user gene list to be tested and which contains all id which are to be annotated. }
  \item{scale}{ Scale to format the genome. }
  \item{source}{ The source folder where would be placed all the Rdata files. }
  \item{res}{ The results folder where would be saved all the resulint files. }
}
\details{
This function produces a matrix with annotations and genomic data for all the genome.
}
\value{
The columns of the matrix produced by this function are:
	\item{chr }{ Chromosome names. }
	\item{m }{ Middle of the cytobands. }
	\item{cytoband }{ Name of the cytobands. }
	\item{name }{ Gene names. }
	\item{ensembl }{ ensembl IDs. }
	\item{GO }{ GO terms. }
	\item{expression }{ Regulation type of the gene: 1 for up regulation and -1 for down regulation. }
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
