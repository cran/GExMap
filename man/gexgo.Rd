\name{gexgo}
\Rdversion{1.1}
\alias{gexgo}
\title{
GO terms analysis of the tested gene list
}
\description{
Search and analyze the GO identifiers for each idetified gene from the tested gene list and produces 3 scored lists, one for each GO (Biological processes, Cellular component and Biological function. )
This function is called by the main function gexmap() if the variable isGO=TRUE (by default).
}
\usage{
gexgo(list.ENS, source, res)
}
\arguments{
  \item{list.ENS}{ List of all identified genes }
  \item{source}{ Source data folder where is located the go.Rdata file }
  \item{res}{ Results folder where will be saved all the data files producede by the function }
}
\details{
This function is included in the GExMap package and could not be used independently.
}
\value{
  The GO.results folder created contains 3 files for each GO type:\\
  \item{detail. .txt }{Complete list of genes for each GO ID sorted by decreasing score}
  \item{GO. .txt }{List of scored GO id with scores and relative proportion of each ID}
  \item{GO_ .pdf }{Pie chart of the GO ID which score is > 1 percent }
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
