\name{gexgo}
\alias{gexgo}
\title{ Perform a GO analysis of the tested gene list }
\description{
  Search and analyze the GO identifiers for each idetified gene from the tested gene list and produces 3 scored lists, one for each GO (Biological processes, Cellular component and Biological function. )
  This function is used or not by GExMap, called by the main function gexmap() if the variable isGO=TRUE (by default).
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
  A GO.results folder is created which contains 3 files for each GO:\\
  \item{detail. .txt }{Complete list of genes for each GO ID sorted by decreasing score}
  \item{GO. .txt }{List of scored GO id with scores and relative proportion of each ID}
  \item{GO_ .pdf }{Pie chart of the GO ID which score is > 1 percent }
}
\references{ }
\author{ Nicolas Cagnard <nicolas.cagnard@inserm.fr> }
\note{
}
\seealso{ GExMap }
\keyword{ misc }

