\name{list}
\alias{list}
\docType{data}
\title{ Gene list used as an example by GExMap }
\description{
 Text file tab separated with two columns.
 The first columns is for the probes identifiers  and the second for the variation code (1 or -1).
 The title of the first column must be formated as follow: [type of mocroarray], [microarray name].
 The title of the second column must be "expression".
}
\usage{data(list)}
\format{
  A data frame with 3856 observations on the following variable.
  \describe{
    \item{\code{V2}}{a factor with levels \code{-1} \code{1} \code{expression}}
  }
}
\details{
  This dat set included in the GExMap package as an example.
}
\source{
http://gexmap.voila.net/index.html
}
\keyword{datasets}
