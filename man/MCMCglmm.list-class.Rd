\name{MCMCglmm.list-class}
\Rdversion{1.1}
\docType{class}
\alias{MCMCglmm.list-class}
\alias{initialize,MCMCglmm.list-method}
\alias{subsetMCMCglmm.list,MCMCglmm.list-method}

\title{Class "MCMCglmm.list"}
\description{
Class that stores MCMC runs
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MCMCglmm.list", b)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "MCMCglmm.list")}: ... }
    \item{subsetMCMCglmm.list}{\code{signature(x = "MCMCglmm.list")}: ... }
	 }
}

\author{
Lucia Spangenberg, Hugo Naya
}

\keyword{classes}
