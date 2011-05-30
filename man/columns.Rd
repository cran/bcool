\name{columns}
\alias{columns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Retrieves the column slot of the APPT.list object 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Retrieves the column slot of the object APPT.list, so that the user doesn0x00B4t need to access the object directly.

}
\usage{
columns(value)
}

\arguments{
  \item{value}{
%%     ~~Describe \code{value} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
numeric vector
}

\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}

\examples{
data("rpoS")
data("aaindex")

# select the properties to consider from "aaindex"
prop<-c("JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)
columns(myAPPT.list)
}
