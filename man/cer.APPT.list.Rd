\name{cer.APPT.list}
\alias{cer.APPT.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Conditional entropy calculation for each column in the alignment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the reduction of the entropy in the aligment, when the class of the organisms are known. The more informative the label the smaller the conditional entropy. Columns with small values are potentially interesting and are further analyzed.
}
\usage{
cer.APPT.list(x, class.var, which.columns, maxngaps, nummin, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
APPT.list object}
  \item{class.var}{
character, class considered}
  \item{which.columns}{
numeric
}
  \item{maxngaps}{
numeric, default 0
}
  \item{nummin}{
numeric, default 2
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~

}
}

\value{
numeric vector, containing the conditional entropy values for all columns considered
}
\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}

\examples{
data("rpoS")
data("aaindex")

# select the properties to consider from "aaindex"
prop<-c("JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)

# CER fast approach
my.cer<-cer.APPT.list(myAPPT.list,class.var="pathogenicity")
}
