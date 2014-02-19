\name{anovaAPPT.list}
\alias{anovaAPPT.list}
\title{
ANOVA on each column
}
\description{
ANOVA is applied as a fast scanning of the columns. For each column an ANOVA is determined, and a p-value for each is calculated. Significant columns (according to a determined significance level) are kept for further analysis.
}
\usage{
anovaAPPT.list(x, class.var, which.columns, maxngaps, nummin, ...)
}
\arguments{
  \item{x}{
  	APPT.list object
}
  \item{class.var}{
character, the classes considered
}
  \item{which.columns}{
numeric, which columns are considered
}
  \item{maxngaps}{
numeric, maximal number of gaps in order to consider a column for calculations. If gaps > maxngaps the column is not considered.
}
  \item{nummin}{
numeric, minimal number of properties considered
}
  \item{\dots}{
}
}
\value{
 \item{cuantos}{how many of the properties considered are significant at a level of 0.01}
 \item{median }{the median (log-transformed) of the result of the ANOVA for each property}
 \item{mean }{the mean (log-transformed) of the result of the ANOVA for each property}

}

\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}


\examples{
data("rpoS")

env <- new.env()
utils::data("aaindex", package = "seqinr", envir = env)
aaindex <- env$aaindex


# select the properties to consider from "aaindex"
prop<-c("JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,
properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)

# ANOVA fast approach
my.anova<-anovaAPPT.list(myAPPT.list,class.var=~pathogenicity,
which.columns=NULL,nummin=2,maxgaps=10)
# sum(my.anova$hm.signif==2,na.rm=TRUE)
# which(my.anova$hm.signif==2)

}



