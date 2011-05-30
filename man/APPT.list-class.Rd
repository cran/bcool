\name{APPT.list-class}
\Rdversion{1.1}
\docType{class}
\alias{APPT.list-class}
\alias{MCMCglmm.APPT.list,APPT.list,formula,character,numeric,numeric,list,logical,logical-method}
\alias{RGBcolor.APPT.list,APPT.list,numeric,logical-method}
\alias{[,APPT.list,logical,ANY,ANY-method}
\alias{[,APPT.list,numeric,ANY,ANY-method}
\alias{alignment,APPT.list-method}
\alias{anovaAPPT.list,APPT.list,formula-method}
\alias{bootMCMCglmm.APPT.list,APPT.list,formula,numeric,character,character,numeric,numeric,list,logical,logical-method}
\alias{cer.APPT.list,APPT.list,character-method}
\alias{columns,APPT.list-method}
\alias{initialize,APPT.list-method}
\alias{length,APPT.list-method}
\alias{levelsMCMCglmm,APPT.list-method}
\alias{multiMCMCglmm,APPT.list-method}
\alias{pheno,APPT.list-method}
\alias{properties,APPT.list-method}
\alias{summaryAPPT.list,APPT.list,character,character-method}
\alias{tree,APPT.list-method}

\title{Class "APPT.list"}
\description{
This class join all the elements needed for the detection of interesting columns in an aligment with labels using the MCMCglmm. It stores an alignment of a protein of interest in several organisms and a corresponding phylogenetic tree. The classes of interest of the organisms, e. g. pathogenicity, aerobiosis, etc.. In addition, amino acid properties, which will hopefully differ in the classes considered. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("APPT.list", alignment, pheno, properties, tree, columns, multiMCMCglmm.list, nowarns)}.

}
\section{Slots}{
  \describe{
    \item{\code{alignment}:}{Object of class \code{"matrix"} ~~ An alignment matrix, obtained from the alignment object from the seqinr package}
    \item{\code{pheno}:}{Object of class \code{"data.frame"} ~~ A data frame storing the phenotypic data, the classes considered}
    \item{\code{properties}:}{Object of class \code{"list"} ~~ The amino acid properties considered}
    \item{\code{tree}:}{Object of class \code{"phylo"} ~~ A tree of the class phylo (ape package)}
    \item{\code{columns}:}{Object of class \code{"numeric"} ~~ The indices of the columns of the alignment, which are considered for the calculations}
    \item{\code{multiMCMCglmm.list}:}{Object of class \code{"list"} ~~ A list of MCMCglmm runs}
  }
}
\section{Methods}{
  \describe{
    \item{MCMCglmm.APPT.list}{\code{signature(x = "APPT.list", class.var = "formula", random.eff = "character", nitt = "numeric", burnin = "numeric", prior = "list", scale = "logical", parallel = "logical")}: Loads the MCMCglmm list into the APPT.list object, in the corresponding multiMCMCglmm slot}
    \item{RGBcolor.APPT.list}{\code{signature(x = "APPT.list", which.prop = "numeric", write.score = "logical")}: Generates a RGB color code for the amino acid properties in the alignment, in order to be able to visually detect some tendencies and interesting columns}
    \item{[}{\code{signature(x = "APPT.list", i = "logical", j = "ANY", drop = "ANY")}: Subsets the APPT.list object }
    \item{[}{\code{signature(x = "APPT.list", i = "numeric", j = "ANY", drop = "ANY")}: Subsets the APPT.list object }
    \item{alignment}{\code{signature(value = "APPT.list")}: Retrieves the alignment slot of the APPT.list object}
    \item{anovaAPPT.list}{\code{signature(x = "APPT.list", class.var = "formula")}: Calculates for each column an ANOVA in order to detect interesting column in a first and fast approach. The significant columns will be further analyzed}
    \item{bootMCMCglmm.APPT.list}{\code{signature(x = "APPT.list", class.var = "formula", boot = "numeric", contrast = "character", random.eff = "character", nitt = "numeric", burnin = "numeric", prior = "list", scale = "logical", parallel = "logical")}: In order to assess the significance of a column, a bootstrap method has been implemented. Randomizing the observations and in each iteration recalculating the integrated score of the columns, gives an idea of the empirical random distribution of the score. Columns scores above the random distribution are considered significant.}
    \item{cer.APPT.list}{\code{signature(x = "APPT.list", class.var = "character")}: For each column the conditional entropy is calculated in order to scan potentially interesting columns in a fast manner. Columns with low conditional entropy values are interesting and will be firther analyzed.}
    \item{columns}{\code{signature(value = "APPT.list")}: Retrieves the column slot of the APPT.list object }
    \item{initialize}{\code{signature(.Object = "APPT.list")}:initializes the APPT.list object}
    \item{length}{\code{signature(x = "APPT.list")}: Retrieves the lenght of the APPT.list object}
    \item{levelsMCMCglmm}{\code{signature(value = "APPT.list")}: Retrieves the levels of the classes considered}
    \item{multiMCMCglmm}{\code{signature(value = "APPT.list")}: Retrieves the multiMCMCglmm slot of the APPT.list object }
    \item{pheno}{\code{signature(value = "APPT.list")}: Retrieves the pheno slot of the APPT.list object}
    \item{properties}{\code{signature(value = "APPT.list")}: Retrieves the properties slot of the APPT.list object}
    \item{summaryAPPT.list}{\code{signature(x = "APPT.list", contrast = "character", class.var = "character")}: Retrieves the summary of the whole run. All the relevant information is given in the summary function. The table of the scores for each property and each column, the effect sizes, the integrated schores, the frequency of the amino acids in each class, etc..}
    \item{tree}{\code{signature(value = "APPT.list")}: Retrieves the tree slot of the APPT.list object}
	 }
}

\author{Lucia Spangenberg \email{lucia@pasteur.edu.uy} and Hugo Naya \email{naya@pasteur.edu.uy}}

\examples{
data("rpoS")
data("aaindex")

# select the properties to consider from "aaindex"
prop<-c("LEVM780101","LEVM780102","CHOC760102","JOND920102","KYTJ820101")

# instantiate the APPT.list object
myAPPT.list<-new("APPT.list",alignment=rpoSalign,pheno=labels,properties=lapply(aaindex[prop],function(x) x$I),tree=tree7)
}

\keyword{classes}
