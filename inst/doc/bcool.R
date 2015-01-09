##
##
## Copyright (c) 2011 Hugo Naya and Lucia Spangenberg
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 


#############################################################################
## setOldClass("phylo")
#############################################################################

setOldClass("phylo","pyhlo")


#############################################################################
## setClass("MCMCglmm.list")
#############################################################################

setClass("MCMCglmm.list",representation("list"),contains="list",prototype(list()))

setMethod("initialize","MCMCglmm.list",function(.Object,b) {
		if (sum(unlist(lapply(b,function(x)return(is(x)=="MCMCglmm"))))==length(b)) {
			.Object@.Data<-b
		} else {
			.Object@.Data<-b[unlist(lapply(b,function(x)return(is(x)=="MCMCglmm")))]
			warning(paste("length of MCMCglmm.list:",sum(unlist(lapply(b,function(x)return(is(x)=="MCMCglmm")))),", length of list:",length(b)))
		}	
	.Object
}
)



#############################################################################
## setClass("APPT.list")	Alignment Phenotype Properties Tree
#############################################################################

setClass("APPT.list",representation(alignment="matrix",pheno="data.frame",properties="list",tree="phylo",columns="numeric",multiMCMCglmm.list="list"))

setMethod("initialize","APPT.list",function(.Object,alignment,pheno,properties,tree,columns=vector("numeric"),multiMCMCglmm.list=list(),nowarns=FALSE) {
	
	if (length(tree$tip.label)==dim(alignment)[1]&length(tree$tip.label)==dim(pheno)[1]) {
		.Object@alignment<-alignment
		.Object@pheno<-pheno
		.Object@properties<-properties
		.Object@tree<-tree	
		if (length(columns)<1) {.Object@columns<-1:dim(alignment)[2]} else {.Object@columns<-columns}							
		.Object@multiMCMCglmm.list<-multiMCMCglmm.list

	}  else {
		if (nowarns==FALSE) warning (paste("Species in alignment:",dim(alignment)[1],", phenotype table:",dim(pheno)[1],"or tree:",length(tree$tip.label)," differs\nEmpty object!!"))	
		
	}	
	.Object
}
)



#############################################################################
## setGeneric("alignment")
#############################################################################

setGeneric("alignment",function(value)standardGeneric("alignment"))

setMethod("alignment","APPT.list",function(value){
		x<-value@alignment
		x
}
)


#############################################################################
## setGeneric("pheno")
#############################################################################

setGeneric("pheno",function(value)standardGeneric("pheno"))

setMethod("pheno","APPT.list",function(value){
		x<-value@pheno
		x
}
)

#############################################################################
## setGeneric("properties")
#############################################################################

setGeneric("properties",function(value)standardGeneric("properties"))

setMethod("properties","APPT.list",function(value){
		x<-value@properties
		x
}
)

#############################################################################
## setGeneric("tree")
#############################################################################

setGeneric("tree",function(value)standardGeneric("tree"))

setMethod("tree","APPT.list",function(value){
		x<-value@tree
		x
}
)

#############################################################################
## setGeneric("columns")
#############################################################################

setGeneric("columns",function(value)standardGeneric("columns"))

setMethod("columns","APPT.list",function(value){
		x<-value@columns
		x
}
)

#############################################################################
## setGeneric("multiMCMCglmm")
#############################################################################

setGeneric("multiMCMCglmm",function(value)standardGeneric("multiMCMCglmm"))

setMethod("multiMCMCglmm","APPT.list",function(value){
		x<-value@multiMCMCglmm.list
		x
}
)

#############################################################################
## setGeneric("levelsMCMCglmm")
#############################################################################

setGeneric("levelsMCMCglmm",function(value)standardGeneric("levelsMCMCglmm"))

setMethod("levelsMCMCglmm","APPT.list",function(value){
		x<-colnames(multiMCMCglmm(value)[[1]][[1]]$Sol)
		x
}
)



#############################################################################
## setMethod("length")
#############################################################################

setMethod("length", c("APPT.list"), function(x){
	if (length(tree(x)$tip.label)==dim(pheno(x))[1]&&length(tree(x)$tip.label)==dim(alignment(x))[1]) {
		length(tree(x)$tip.label)	
	} else {
		print(paste("Error: species in alignment:",dim(alignment)[1],", phenotype table:",dim(pheno)[1],"or tree:",length(tree$tip.label)," differs"))
	}
	}
)

#############################################################################
## setGeneric("subsetMCMCglmm.list")
#############################################################################

setGeneric("subsetMCMCglmm.list",function(x,value)standardGeneric("subsetMCMCglmm.list"))

setMethod("subsetMCMCglmm.list","MCMCglmm.list",function(x,value){
		x@.Data<-x@.Data[value]
		x
}
)



#############################################################################
## setMethod("[")
## WARNING this function subset only columns and multiMCMCglmm.list slots  
#############################################################################

setMethod("[", c("APPT.list", "numeric","ANY", "ANY"),
          function(x, i,..., drop=TRUE)
{
	new.multiMCMCglmm<-lapply(multiMCMCglmm(x),function(z) return(z[i]))	
   	y<-new("APPT.list",alignment=alignment(x),pheno=pheno(x),properties=properties(x),tree=tree(x),columns=columns(x)[i],multiMCMCglmm.list=new.multiMCMCglmm)
	y	
})

setMethod("[", c("APPT.list", "logical","ANY", "ANY"),
          function(x, i,..., drop=TRUE)
{
	new.multiMCMCglmm<-lapply(multiMCMCglmm(x),function(z) return(z[i]))	
   	y<-new("APPT.list",alignment=alignment(x),pheno=pheno(x),properties=properties(x),tree=tree(x),columns=columns(x)[i],multiMCMCglmm.list=new.multiMCMCglmm)
	y
})


#############################################################################
## setGeneric("cer.APPT.list")
#############################################################################

setGeneric("cer.APPT.list",function(x,class.var,which.columns=NULL,maxngaps=0,nummin=2,...) standardGeneric("cer.APPT.list"))

setMethod("cer.APPT.list", c("APPT.list","character","ANY","ANY","ANY"),
          function(x, class.var,which.columns=NULL,maxngaps=0,nummin=2,...)
{
	
	####	CONDITIONAL ENTROPY REDUCTION		####
	prop<-properties(x)
	if (is.null(which.columns)){which.columns<-columns(x)}

	ma.red<-alignment(x)[,which.columns]
	
	col2anal<-which(unlist(lapply(apply(ma.red,2,table),length))>=nummin & apply(ma.red,2,function(x)return(sum(x=="-")))<=maxngaps)

	####	CONDITIONAL ENTROPY	####
	ce<-function(x,y){
		t1<-table(y,x)/length(x)
		px<-apply(t1,2,sum)
		return( -sum(t(log(t(t1)/px,2))*t1,na.rm=TRUE))
	}
	
	vect.ce<-vector("numeric")
	hy<-vector("numeric")
	for(i in 1:length(col2anal)){
		vect.ce[i]<-ce(pheno(x)[,class.var],toupper(ma.red[,col2anal[i]]))
		hy[i]<- -sum((table(ma.red[,col2anal[i]])/dim(ma.red)[1])*(log(table(ma.red[,col2anal[i]])/dim(ma.red)[1],2)))
	}
	
	reduct.H<- 100*(1-vect.ce/hy)
	names(reduct.H)<-colnames(ma.red[,col2anal])
	reduct.H		
})




#############################################################################
## setGeneric("anovaAPPT.list")
#############################################################################

setGeneric("anovaAPPT.list",function(x,class.var,which.columns=NULL,maxngaps=0,nummin=2,...) standardGeneric("anovaAPPT.list"))

setMethod("anovaAPPT.list", c("APPT.list","formula","ANY","ANY","ANY"),
          function(x, class.var,which.columns=NULL,maxngaps=0,nummin=2,...)
{
	
	####	ANOVA				####
	prop<-properties(x)
	if (is.null(which.columns)){which.columns<-columns(x)}

	ma.red<-alignment(x)[,which.columns]
	
	col2anal<-which(unlist(lapply(apply(ma.red,2,table),length))>=nummin & apply(ma.red,2,function(x)return(sum(x=="-")))<=maxngaps)
	tab.anova<-matrix(0,length(col2anal),length(prop))
	
	for (m in 1:length(prop)) {
	
		listaDep<-lapply(col2anal,function(z) prop[[m]][aaa(toupper(ma.red[,z]))])
		
		for (i in 1:length(listaDep)) {
		      Y<-data.frame(DV=listaDep[[i]],pheno(x))
		      modelD<-lm(formula(paste("DV ~", as.character(class.var)[2])),data=Y)
		      tab.anova[i,m]<-summary(modelD)$coefficients[2,4]
		}
	}
	colnames(tab.anova)<-names(properties(x))
	rownames(tab.anova)<-names(listaDep)
	
	tab.anova<-as.data.frame(tab.anova)
	tab.anova$hm.signif<-apply(tab.anova,1,function(x)return(sum(x<1e-2)))
	tab.anova$median<-apply(tab.anova[,1:length(prop)],1,function(x)return(median(-log(x,10))))
	tab.anova$mean<-apply(tab.anova[,1:length(prop)],1,function(x)return(mean(-log(x,10))))

	tab.anova
			
})



#############################################################################
## setGeneric("MCMCglmm.APPT.list")
#############################################################################

setGeneric("MCMCglmm.APPT.list",function(x, class.var,random.eff,nitt=5.5e4,burnin=5e3,prior,scale=FALSE,parallel=FALSE,maxngaps=0,nummin=2,count=4,pr=FALSE, which.columns=NULL) 
standardGeneric("MCMCglmm.APPT.list"))

setMethod("MCMCglmm.APPT.list", c("APPT.list","formula","character","numeric","numeric","list","logical","logical","ANY","ANY","ANY","ANY","ANY"),
          function(x, class.var,random.eff,nitt=5.5e4,burnin=5e4,prior,scale=FALSE,parallel=FALSE,maxngaps=0,nummin=2,count=4,pr=FALSE,which.columns=NULL)
{

	prop<-properties(x)
	if (is.null(which.columns)){which.columns<-columns(x)}

	ma.red<-alignment(x)[,which.columns]
	
	col2anal<-which(unlist(lapply(apply(ma.red,2,table),length))>=nummin & apply(ma.red,2,function(x)return(sum(x=="-")))<=maxngaps)

	
	pheno<-pheno(x)
	tree<-tree(x)
	
	toto<-list()
	for (m in 1:length(prop)) {
	
		listaDep<-lapply(col2anal,function(z) prop[[m]][aaa(toupper(ma.red[,z]))])
		
		titi<-list()
		testObject <- function(object){exists(as.character(substitute(object)))}
		
		if (parallel==TRUE) {
			
			cl <- startMPIcluster(count=count,comm=7)
			registerDoMPI(cl)
			n<-1			
			titi<-foreach(n=1:length(listaDep),.packages=c('MCMCglmm'))%dopar%{
				Y<-data.frame(DV=listaDep[[n]],pheno)
			    Y$animal<-Y[,random.eff]
			    MCMCglmm(formula(paste("DV ~", as.character(class.var)[2])),data=Y,random=~animal,family="gaussian",rcov=~units,pedigree=tree,nitt=nitt,prior=prior[[m]],scale=scale,burnin=burnin,verbose=FALSE,pr=pr)
			     
			}
			closeCluster(cl)
		} else {
			
			for (i in 1:length(listaDep)) {
			      Y<-data.frame(DV=listaDep[[i]],pheno(x))
			      Y$animal<-Y[,random.eff]
			      modelD<-MCMCglmm(formula(paste("DV ~", as.character(class.var)[2])),data=Y,random=~animal,family="gaussian",rcov=~units,pedigree=tree(x),nitt=nitt,prior=prior[[m]],scale=scale,burnin=burnin,verbose=FALSE,pr=pr)
			      titi[[i]]<-modelD
			}
		}
				
		
		
		names(titi)<-names(listaDep)
		x@columns<-as.numeric(names(listaDep))
		x@multiMCMCglmm.list[[m]]<-titi
	}
	names(x@multiMCMCglmm.list)<-names(properties(x))
	x			
})


#############################################################################
## setGeneric("summaryAPPT.list")
#############################################################################

setGeneric("summaryAPPT.list",function(x,contrast,class.var,what.prop) standardGeneric("summaryAPPT.list"))

setMethod("summaryAPPT.list", c("APPT.list","character","character","ANY"),
          function(x,contrast,class.var,what.prop=NULL)
{
	
	if(is.null(what.prop)){what.prop<-1:length(properties(x))}
	
	tabla<-matrix(0,length(multiMCMCglmm(x)[[1]]),length(what.prop))
	tabSize<-matrix(0,length(multiMCMCglmm(x)[[1]]),length(what.prop))
	lrange<-lapply(properties(x)[what.prop],function(z) range(z)[2]-range(z)[1])
	toto<-multiMCMCglmm(x)[what.prop]
	prop<-names(toto)
	
	# gt0: tabla
	# Effect Size:	tabSize
	for (m in 1:length(prop)) {
		for (i in 1:length(multiMCMCglmm(x)[[1]])) {
			tabla[i,m]<-sum((toto[[m]][[i]]$Sol[,contrast[1]]-toto[[m]][[i]]$Sol[,contrast[2]])>0)/length(toto[[m]][[i]]$Sol[,contrast[2]])	
			
			tabSize[i,m]<-median(toto[[m]][[i]]$Sol[,contrast[1]]-toto[[m]][[i]]$Sol[,contrast[2]])/lrange[[m]]	
		}
	}
	
	colquevan<-columns(x)
	
	colnames(tabla)<-prop
	rownames(tabla)<-colquevan
	
	# CHI-SQUARE
	tabcor<-2*(tabla-0.5)
	tabCS<-(tabcor*sqrt(dim(toto[[1]][[1]]$Sol)[1]))^2

	tabEff<-matrix(0,length(multiMCMCglmm(x)[[1]]),length(what.prop))
	for (i in 1:length(what.prop)){tabEff[,i]<-unlist(lapply(multiMCMCglmm(x)[[what.prop[i]]],function(z)return(median(summary(z)$solutions[,4]))))}

	tabCSeff<-(tabcor^2)*tabEff

	resuCS<-apply(tabCS,1,sum)
	resuCSeff<-apply(tabCSeff,1,sum)
	
	colnames(tabSize)<-prop
	rownames(tabSize)<-colquevan

	
	sumtr<-apply(tabcor^2,1,sum)/dim(tabcor)[2]
	names(sumtr)<-rownames(tabcor)
	

	liAAs<-list()
	for (i in 1:length(colquevan)) { 	
		mdf<-data.frame(seq=row.names(alignment(x)),pos=alignment(x)[,colquevan[i]],patho=pheno(x)[,class.var])
		liAAs[[i]]<-table(mdf$patho,mdf$pos)/apply(table(mdf$patho,mdf$pos),1,sum)	
	}
	names(liAAs)<-names(sumtr)
	rm(toto)

	out<-list()
	out[[1]]<-tabcor
	out[[2]]<-tabSize
	out[[3]]<-sumtr
	out[[4]]<-resuCS
	out[[5]]<-resuCSeff
	out[[6]]<-liAAs
	names(out)<-c("tabcor","tabSize","SumTr","ChiSq","ChiSq.eff","AAlist")
	out
	
})



#############################################################################
## setGeneric("RGBcolor.APPT.list")
#############################################################################

setGeneric("RGBcolor.APPT.list",function(x,which.prop=0,write.score=FALSE,...) standardGeneric("RGBcolor.APPT.list"))

setMethod("RGBcolor.APPT.list", c("APPT.list","numeric","logical"),
          function(x,which.prop=0,write.score=FALSE,...)
{
	
	std01<-function(ms1){return((ms1-(min(ms1)))/(max(ms1)-min(ms1)))}
	####	COLORS		####
	if (which.prop==0) {
		listprop<-properties(x)
	} else {
		listprop<-properties(x)[which.prop]	
	}
	
	matstd01<-matrix(0,20,3)
	rownames(matstd01)<-names(listprop[[1]])
	colnames(matstd01)<-c("R","G","B")	
	
	
	if (length(listprop)>3) {
	
		matprop<-matrix(0,20,length(listprop))
		for (i in 1:length(listprop)) {matprop[,i]<-listprop[[i]]}
		rownames(matprop)<-names(listprop[[1]])
		colnames(matprop)<-names(listprop)
		mipca.prop<-princomp(matprop,cor=TRUE)
		if (write.score==TRUE) {write.table(mipca.prop$score,file="colors_pca_score.txt")}
		
		for (i in 1:3){matstd01[,i]<-std01(mipca.prop$score[,i])}

	
	} else {
	
		for (i in 1:length(listprop)) {
			matstd01[,i]<-std01(listprop[[i]])
		}
	}
	colAA<-rgb(matstd01[,1],matstd01[,2],matstd01[,3])
	names(colAA)<-names(listprop[[1]])
	
	ma.red<-alignment(x)[,columns(x)]
	ma.col<-matrix("",dim(ma.red)[1],dim(ma.red)[2])
	for (j in 1:dim(ma.red)[2]) {
		
		ma.col[,j]<-colAA[aaa(ma.red[,j])]	
	}
	rownames(ma.col)<-rownames(ma.red)
	colnames(ma.col)<-columns(x)
	ma.col
	#barplot(1:20,col=colAA,names.arg=names(colAA),las=2)
	
})


#############################################################################
## setGeneric("bootMCMCglmm.APPT.list")
#############################################################################

setGeneric("bootMCMCglmm.APPT.list",function(x, class.var,boot=100,contrast,random.eff,nitt=1e5,burnin=5e4,prior,scale=FALSE,parallel=FALSE,what.prop=NULL,which.column=NULL,maxngaps=0,nummin=2,count=4,...) standardGeneric("bootMCMCglmm.APPT.list"))

setMethod("bootMCMCglmm.APPT.list", c("APPT.list","formula","numeric","character","character","numeric","numeric","list","logical","logical","ANY","ANY","ANY","ANY","ANY"),
          function(x, class.var,boot=100,contrast,random.eff,nitt=1e5,burnin=5e4,prior,scale=FALSE,parallel=FALSE,what.prop=NULL,which.column=NULL,maxngaps=0,nummin=2,count=4,...)
{

	if(is.null(what.prop)){what.prop<-1:length(properties(x))}
	if (is.null(what.prop)) {prop<-properties(x)} else {prop<-properties(x)[what.prop]}
	ma.red<-alignment(x)[,which.column]
	
	col2anal<-which.column*(length(table(ma.red))>=nummin & sum(ma.red=="-")<=maxngaps)	
	
	pheno<-pheno(x)
	tree<-tree(x)

	lisa<-list()
	for (r in 1:boot){lisa[[r]]<-sample(1:length(ma.red))}

	
	toto<-list()
	for (m in 1:length(prop)) {
	
		listaDep<-prop[[m]][aaa(toupper(ma.red))]

		titi<-list()
		testObject <- function(object){exists(as.character(substitute(object)))}
		
				
		if (parallel==TRUE) {
			
			cl <- startMPIcluster(count=count,comm=7)
			registerDoMPI(cl)
			n<-1			
			titi<-foreach(n=1:boot,.packages=c('MCMCglmm'))%dopar%{
				Y<-data.frame(DV=listaDep[lisa[[n]]],pheno)
			    Y$animal<-Y[,random.eff]
			    MCMCglmm(formula(paste("DV ~", as.character(class.var)[2])),data=Y,random=~animal,family="gaussian",rcov=~units,pedigree=tree,nitt=nitt,prior=prior[[m]],scale=scale,burnin=burnin,verbose=FALSE)
			    
			}
			closeCluster(cl)
		} else {
			
			for (i in 1:boot) {
			      Y<-data.frame(DV=sample(listaDep),pheno(x))
			      Y$animal<-Y[,random.eff]
			      modelD<-MCMCglmm(formula(paste("DV ~", as.character(class.var)[2])),data=Y,random=~animal,family="gaussian",rcov=~units,pedigree=tree(x),nitt=nitt,prior=prior[[m]],scale=scale,burnin=burnin,verbose=FALSE)
			      titi[[i]]<-modelD
			}
		}
						
		names(titi)<-1:boot
		toto[[m]]<-titi
	}
	names(toto)<-names(prop)
	
	#SUMMARY
	tabla<-matrix(0,boot,length(what.prop))
	tabSize<-matrix(0,boot,length(what.prop))
	lrange<-lapply(properties(x)[what.prop],function(z) range(z)[2]-range(z)[1])
	
	# gt0: tabla
	# Effect Size:	tabSize
	for (m in 1:length(prop)) {
		for (i in 1:boot) {
			tabla[i,m]<-sum((toto[[m]][[i]]$Sol[,contrast[1]]-toto[[m]][[i]]$Sol[,contrast[2]])>0)/length(toto[[m]][[i]]$Sol[,contrast[2]])	
			
			tabSize[i,m]<-median(toto[[m]][[i]]$Sol[,contrast[1]]-toto[[m]][[i]]$Sol[,contrast[2]])/lrange[[m]]	
		}
	}
	
	colquevan<-which.column
	
	colnames(tabla)<-names(prop)
	rownames(tabla)<-1:boot
	
	# CHI-SQUARE
	tabcor<-2*(tabla-0.5)
	tabCS<-(tabcor*sqrt(dim(toto[[1]][[1]]$Sol)[1]))^2
	resuCS<-apply(tabCS,1,sum)
	
	colnames(tabSize)<-names(prop)
	rownames(tabSize)<-1:boot

	sumtr<-apply(tabcor^2,1,sum)/dim(tabcor)[2]
	names(sumtr)<-1:boot
	
	rm(toto)

	out<-list()
	out[[1]]<-tabcor
	out[[2]]<-tabSize
	out[[3]]<-sumtr
	out[[4]]<-resuCS
	names(out)<-c("tabcor","tabSize","SumTr","ChiSq")
	out	
				
})



