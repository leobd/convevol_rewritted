############################################################################################
## Rewritting of the convratperso function to work with more than two couples of species  ##
##      Adapted from Tristan Stayton's function convrat (convevol R package; 2014), 	    ##
############################################################################################
	
convratpar<-function(phyl,phendata,convtips)
{
## Required for parallel computing
require(doMC)
registerDoMC(detectCores())
## Computation of phenotypic values for internal nodes
  phentot<-multiancperso(phyl,phendata)#gros calcul pas necessaire d'avoir tout les ancestral state, pourrais être realisé sur le sous arbre defini par MRCAtips, déjà un peu amelioré en parallelisant multianc
  rownames(phentot)<-c(rownames(phentot)[1:length(phyl$tip.label)],c((1+length(phyl$tip.label)):length(rownames(phentot))))

## Single run of uniconvrat if convtips contains only one pair of convergent species
  if(length(convtips)==2)
    {
      out<-uniconvrat(phyl,phentot,convtips)
    }
    
## Parallel run of uniconvrat if convtips contains more than two species
  if(length(convtips)>2)
  {
   require(utils)
   paires<-combn(convtips,2)
   results<-matrix(data =unlist(foreach(i=1:dim(paires)[[2]]) %dopar% {uniconvrat(phyl, phentot, paires[,i])}),nrow = dim(paires)[[2]], ncol = 4, byrow = TRUE)
   rownames(results)<-paste(paires[1,],paires[2,],sep='-')
   moy<-colMeans(results)
   names(moy)<-c("C1","C2","C3","C4")
   out<-list(mean=moy,results=results)
  }
out
}
