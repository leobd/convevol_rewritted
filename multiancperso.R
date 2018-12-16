##########################################################################################
##########################################################################################
###			            Stayton's multianc (convevol R package; 2014) parallelized		      ##
###		                      	from Botton-Divet et al. 2016 				                   	##
###			                   please contact lbottondivet@mnhn.fr if needed		           	##
###			                      or visit https://github.com/leobd 				                ##
##########################################################################################
##########################################################################################
## Required for parallel computing
require(doMC)
registerDoMC(detectCores())
##
multiancperso<-function (phyl, phendata) 
{
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (nrow(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    if (is.null(rownames(phendata))) {
        warning("no row names for data.  Assuming that the rows are in the same order as tips.")
        rownames(X) <- phyl$tip.label
    }
    firstvar <- fastAnc(phyl, phendata[, 1])
    
  allancstates<-matrix(data=unlist(foreach(i=1:ncol(phendata))%dopar%{fastAnc(phyl,phendata[,i])}),nrow=length(firstvar),ncol=ncol(phendata),byrow=F)
    colnames(allancstates) <- colnames(phendata)
    alldata <- rbind(phendata, allancstates)
}
