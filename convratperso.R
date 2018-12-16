##########################################################################################
##########################################################################################
###		rewriting of Stayton's function convrat (convevol R package; 2014), 	##
###		implement some Liam Revell code http://blog.phytools.org/ 		##
###		Analysis from Botton-Divet et al. 2016 					##
###		please contact lbottondivet@mnhn.fr if needed				##
###		or visit https://github.com/leobd 					##
##########################################################################################
##########################################################################################
	
convratperso<-function(phyl,phendata,convtips)
{
## Getting numerical indexes
conv1<-convtips[1]
conv2<-convtips[2]
ifelse(is.character(conv1)==T,tipvalue1<-which(phyl[["tip.label"]]==conv1),tipvalue1<-conv1)
ifelse(is.character(conv2)==T,tipvalue2<-which(phyl[["tip.label"]]==conv2),tipvalue2<-conv2)



## Compute lineages
require(ape)
mrcatips<-getMRCA(phyl, convtips)

	i=tipvalue1
	lineage1<-tipvalue1
	while(i!=mrcatips)
		{
			ancestor<-phyl$edge[which(phyl$edge[,2]==i),1]
			lineage1<-c(lineage1,ancestor)
			i<-ancestor
		}

	j=tipvalue2
	lineage2<-tipvalue2
	while(j!=mrcatips)
		{
			ancestor<-phyl$edge[which(phyl$edge[,2]==j),1]
			lineage2<-c(lineage2,ancestor)
			j<-ancestor
		}
		
## Computation of phenotypic values for internal nodes
	phentot<-multiancperso(phyl,phendata)# huge computation, could be improved by computing the values of the nodes under the MRCA of considered tips only. Already slightly improved by parallelisation in multiancperso
	rownames(phentot)<-c(rownames(phentot)[1:length(phyl$tip.label)],c((1+length(phyl$tip.label)):length(rownames(phentot))))

## Computation of the phenotypic distances
	dt<-dist(phentot[unique(c(lineage1,lineage2)),])

## Extraction of maximal values
	maxval<-dt[which.max(dt)]
	tipdist<-dist(phendata[convtips,])

	C1<-1-(tipdist/maxval)
	C2<-maxval-tipdist
## Computation of phenotypic distances branches along lineages
distlineage1<-0
for (i in 1:(length(lineage1)-1))
	{
	distlineage1<-distlineage1+dist(phentot[c(lineage1[i],lineage1[i+1]),])
	}
distlineage2<-0
for (i in 1:(length(lineage2)-1))
	{
	distlineage2<-distlineage2+dist(phentot[c(lineage2[i],lineage2[i+1]),])
	}


C3<-C2/(distlineage1+distlineage2)

## Computation of the total phenotypic distances in the lineages defined by MRCA of convergent tips

	getDescendants<-function(tree,node,curr=NULL){				##Fonction from Liam Revell's blog http://blog.phytools.org/
	if(is.null(curr)) curr<-vector()
	daughters<-tree$edge[which(tree$edge[,1]==node),2]
	curr<-c(curr,daughters)
	w<-which(daughters>=length(tree$tip))
	if(length(w)>0) for(i in 1:length(w))
	curr<-getDescendants(tree,daughters[w[i]],curr)
	return(curr)
	}

subtr<-getDescendants(phyl,mrcatips)

a<-list()
for (i in 1:length(subtr))	
	{
	a[[i]]<-(which(phyl[["edge"]][,1]==subtr[i]|phyl[["edge"]][,2]==subtr[i]))
	}
phy2dist<-unique(unlist(a))

## Computation of the total phenotypic distances in the subtree defined by MRCA, including side branches that are not part of tips lineage
distsub<-NULL
for(i in 1:length(phy2dist))
	{
	distsub[i]<-dist(rbind(phentot[phyl$edge[phy2dist,][i,1],],phentot[phyl$edge[phy2dist,][i,2],]))
	}
totalchanges<-sum(distsub)

C4<-C2/totalchanges

output<-c(C1,C2,C3,C4)
names(output)<-c("C1","C2","C3","C4")
output
}
