#takes a given consensus topology and assigns branch lengths based on the average of those source trees that have that branch
library(phylobase)

isClade<-function(phy,taxonlist) {
	if (length(taxonlist)==1) {
		return(TRUE)
	}
	if(length(descendants(phy,MRCA(phy,taxonlist),type="tips"))==length(taxonlist)) {
		return(TRUE)
	}
	else {
		return(FALSE)
	}
}

edgeLengthTaxset<-function(phy,taxonlist) {
	return(edgeLength(phy,MRCA(phy,taxonlist)))
}

summarizeNode<-function(nodeId,focalTree,sourceTreeList,print.progress) {
	matchingVector<-unlist(lapply(sourceTreeList,isClade,descendants(focalTree,nodeId,type="tips")))
	proportion<-sum(matchingVector) / length(matchingVector)
	lengths<-unlist(lapply(sourceTreeList[matchingVector],edgeLengthTaxset,descendants(focalTree,nodeId,type="tips")))
	result<-c(nodeId,NA,NA,NA,NA)
	if (sum(is.na(lengths))<length(lengths)) {
		result<-c(nodeId,proportion,mean(lengths,na.rm=TRUE),median(lengths,na.rm=TRUE),sd(lengths,na.rm=TRUE))
	}
	if (print.progress) {
		print(result)
	}
	names(result)<-c("nodeId","proportion","mean_brlen","median_brlen","sd_brlen")
	return(result)
}

#trees should be rooted
consensusBrlen<-function(focalTree,sourceTreeList,type=c("proportion","mean_brlen","median_brlen","sd_brlen"),print.progress=TRUE,return.val="tree") {
	type<-match.arg(type)
	if (class(focalTree)!="phylo4") {
		focalTree<-as(focalTree,"phylo4")
	}
	if (class(sourceTreeList[[1]])!="phylo4") {
		sourceTreeList<-lapply(sourceTreeList,as,"phylo4")
	}
	allNodes<-nodeId(focalTree,"all")
	allNodes<-allNodes[which(allNodes!=nodeId(focalTree,"root"))] #do not care about root edge
	if (print.progress) {
		print(c("nodeId","proportion","mean_brlen","median_brlen","sd_brlen"))
	}
	allResults<-sapply(allNodes,summarizeNode,focalTree,sourceTreeList,print.progress)
	if (return.val=="tree") {
		newEdgeLengths<-edgeLength(focalTree)
		newNodeLabels<-nodeLabels(focalTree)
		for (nodeIndex in 1:length(allNodes)) {
			newLength<-allResults[which(row.names(allResults)==type),nodeIndex]
			if (is.na(newLength)) {
				newLength=0
			}
			newEdgeLengths[ which(names(newEdgeLengths)==getEdge(focalTree,allNodes[nodeIndex])) ]<-newLength
			newNodeLabels[ which(names(newNodeLabels)==allNodes[nodeIndex]) ] <- round(allResults[which(row.names(allResults)=="proportion"),nodeIndex],2)
		}
		edgeLength(focalTree)<-newEdgeLengths
		nodeLabels(focalTree)<-newNodeLabels
		return(focalTree)
	}
	else {
		return(allResults)
	}
}
