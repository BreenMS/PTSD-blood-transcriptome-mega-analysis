library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads(n=34)

###################################################################################################################################################
#INPUT Gene expression matrix
mat <- read.delim("Matrix_Consensus.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(mat) =names(mat)
dataExpr0<-as.data.frame(t(mat))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

#INPUT meta-data file
traitData <- read.delim("Matrix_Consensus_targets.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)
head(traitData)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$ID.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) 
names(datTraits)

powers=c(1:30) # in practice this should include powers up to 30
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

pdf("WGCNA_Consensus_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.80,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyPre = adjacency((dataExpr0),power=13 ,type="signed") #A Power of 13 was used to build a consensus network
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

#MODULE ASSIGNMENTS
mColorh=NULL
for (ds in 0:4){
 tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
   minClusterSize = (25), cutHeight = 0.9999, 
   deepSplit = ds, distM = dissTOMPre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

pdf("WGCNA_Consensus_DeepSplit.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,4]
table(modulesPRE)

PTSD = as.data.frame(datTraits$PTSD)
names(PTSD)="PTSD"
GS.PTSD=as.numeric(cor(dataExpr0,PTSD,use="p"))
#GS.PTSD_abs=abs(GS.PTSD) #absolute Value
GS.PTSDColor=numbers2colors(GS.PTSD,signed=T)

Control = as.data.frame(datTraits$Control)
names(Control)="Control"
GS.Control=as.numeric(cor(dataExpr0,Control,use="p"))
#GS.Control_abs=abs(GS.Control) #absolute Value
GS.ControlColor=numbers2colors(GS.Control,signed=T)

Combat_PTSD = as.data.frame(datTraits$Combat_PTSD)
names(Combat_PTSD)="Combat_PTSD"
GS.Combat_PTSD=as.numeric(cor(dataExpr0,Combat_PTSD,use="p"))
#GS.Combat_PTSD_abs=abs(GS.Combat_PTSD) #absolute Value
GS.Combat_PTSDColor=numbers2colors(GS.Combat_PTSD,signed=T)

Combat_Cntl = as.data.frame(datTraits$Combat_Cntl)
names(Combat_Cntl)="Combat_Cntl"
GS.Combat_Cntl=as.numeric(cor(dataExpr0,Combat_Cntl,use="p"))
#GS.Combat_Cntl_abs=abs(GS.Combat_Cntl) #absolute Value
GS.Combat_CntlColor=numbers2colors(GS.Combat_Cntl,signed=T)

IP_PTSD = as.data.frame(datTraits$IP_PTSD)
names(IP_PTSD)="IP_PTSD"
GS.IP_PTSD=as.numeric(cor(dataExpr0,IP_PTSD,use="p"))
#GS.IP_PTSD_abs=abs(GS.IP_PTSD) #absolute Value
GS.IP_PTSDColor=numbers2colors(GS.IP_PTSD,signed=T)

IP_Cntl = as.data.frame(datTraits$IP_Cntl)
names(IP_Cntl)="IP_Cntl"
GS.IP_Cntl=as.numeric(cor(dataExpr0,IP_Cntl,use="p"))
#GS.IP_Cntl_abs=abs(GS.IP_Cntl) #absolute Value
GS.IP_CntlColor=numbers2colors(GS.IP_Cntl,signed=T)

Gender = as.data.frame(datTraits$Gender)
names(Gender)="Gender"
GS.Gender=as.numeric(cor(dataExpr0,Gender,use="p"))
#GS.Gender_abs=abs(GS.Gender) #absolute Value
GS.GenderColor=numbers2colors(GS.Gender,signed=T)

CntlClusters = as.data.frame(datTraits$CntlClusters)
names(CntlClusters)="CntlClusters"
GS.CntlClusters=as.numeric(cor(dataExpr0,CntlClusters,use="p"))
#GS.CntlClusters_abs=abs(GS.CntlClusters) #absolute Value
GS.CntlClustersColor=numbers2colors(GS.CntlClusters,signed=T)

PTSDCluster1 = as.data.frame(datTraits$PTSDCluster1)
names(PTSDCluster1)="PTSDCluster1"
GS.PTSDCluster1=as.numeric(cor(dataExpr0,PTSDCluster1,use="p"))
#GS.PTSDCluster1_abs=abs(GS.PTSDCluster1) #absolute Value
GS.PTSDCluster1Color=numbers2colors(GS.PTSDCluster1,signed=T)

PTSDCluster2 = as.data.frame(datTraits$PTSDCluster2)
names(PTSDCluster2)="PTSDCluster2"
GS.PTSDCluster2=as.numeric(cor(dataExpr0,PTSDCluster2,use="p"))
#GS.PTSDCluster2_abs=abs(GS.PTSDCluster2) #absolute Value
GS.PTSDCluster2Color=numbers2colors(GS.PTSDCluster2,signed=T)


datColors0=data.frame(modulesPRE, GS.PTSDColor,GS.ControlColor,GS.Combat_PTSDColor, GS.Combat_CntlColor, GS.IP_PTSDColor, GS.IP_CntlColor, GS.GenderColor)

pdf("WGCNA_Consensus_Colored_Modules.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Case Gene Dendrogram and Module Colors", groupLabels=c("Module colors","PTSD", "Control","Combat PTSD", "Combat Control", "IP PTSD", "IP Control", "Gender"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()


#Check to see if network modules can be cut and merged...
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")

pdf("Module_RelationshipsIP.pdf")
plot(METree, main ="Clustering of Module Eigengenes",xlab ="",sub="")
abline(h=0.15,col="red")
abline(h=0.25,col="red")
abline(h=0.2,col="red")
abline(h=0.3,col="red")
abline(h=0.35,col="red")
dev.off()
####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
MEDissThres=0.25
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PreDeploy_Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreeControlPre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")

##########################################################################################################################

#CALCULATE PC FOR VISUALIZATION FOR CASE PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names = row.names((dataExpr0))

pdf("WGCNA_Consensus_Module_Characters.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesPRE)))
{
 par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
 plotMat(t(scale(dataExpr0[,modulesPRE==which.module])),
,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))

  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.5, font=2)
};
dev.off();



##Module enrichment analysis based on userListEnrichment function
Gene  = colnames(dataExpr0)

enrichments = userListEnrichment(Gene, modulesPRE,
fnIn = NULL,
catNmIn = NULL,
#fnIn = c("GeneList","ModuleColors"),
#catNmIn =  c("Genes","Modules"),
nameOut = "WGCNA_Consensus_Module_Enrichment.csv", 
useBrainLists = TRUE,
useBloodAtlases = TRUE,
omitCategories = "grey", 
outputCorrectedPvalues = TRUE,
useStemCellLists = FALSE, 
outputGenes = TRUE, 
minGenesInCategory = 5, 
useBrainRegionMarkers = FALSE,
useImmunePathwayLists = TRUE,
usePalazzoloWang = TRUE)

# More detailed overlap information is in the pValue output. 
head(enrichments$pValue)

# To see a list of all significant enrichments
enrichments$sigOverlaps
 


#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs,PTSD))
pdf("WGCNA_Consensus_Correlations_Heatmap.pdf", h=16, w=15)
plotEigengeneNetworks(MET,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
write.table(MEs, "WGCNA_Consensus_MEs.txt", sep="\t") #Write out ME values for subsequent analyses. 

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
data = cbind(moduleTraitCor,moduleTraitPvalue)
pdf("WGCNA_Consensus_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("ME-Trait Relationships"))
dev.off()


#MODULE MEMBERSHIP (kME) KME is defined as correlation between expression and modules
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
geneModuleMembership1 = signedKME((dataExpr0), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"WGCNA_Consensus_Network.csv",row.names=FALSE)

topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 maxKMErank  = rank(apply(cbind(kMErank1+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
topGenesKME

#WRITE OUTPUT
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="WGCNA_Consensus_Network.csv") #Network output. 

