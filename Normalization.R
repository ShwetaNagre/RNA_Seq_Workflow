## Filtering the counts file and renaming the samples

data=read.delim("/home/shweta//Downloads/Covid_test/covid_alignment/ReadCountsFinal.txt")

dataFilt <- subset(data, select = -c(Chr, Start, End, Strand,Length))

colnames(dataFilt)[which(names(dataFilt) == "SRR18829396.bam")] <- "COVID19_8_HC"
colnames(dataFilt)[which(names(dataFilt) == "SRR18829405.bam")] <- "COVID19_021_Severe"
colnames(dataFilt)[which(names(dataFilt) == "SRR18829438.bam")] <- "COVID19_034_Severe"
colnames(dataFilt)[which(names(dataFilt) == "SRR18829439.bam")] <- "COVID19_033_Severe"
colnames(dataFilt)[which(names(dataFilt) == "SRR18829462.bam")] <- "COVID19_2_HC"
colnames(dataFilt)[which(names(dataFilt) == "SRR18829463.bam")] <- "COVID19_1_HC"


#####################################################

## Extracting gene biotype and hgnc symbols form ensembl dataset

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
GeneInfo <- getBM(filters= "ensembl_gene_id",
                  attributes= c("ensembl_gene_id","gene_biotype","hgnc_symbol"),
                  values=dataFilt$Geneid,mart= ensembl)
dataFilt_Modified=merge(dataFilt,GeneInfo,by.x = "Geneid", by.y = "ensembl_gene_id")
write.table(dataFilt_Modified,file="/home/shweta//Downloads/Covid_test/covid_alignment/FinalData.txt",sep="\t",col.names = NA,quote = FALSE)


########################################################

## Extract protein coding genes. 

ProtCode <- subset(dataFilt_Modified, gene_biotype =="protein_coding")
rownames(ProtCode)=paste(ProtCode$Geneid,ProtCode$hgnc_symbol,sep="_")
ProtCode <- subset(ProtCode, select = -c(gene_biotype, hgnc_symbol, Geneid))

write.table(ProtCode,file="/home/shweta/Downloads/Covid_test/covid_alignment/ProteinCodingCount.txt",
            sep="\t",col.names = NA,quote = FALSE)


########################################################

## Normalization 

TotalCounts <- colSums(ProtCode)
CPM <- t(t(ProtCode) / TotalCounts * 1e6)

CPM <- CPM[rowSums(CPM != 0) > 0, ]

########################################################

## PCA

library(PCAtools)
library(ggplot2)
library(ggrepel)

m=read.delim("/home/shweta/Downloads/Covid_test/Metadata.csv",row.names = 1)

p <- pca(CPM, metadata = m, removeVar = 0.1,scale = TRUE)
p$rotated$Group=m$Group
write.table(p$rotated,file = "/home/shweta/Downloads/Covid_test/PCA_covid_test.txt",sep = "\t",col.names = NA,quote = FALSE)
head(p$variance)



data_pca=read.csv("/home/shweta/Downloads/Covid_test/PCA_covid_test.txt",row.names = 1,check.names = FALSE)

pdf("/home/shweta/Downloads/Covid_test/PCA_covid_test.pdf")

ggplot(data_pca, aes(PC1,PC2,color=Group, fill=Group))+geom_point(size=2,shape=21)+
  scale_color_manual(values=c(HC="#c9a711",Covid="#005485"))+
  scale_fill_manual(values=c(HC="#fcd116",Covid="#006aa7"))+
  labs(x="PC1",y="PC2")+
  theme(axis.title = element_text(size=10),legend.position = "right",plot.margin = margin(4.5,3,4.5,3, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(ncol =1,override.aes=list(fill="grey",color="grey")),color = guide_legend(ncol =1),fill = guide_legend(ncol =1))

dev.off()
