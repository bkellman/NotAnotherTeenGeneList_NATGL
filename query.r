source('query.functions.r')
library(biomaRt)
library(xlsx)
library(KEGGREST)
library(gplots)
library(ggplot2)

override=FALSE
signed=FALSE

#ensembl=useMart("ensembl",'hsapiens_gene_ensembl')
ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

r1 <- read.xlsx('geneList/Courchesne_summary_10062015.xlsx',sheetIndex=1)
r2 <- read.csv('geneList/tito_list_clinical_data.txt')[,1:5]
gene_ids <- unique(as.character(r1$accession))
gene_ids <- unlist(lapply(strsplit(gene_ids,'\\.'),function(x) x[1]))
r <- merge(r1,r2,by.x='Individual',by.y='subjectid')

gene_type = 'refseq_mrna'
attributes <- c(gene_type)

# get gene_ids
attributes <- c(attributes,'hgnc_symbol','ensembl_gene_id','entrezgene')#,'hsa ')

# get phenotype (OMIN)
#attributes <- c(attributes,'mim_morbid_description') # ,'mim_gene_accession','mim_gene_description')

# get ontology (GO)
attributes <- c(attributes,'go_id')

# run BM query
BM = getBM(attributes=attributes, filters=gene_type, values=gene_ids, mart=ensembl)
gene_ids <- unique(BM$refseq_mrna)

# get KEGG pathwasy
#p = keggLink('hsa','pathway')
#map = split(names(p),unname(p))

# cluster by tissue expression (BioGPS, GNF, GeTX)
GTEx_file <- '../data/DATA_Autism_Genomic_Varients/GTEx/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
if(!file.exists('../data/DATA_Autism_Genomic_Varients/tmp.txt') | override){
	system(paste("grep -E '",paste(c('Name',unique(BM$ensembl_gene_id)),collapse="|"),GTEx_file,"' > ../data/DATA_Autism_Genomic_Varients/tmp.txt",sep=''))
}
expression_annot <- read.csv('../data/DATA_Autism_Genomic_Varients/GTEx/GTEx_Data_V6_Annotations_SampleAttributesDS.txt',sep='\t')
expression <- read.csv('../data/DATA_Autism_Genomic_Varients/tmp.txt',sep='\t')
rownames(expression) <- unlist(lapply(strsplit(as.character(expression$Name),'\\.'),function(x) x[1]))
expression <- expression[,-c(1:2)]
expression <- data.matrix(expression) # make data matrix
expression <- log(expression+1) # approximate gaussian
expression <- scale(expression) # normalize values
#tissue_site <- data.frame(do.call(rbind,strsplit(colnames(expression),'\\.')))$X3
tissue_site_df <- expression_annot[match(gsub('\\.','-',colnames(expression)),as.character(expression_annot$SAMPID)),c('SMTS','SMTSD')]
tissue_site <- tissue_site_df$SMTSD

#heatmap.2(data.matrix(expression),
#	ColSideColors=rainbow(length(levels(tissue_site)))[tissue_site],trace='none',
#	labRow=convert(rownames(expression),'ensembl_gene_id','hgnc_symbol',BM ))
tissue_site <- tissue_site[order(tissue_site)]
tissue_site_nam <- ifelse( duplicated(tissue_site) , '' , as.character(tissue_site))
pdf('clusters_out/tissue_expression.pdf',width=15)
heatmap.2(data.matrix(expression)[,order(tissue_site)],Colv=FALSE,
	ColSideColors=c(rainbow(length(levels(tissue_site)))[tissue_site]),trace='none',
	labRow=convert(rownames(expression),'ensembl_gene_id','hgnc_symbol',BM ),
	labCol=tissue_site_nam,mar=c(10,5))
dev.off()


######### compute distance
dL <- list()

# GO distance
dL[['GO_jaccard']] <- 1 - GO_jaccard(BM,gene_ids,gene_type)
dL[['GO_jaccard_ramigo_flat']] <- 1 - GO_jaccard(BM,gene_ids,gene_type,ramigo_search=TRUE,ramigo_paths=FALSE)
dL[['GO_jaccard_ramigo_paths']] <- 1 - GO_jaccard(BM,gene_ids,gene_type,ramigo_search=FALSE,ramigo_paths=TRUE)


# Tissue expression
dL[['abs_tissue_cor_pearson']] <- 1-abs( cor(t(expression[,-c(1:2)]+1)) )


# GeneMania combined 
# wget -r --no-parent http://genemania.org/data/current/Homo_sapiens/
GM_comb <- load_GeneMania(combined=TRUE)
normalize <- function(x) x/max(x)
dL[['geneMania_combined_net_weightedShortestPath']] <- Network_ShortestPaths_distance(unique(BM$ensembl_gene_id),GM_comb,weighted=TRUE,func=normalize)
dL[['geneMania_combined_net_unweightedShortestPath']] <- Network_ShortestPaths_distance(unique(BM$ensembl_gene_id),GM_comb,weighted=FALSE,func=normalize)


#types<-c('Physical_interactions','Predicted','Genetic_Interactions','Pathway','Co-localization','Co-expression','Shared_protein_domains')
types<-c("Pathway","Co-expression","Co-localization","Genetic_interactions","Physical_interactions","Predicted","Shared_protein_domains")
# GeneMania individual
for(t_i in types){
	GM_df <- load_GeneMania(combined=FALSE,type=t_i)
	dL[[paste('geneMania_',t_i,'_net_weightedShortestPath',sep='')]] <- Network_ShortestPaths_distance(unique(BM$ensembl_gene_id),GM_df,weighted=TRUE,func=normalize)
	dL[[paste('geneMania_',t_i,'_net_unweightedShortestPath',sep='')]] <- Network_ShortestPaths_distance(unique(BM$ensembl_gene_id),GM_df,weighted=FALSE,func=normalize)	
}


### save
save(dL,file='../data/DATA_Autism_Genomic_Varients/dL.rda')

### combine distance matrixes
dL[['final']] <- combine(dL,sum,BM,gene_ids,gene_type)
for(i in 1:length(dL)){
	d_i = dL[[i]]
	d_i_id <- names(which.max(apply(BM,2,function(x) sum(rownames(d_i) %in% x))))
	if(!all(rownames(d_i)==colnames(d_i))){stop('distance matrix is not symetrical')}
	rownames(d_i)=colnames(d_i)=convert(rownames(d_i),d_i_id,'hgnc_symbol',BM)
	dL[[i]] = d_i
}

#rownames(d)=colnames(d)=convert(rownames(d),gene_type,'hgnc_symbol',BM)

### cluster
par(mfrow=c(1,3))
plot(hclust(as.dist(dL[[1]])),main='jaccard GO cluster')
plot(hclust(as.dist(dL[[2]])),main='jaccard GO+Ramigo cluster')
plot(hclust(as.dist(dL[[3]])),main='shortest paths GO+Ramigo cluster')

plot(hclust(as.dist(dL[[4]])),main='tissue expression correlation cluster')

plot(hclust(as.dist(d)),main='combined cluster')

# PCoA
library(ape)
x<-pcoa(as.dist(d), correction="none", rn=NULL)
 
## S3 method for class 'pcoa':
biplot(x, Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL)



