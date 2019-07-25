##################################################################################################################################
##### Transcriptomic variations along introduction pathway" ###
##################################################################################################################################

##### Working directory, libraries and functions #################################################################################

rm(list = ls())
DataDir = "/Users/antoinefelden/Documents/Research/Manuscripts/X-TransPathway/Analysis/01_data"
FigDir = "/Users/antoinefelden/Documents/Research/Manuscripts/X-TransPathway/Analysis/03_figures"

library("limma")
library("edgeR")
library("RColorBrewer")
library("gplots")
library("rtracklayer")
library("ggplot2")
library("gridExtra")
library("reshape")
library("pgirmess")
library("WGCNA")
library("vegan")
library("directlabels")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_palette <- colorRampPalette(c("blue", "white", "red3"))(n = 24)

##### Load data ##################################################################################################################

# Raw counts matrix
data_rnaseq = read.csv(paste(DataDir,"/gene_count_matrix.csv",sep=""))

# Sample information
pheno_data <- read.csv(paste(DataDir,"/phenotypes_pathway.csv",sep=""))

# Sample filtering, i.e. samples treated with low sugar + octopamine that are not included in the DGE analysis
samples_to_remove <- c("pathway_AR_LSOA_A","pathway_AU_LSOA_A","pathway_AU_LSOA_B","pathway_CA_LSOA_A")

# L. humile RefSeq
gene_set = read.csv(paste(DataDir,"/Lhum_RefSeq.csv",sep=""))

# Tables of candidate genes
ref_immune <- read.csv(paste(DataDir,"/lhum_immune_genes.csv",sep=""),na=TRUE,h=T); ref_immune <- unique(ref_immune[,c("LOC","Name","Pathway","Type","Comment")])

# creates DGElist object from gene count matrix (raw counts)
x_all <- DGEList(as.matrix(data_rnaseq[,c(2:length(data_rnaseq))]),remove.zeros=F)

##### Organising sample info #####################################################################################################

# Match the order of samples in pheno_data with the order of samples in data_rnaseq
pheno_data <- pheno_data[order(match(pheno_data$ids,colnames(data_rnaseq)[2:length(colnames(data_rnaseq))])),]

# Add sample information into x (DGEList)
x_all$samples$region <- pheno_data$region
x_all$genes <- data.frame("gene_id"=data_rnaseq[,1],"line"=seq(1:nrow(data_rnaseq)))

x <- x_all[,setdiff(colnames(data_rnaseq[,c(2:length(data_rnaseq))]),samples_to_remove),keep.lib.sizes=TRUE]

##### Organising gene info #######################################################################################################

# Import GTF file, turn it into a data frame and subset the transcripts only (no exons or tRNA)
gtf <- rtracklayer::import('~/Documents/Research/PhD/experiments/exp_RNA_pathway/New_Tuxedo/01_data/xxx20_simple_stringtie_for_DGE_755987/pathway_AR_A.gtf')
gtf_df=as.data.frame(gtf)
gtf_df_transcripts <- droplevels(subset(gtf_df, type == "transcript"))

# Get the unique matches between gene_id (MSTRG transcripts) and gene_name (LOC RefSeq numbers)
unique_gtf_df_transcripts <- unique(gtf_df_transcripts[,c("gene_id","ref_gene_name","gene_name")])

# Fill out ref_gene_name with LOC number when description not available (^^^)
unique_gtf_df_transcripts$ref_gene_name[is.na(unique_gtf_df_transcripts$ref_gene_name)] <- unique_gtf_df_transcripts[is.na(unique_gtf_df_transcripts$ref_gene_name),"gene_name"]
unique_gtf_df_transcripts$ref_gene_name <- as.factor(unique_gtf_df_transcripts$ref_gene_name)
unique_gtf_df_transcripts <- subset(unique_gtf_df_transcripts,select=c("gene_id","ref_gene_name"))

# Calculate gene-level length (i.e. average from isoforms)
mean_length <- aggregate(Gene...Transcripts...Length~Gene...Symbol, data=unique(gene_set[,c("Gene...Symbol","Gene...Transcripts...Length")]),FUN="mean")
# Retrive unique description per gene
unique_description <- unique(gene_set[,c("Gene...Symbol","Gene...Description")])
# Merge description and length
unique_description_length <- merge(unique_description,mean_length,by="Gene...Symbol")

# Merge into final object
all_genes_info <- na.omit(merge(unique_gtf_df_transcripts,unique_description_length,by.x = "ref_gene_name", by.y = "Gene...Symbol",all.x=TRUE))
summary(all_genes_info)

##### Data pre-processing #######################################################################################################

# Calculate CPMs before filtering and normalisation (i.e. not for the actual analysis)
cpm <- cpm(x)
log_cpm <- cpm(x,log=TRUE)
# Remove low-expressed transcripts
keep.exprs <- rowSums(cpm>1) >= 3 # Keep transcripts that are expressed at least 1 CPM in 3 samples
x_filt <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
dim(x_filt)
# Calculate CPMs post-filtering out low-expressed transcripts
cpm_filt <- cpm(x_filt)
log_cpm_filt <- cpm(x_filt,log=TRUE)

# Density plots before and after filtering
nsamples <- ncol(x)
colpal <- colorRampPalette(brewer.pal(8, "Set1"))
col <- colpal(nsamples)
par(mfrow=c(1,2))
plot(density(log_cpm[,1]), col=col[1], lwd=2, las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

plot(density(log_cpm_filt[,1]), col=col[1], lwd=2, las=2,
          main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_filt[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

dev.copy2pdf(file=paste(FigDir,"/sample_filtering.pdf",sep=""),
             width=12, 
             height=8)

##### Normalising gene expression distribution ###################################################################################

# Create object x2 (copy of x_filt) to first have a look at the effect of normalisation
x2 <- x_filt
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
cpm_x2 <- cpm(x2,log=TRUE)
boxplot(cpm_x2, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors
cpm_x2 <- cpm(x2,log=TRUE)
boxplot(cpm_x2, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Go ahead with TMM normalisation on x_filt
x_filt_tmm = calcNormFactors(x_filt, method='TMM')
cpm_filt_tmm <- cpm(x_filt_tmm,normalized.lib.sizes=TRUE)
rownames(cpm_filt_tmm) <- x_filt_tmm$genes$gene_id
log_cpm_filt_tmm <- cpm(x_filt_tmm,log=TRUE)

dev.copy2pdf(file=paste(FigDir,"/sample_normalisation.pdf",sep=""),
             width=12, 
             height=8)

##### Figure 1: Unsupervised clustering of samples ###############################################################################

# Define conditions (e.g region)
group <- as.factor(pheno_data[-which(pheno_data$ids %in% samples_to_remove),]$label)

# Plot sample clustering
par(mfrow=c(1,3),mar=c(4,4,4,4))
pheno_data$col.group <- ifelse(pheno_data$region=="AR","cornflowerblue",
                               ifelse(pheno_data$region=="CA","forestgreen",
                                      ifelse(pheno_data$region=="EU","gold2",
                                             ifelse(pheno_data$region=="AU","darkorange","red1"))))
col.group <- pheno_data[-which(pheno_data$ids %in% samples_to_remove),]$col.group

# labels
plotMDS(log_cpm_filt_tmm, col=col.group, label=group, dim=c(1,2), xlab = "MDS 1",ylab= "MDS 2",bty="l",xlim=c(-3.5,3.5))
plotMDS(log_cpm_filt_tmm, col=col.group, label=group, dim=c(3,2), xlab = "MDS 3",ylab= " ",bty="l",xlim=c(-2,3))

layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE))

# points
plotMDS(log_cpm_filt_tmm, col=col.group, pch=16, cex=2, cex.lab=2,cex.axis=2,dim=c(1,2), xlab = "MDS 1",ylab= "MDS 2",bty="l",xlim=c(-3.5,3.5))
plotMDS(log_cpm_filt_tmm, col=col.group, pch=16, cex=2, cex.lab=2,cex.axis=2, dim=c(3,2), xlab = "MDS 3",ylab= " ",bty="l",xlim=c(-2,3))

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Argentina', 'California', 'Europe',
                            'Australia', 'New Zealand'), pch=16, pt.cex=2, cex=2, bty='n',
       col = c('cornflowerblue', 'forestgreen', 'gold2', 'darkorange', 'red1'))
mtext("Region", at=0.2, cex=1.5)

dev.copy2pdf(file=paste(FigDir,"/F1_sample_clustering.pdf",sep=""),
             width=15, 
             height=5)

##### Differential Gene Expression analysis ######################################################################################

# Define experimental design
region <- x$sample$region
design = model.matrix(~0+region)

colnames(design) <- gsub("region","",colnames(design))
contr.matrix <- makeContrasts(ARvsCA = AR-CA,
                              ARvsEU = AR-EU,
                              ARvsAU = AR-AU,
                              ARvsNZ = AR-NZ,
                              levels = colnames(design))

# Prepare data for linear modelling
par(mfrow=c(1,2))
v = voom(x_filt_tmm, design, plot=T)
v$genes$gene_name <- all_genes_info[match(v$genes$gene_id,all_genes_info$gene_id),1]
v$genes$gene_descr <- all_genes_info[match(v$genes$gene_id,all_genes_info$gene_id),3]

# Fitting linear models
vfit <- lmFit(v,design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
plotSA(vfit, main="Final model")

# Number of DE genes (empirical Bayes stats)
tfit <- treat(vfit,lfc=log2(1.1))
dt <- (decideTests(tfit,p.value=0.05))
summary(dt)

par(mfrow=c(1,1))
de.common <- which(dt[,1]!=0&dt[,2]!=0&dt[,3]!=0&dt[,4]!=0)
length(de.common)
vennDiagram(dt,cex=c(2,2,2))
dev.copy2pdf(file=paste(FigDir,"/venn_diag.pdf",sep=""),
             width=10, 
             height=10)

tfit$genes$test <- "no"
tfit$genes$test[de.common] <- "yes"
de.common_names <- as.vector(tfit$genes$gene_id[de.common])
de.common_names <- as.vector(tfit$genes[de.common,])

as.vector(apply(tfit$coefficients[de.common,],1,mean))
DE_for_GO <- cbind(as.character(tfit$genes$gene_name[de.common]),as.vector(apply(tfit$coefficients[de.common,],1,mean)))
write.table(DE_for_GO,file=paste(DataDir,"/DE_for_GO_coeff.txt",sep=""),row.names = FALSE,col.names=FALSE)

##### DGE

# Figure 2. Volcano plots
ARvsCA = topTreat(tfit,coef=1,n=Inf)
ARvsEU = topTreat(tfit,coef=2,n=Inf)
ARvsAU = topTreat(tfit,coef=3,n=Inf)
ARvsNZ = topTreat(tfit,coef=4,n=Inf)

pairwise_comps <- list("AR versus CA" = ARvsCA,"AR versus EU" = ARvsEU,"AR versus AU" = ARvsAU,"AR versus NZ" = ARvsNZ)

plot_list_volcanos = list()
iteration = c(1:length(pairwise_comps))
for (comparison in iteration) {
  volcano_dat <- as.data.frame(pairwise_comps[[comparison]])
  title <- names(pairwise_comps)[comparison]
  volcano_dat$color <- "FDR > 0.05"
  volcano_dat$color[volcano_dat$adj.P.Val<0.05] <- "FDR < 0.05 & FC > 1.1"
  volcano_dat$color[which(volcano_dat$test == "yes")] <- "FDR < 0.05 & FC > 1.1 in all regions"
  
  volcano_plot = ggplot(aes(y=-1*log10(adj.P.Val),x=logFC,colour=color),data= volcano_dat) + theme_bw() + geom_point(size=2.5) +
    scale_color_manual(values=c("darkorange","red","black"),name="") +
    theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
    xlim(-12,12) + ylim(0,8) +
    labs(title=title,x="log2(FC)",y="-1*log10(FDR)") +
    theme(plot.title=element_text(size=18, vjust=2),legend.position="bottom", legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black"))
  if(comparison == 1) {legend <- get_legend(volcano_plot)}
  volcano_plot <- volcano_plot + theme(legend.position="none")
  plot_list_volcanos[[comparison]] = volcano_plot
}

grid.arrange(plot_list_volcanos[[1]],plot_list_volcanos[[2]],
             plot_list_volcanos[[3]],plot_list_volcanos[[4]],
             legend,ncol=2,nrow=3,layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),heights = c(rep(2.5,2), 0.25))

dev.copy2pdf(file=paste(FigDir,"/F2_volcanos.pdf",sep=""),
             width=9, 
             height=9)

##### DGE: Examining DE genes (heatmap) ##########################################################################################

# Figure 3. Heatmap
data_heatmap <- v$E[de.common,]

par(mfrow=c(1, 1))
heatmap_plot <- heatmap.2(v$E[de.common,],scale="row",col=my_palette,margins=c(8,24),key.title = "Gene expression",key.ylab = NA, key.ytickfun = 0,
          labRow=v$genes$gene_descr[de.common], labCol= rownames(v$targets), ColSideColors = c(rep("cornflowerblue",5),rep("forestgreen",6),rep("gold2",3),rep("darkorange",5),rep("red1",7)),
          distfun=function(x) as.dist(1-cor(t(x))), #looks good, but figure out how it is different from the euclidean distance (default)
          hclustfun=function(x) hclust(x, method="average"),
          keysize=1,cexRow=0.8,cexCol = 1.5,trace="none",density.info="none")

dev.copy2pdf(file=paste(FigDir,"/F3_heatmap_heads_region.pdf",sep=""),
             width=14, 
             height=20)

ordered_ids <- rownames(data_heatmap[rev(heatmap_plot$rowInd), heatmap_plot$colInd])
gene_names_heatmap_order <- all_genes_info[match(ordered_ids,all_genes_info$gene_id),]

##### DGE: Examining candidate genes #############################################################################################

### FPKM FROM DGEList
x_filt_tmm$genes$length <- all_genes_info$Gene...Transcripts...Length[match(x_filt_tmm$genes$gene_id,all_genes_info$gene_id)]
x_filt_tmm$genes$ref_gene_name <- all_genes_info[match(x_filt_tmm$genes$gene_id,all_genes_info$gene_id),1]
rpkm_all <- data.frame(rpkm(x_filt_tmm,normalized.lib.sizes=TRUE,gene.length=x_filt_tmm$genes$length,log=TRUE))
rownames(rpkm_all) <- x_filt_tmm$genes$ref_gene_name
t_data_num_all <- data.frame(scale(t(rpkm_all),center=TRUE,scale=TRUE),check.names=FALSE)

################## start loop based on candidate gene subset here
categories <- unique(ref_immune$Pathway)
plot_gene_region_list = list()
plot_gene_sample_list = list()
for (category in unique(categories)) {
  reference_subset <- ref_immune[ref_immune$Pathway == category,]
  reference_subset$LOC <- gsub("LOC","",reference_subset$LOC)
  if(nrow(reference_subset) <2) next
  t_data_num_subset <- t_data_num_all[,colnames(t_data_num_all) %in% paste("LOC",reference_subset$LOC,sep="")]
  t_data_num_subset_melted <- melt(as.matrix(t_data_num_subset))
  if(length(t_data_num_subset) == nrow(t_data_num_all)) next
  t_data_num_subset_melted$region <- unlist(lapply(rownames(t_data_num_subset),function(x) ifelse(substr(x,9,10) == "AR","AR",
                                                                                           ifelse(substr(x,9,10) == "CA","CA",
                                                                                                  ifelse(substr(x,9,10) == "EU","EU",
                                                                                                         ifelse(substr(x,9,10) == "AU","AU", "NZ"))))))

  names(t_data_num_subset_melted) <- c("sample","gene_name","value","region")
  t_data_num_subset_melted_names <- merge(t_data_num_subset_melted,reference_subset,by.x="gene_name",by.y="LOC",all.x=TRUE)[,c(1:6)]
  
  results_gene_level=NULL
  for (gene in unique(t_data_num_subset_melted$gene_name)) {
    subset_gene <- subset(t_data_num_subset_melted_names,t_data_num_subset_melted_names$gene_name == gene)
    short_name <- reference_subset[match(gene,reference_subset$LOC),2]
    gene_label <- paste(category,gene,"-",short_name,sep=" ")
    print(gene_label)
    mean_values_gene <- aggregate(subset_gene$value,by=list(subset_gene$region),function(x) c(mean = mean(x), sd = sd(x), se = sd(x)/sqrt(length(x))))
    gene_result <- data.frame("LOC" = rep(gene,nrow(mean_values_gene)),
                              "Pathway" = category,
                              "region" = mean_values_gene$Group.1,
                              "value" = mean_values_gene$x[,1],
                              "sd" = mean_values_gene$x[,2],
                              "se" = mean_values_gene$x[,3])
    gene_result_not_reordered <- gene_result
    gene_result$region <- factor(gene_result$region,levels=c("AR","CA","EU","AU","NZ"),ordered=T)
    
    lim_y <- max(abs(gene_result$value))
    
    plot_gene_region <- ggplot(data = gene_result, aes(x = region, y = value)) +
      geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("cornflowerblue","forestgreen","gold2","darkorange","red1")) +
      geom_errorbar(aes(ymax = gene_result$value + gene_result$se, ymin = gene_result$value - gene_result$se),
                    position = position_dodge(0.90), width = 0.25) +
      theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) + ylim(-lim_y-0.6*lim_y,lim_y+0.6*lim_y) +
      labs(title=gene_label,x="",y="Transcript expression (log2-centred\nTMM-normalised RPKM)") +
      theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
            axis.text.x = element_text(size = 14, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.title.y=element_text(size = 16, colour = "black",vjust=1),
            axis.title.x=element_text(size = 16, colour = "black"))
    plot_gene_region_list[[gene]] <- plot_gene_region
    
    subset_gene$shsample <- gsub("pathway_","",subset_gene$sample)
    subset_gene$shsample <- factor(subset_gene$shsample, levels= c("AR_A","AR_B","AR_C","AR_D","AR_LS_A",
                                                                  "CA_A","CA_B","CA_C","CA_D","CA_E","CA_LS_A",
                                                                  "EU_A","EU_B","EU_C",
                                                                  "AU_A","AU_B","AU_C","AU_D","AU_LS_A",
                                                                  "NZ_A","NZ_B","NZ_C","NZ_D","NZ_E","NZ_LS_A","NZ_LS_B"), ordered=T)
    write.csv(subset_gene,file=paste(FigDir,"/individual_genes_detailed/data/",category,"_",gene,".csv",sep=""))
    
    lim_y <- max(abs(subset_gene$value))
    
    plot_gene_sample <- ggplot(data = subset_gene, aes(x = shsample, y = value)) +
      geom_bar(stat = "identity", position = position_dodge(0.90),fill=c(rep("cornflowerblue",5),rep("forestgreen",6),rep("gold2",3),rep("darkorange",5),rep("red1",7))) +
      theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) + ylim(-lim_y-0.1*lim_y,lim_y+0.1*lim_y) +
      labs(title=gene_label,x="",y="Transcript expression\n(log2-centered FPKMs)") +
      theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
            axis.text.x = element_text(size = 14, colour = "black",angle=90),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.title.y=element_text(size = 16, colour = "black",vjust=1),
            axis.title.x=element_text(size = 16, colour = "black"))
    plot_gene_sample_list[[gene]] <- plot_gene_sample
    
    grid.arrange(plot_gene_region,plot_gene_sample,ncol=1,nrow=2)
    dev.copy2pdf(file=paste(FigDir,"/individual_genes_detailed/",category,"_",gene,".pdf",sep=""),width=8, height=12)
    results_gene_level=rbind(gene_result_not_reordered,results_gene_level)
  }
  write.csv(results_gene_level,file=paste(FigDir,"/individual_genes_detailed/data/",category,".csv",sep=""))
}

# At the category level

plot_list_categ = list()
for (categ in unique(categories)) {
  print(categ)
  if(any(list.files(paste(FigDir,"/individual_genes_detailed/data/",sep="")) == paste(categ,".csv",sep="")) ==  FALSE) next
  result_gene_level <- read.csv(paste(FigDir,"/individual_genes_detailed/data/",categ,".csv",sep=""),h=TRUE)
  subset_category_mean <- subset(result_gene_level,result_gene_level$Pathway==categ)
  print(summary(subset_category_mean))
  if (length(unique(subset_category_mean$LOC)) > 10) {
  glm_data <- glm(value~region,data=subset_category_mean)
  summary <- summary(glm_data)
  print(summary)
  
  p_vals_not_ordered <- ifelse(round(as.numeric(summary$coefficients[,4]),4)<0.001,"p < 0.001",
                               ifelse(round(as.numeric(summary$coefficients[,4]),4)<0.01,"p < 0.01",
                                      ifelse(round(as.numeric(summary$coefficients[,4]),4)<0.05,"p < 0.05",
                                             paste("p = ",round(as.numeric(summary$coefficients[,4]),2),sep=""))))
  p_vals <- p_vals_not_ordered[c(3,4,2,5)] # Fix that, no reordering please
  } else {
    kruskal_data <- kruskalmc(value~region,data=subset_category_mean,cont="two-tailed")
    print(kruskal_data)
    
    p_vals_not_ordered <- ifelse(kruskal_data$dif.com$difference == TRUE, "p < 0.05","NS")
    p_vals <- p_vals_not_ordered[c(2,3,1,4)] # Fix that, no reordering please 
  }
  y_min = round_any(min(subset_category_mean$value),.5,f=floor)
  y_max = round_any(max(subset_category_mean$value),.5,f=ceiling)
  lim_y=max(c(abs(y_min),abs(y_max)))
  
  boxplot = ggplot(aes(y=value,x=region,fill=region),data= subset_category_mean) + theme_bw() +
    scale_fill_manual(values=c("cornflowerblue","darkorange", "forestgreen","gold2","red1"),name="") +
    scale_x_discrete(limits=c("AR","CA","EU","AU","NZ")) +
    theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
    labs(title=categ,x="",y="Transcript expression (log2-centred\nTMM-normalised FPKM)") +
    theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          axis.title.y=element_text(size = 18, colour = "black",vjust=1),
          axis.title.x=element_text(size = 18, colour = "black")) +
    annotate("text", x = as.factor(c("CA","EU","AU","NZ")),y=rep(lim_y+0.5*lim_y,4),label = p_vals,family="",fontface = 3,size=5) +
    ylim(-(lim_y+0.5*lim_y),lim_y+0.5*lim_y) + stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=F,outlier.shape=NA) +
    geom_point(size=2,position = position_jitter(width = 0.2)) + stat_summary(fun.y=mean, colour = "white",geom="point", size=4)
  print(boxplot)
  plot_list_categ[[categ]] = boxplot
}

grid.arrange(plot_list_categ$RNAi,plot_list_categ$'JAK-STAT',plot_list_categ$Imd,plot_list_categ$JNK,plot_list_categ$TOLL,ncol=2,nrow=3)
dev.copy2pdf(file=paste(FigDir,"/F5_panel_immune.pdf",sep=""), width=12.8, height=14.4)

##### Biogenic amine genes ###################################################################################################

bioamines <- read.csv(paste(DataDir,"/lhum_bioamine_genes.csv",sep=""))
t_rpkm_bioamines <- t_data_num_all[,colnames(t_data_num_all) %in% paste("LOC",bioamines$gene_name,sep="")]

result_krusk_range = NULL
for (gene in colnames(t_rpkm_bioamines)) {
  subset_gene <- t_rpkm_bioamines[,gene,drop=F]
  subset_gene$region <- substr(rownames(subset_gene),9,10)
  subset_gene$range <- ifelse(subset_gene$region == "AR","native","introduced")
  gene_label <- paste(bioamines[bioamines$gene_name==gene,"short_name"]," - LOC",gene,sep="")
  print(gene_label)
  print(subset_gene)
  mean_values_gene <- aggregate(subset_gene[,1],by=list(subset_gene$region,subset_gene$range),function(x) c(mean = mean(x), sd = sd(x), se = sd(x)/sqrt(length(x))))
  colnames(mean_values_gene) <- c("region","range","gene")
  mean_values_gene$range <- factor(mean_values_gene$range,levels=c("native","introduced"),ordered=TRUE)
  mean_values_gene$region <- factor(mean_values_gene$region,levels=c("AR","CA","EU","AU","NZ"),ordered=T)
  
  ## pairwise kruskal tests
  kruskal_pwse <- kruskalmc(subset_gene[,1]~region,data=subset_gene,cont="two-tailed")
  print(kruskal_pwse)
  significance_pwse <- ifelse(kruskal_pwse$dif.com$difference == "TRUE", "p < 0.05","NS")
  
  ## kruskal test within invaded range
  subset_gene_inv <- subset_gene[subset_gene$range=="introduced",]
  kruskal_inv <- kruskal.test(subset_gene_inv[,1]~as.factor(region),subset_gene_inv)
  print(kruskal_inv)
  significance_inv <- ifelse(kruskal_inv$p.value < 0.05, "p < 0.05", paste("NS, p = ",round(kruskal_inv$p.value,3),sep=""))
  
  ## kruskal test on range for next section (if no difference within invaded range)
  kruskal_range <- kruskal.test(subset_gene[,1]~as.factor(range),subset_gene)
  significance_range <- data.frame("significance" = ifelse(kruskal_range$p.value < 0.001, paste("KW = ",round(kruskal_range$statistic,2),", p < 0.001",sep=""),
                                                           ifelse(kruskal_range$p.value < 0.01, paste("KW = ",round(kruskal_range$statistic,2),", p < 0.01",sep=""),
                                                                  ifelse(kruskal_range$p.value < 0.05, paste("KW = ",round(kruskal_range$statistic,2),", p < 0.05",sep=""), paste("NS, p = ",round(kruskal_range$p.value,3),sep="")))),"gene_name"=gene)
  print(kruskal_range)
  result_krusk_range <- rbind(result_krusk_range, significance_range)
  
  lim_y <- round(max(abs(mean_values_gene$gene[,1])+mean_values_gene$gene[,3]))

  # plot if no differnce within invasive range
if (kruskal_inv$p.value > 0.05) {  
    plot_gene_region <- ggplot(data = mean_values_gene, aes(x = region, y = gene[,1])) +
    geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("cornflowerblue","forestgreen","gold2","darkorange","red1")) +
    geom_errorbar(aes(ymax = mean_values_gene$gene[,1] + mean_values_gene$gene[,3], ymin = mean_values_gene$gene[,1] - mean_values_gene$gene[,3]),
                  position = position_dodge(0.90), width = 0.25) +
    theme_bw() + ylim(-3,4) +
    labs(title=gene_label,x="",y="Transcript expression (log2-centered\nTMM-normalised FPKMs)") +
    theme(plot.title=element_text(size=16, vjust=2, hjust=0),legend.position="", legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black")) +
    geom_segment(aes(x=1,xend=3.5,y=3,yend=3)) +
    annotate("text", x = c(2.25),
                      y=c(3.25),
                      label = c(as.character(significance_range$significance)),fontface = 4,size=5) +
    geom_segment(aes(x=2,xend=5,y=2.25,yend=2.25)) +
      annotate("text", x = c(3.5),
               y=c(2.5),
               label = c(significance_inv),fontface = 3,size=4)
} else {
  # plot if there are differences within invasive range
    plot_gene_region <- ggplot(data = mean_values_gene, aes(x = region, y = gene[,1])) +
    geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("cornflowerblue","forestgreen","gold2","darkorange","red1")) +
    geom_errorbar(aes(ymax = mean_values_gene$gene[,1] + mean_values_gene$gene[,3], ymin = mean_values_gene$gene[,1] - mean_values_gene$gene[,3]),
                  position = position_dodge(0.90), width = 0.25) +
      theme_bw() + ylim(-3,4) +
    labs(title=gene_label,x="",y="Transcript expression (log2-centered\nTMM-normalised FPKMs)") +
    theme(plot.title=element_text(size=16, vjust=2, hjust=0),legend.position="", legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black")) +
      geom_segment(aes(x=2,xend=5,y=2.25,yend=2.25)) +
      annotate("text", x = c(3.5),
               y=c(2.5),
               label = c(significance_inv),fontface = 3,size=4) +
      annotate("text", x = c(2,3,4,5),
               y=rep(c(lim_y+((lim_y/10)*2)*3.5),4),
               label = c(significance_pwse[c(2,3,1,4)]),fontface = 4,size=5)
}
  print(plot_gene_region)
  dev.copy2pdf(file=paste(FigDir,"/individual_genes_detailed/bioamine_",gene,".pdf",sep=""),width=6, height=6)
}

t_rpkm_bioamines_for_pooling <- t_rpkm_bioamines
t_rpkm_bioamines_for_pooling$region <- substr(rownames(t_rpkm_bioamines_for_pooling),9,10)
t_rpkm_bioamines_for_pooling$range <- ifelse(t_rpkm_bioamines_for_pooling$region == "AR","native","introduced")
melted <- melt(t_rpkm_bioamines_for_pooling)
melted_mean <- aggregate(melted$value,by=list("range"=melted$range,"gene_name"=melted$variable),mean)
melted_mean$short_name <- bioamines[match(melted_mean$gene_name,paste("LOC",bioamines$gene_name,sep="")),2]
melted_mean$range <- factor(melted_mean$range,levels=c("native","introduced"),ordered=TRUE)

# create facet labels with gene name and test result
result_krusk_range$short_name <- bioamines[match(result_krusk_range$gene_name,paste("LOC",bioamines$gene_name,sep="")),2]
facet_labels=NULL
for (iteration in seq_along(1:nrow(result_krusk_range))) {
  result <- paste(paste(result_krusk_range[iteration,3],result_krusk_range[iteration,2],sep=" "),result_krusk_range[iteration,1],sep="\n")
  facet_labels <- c(facet_labels,gene_ref=result)
}
names(facet_labels) <- result_krusk_range$short_name

melted_mean$short_name <- factor(melted_mean$short_name, levels=c("5-HTR1","5-HTR2A","5-HTR2A-like","5-HTR2C",
                                                                  "octB1R","octB2R","octB3R","octR1",
                                                                  "dopD2-like","dopR1","dopR2",
                                                                  "TyrBH","TyrR1","TyrR1-like"))

facet_genes <- ggplot(data = melted_mean, aes(x = range, y = x)) +
  geom_bar(stat="identity", position = position_dodge(0),fill=rep(c("cornflowerblue","red1"),length(unique(melted_mean$short_name)))) +
  facet_wrap(~factor(short_name), labeller = as_labeller(facet_labels)) +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) + ylim(-lim_y-0.6*lim_y,lim_y+0.6*lim_y) +
  labs(title="",x="",y="Transcript expression (log2-centered\nTMM-normalised RPKMs)") +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
facet_genes
dev.copy2pdf(file=paste(FigDir,"/F4_facet_plot_amines.pdf",sep=""),width=10, height=12)

##### Viral loads ################################################################################################################

BLASTx_all <- rbind(read.csv(paste(DataDir,"/viruses_pathway_refseq.blastx.outfmt6_clean.csv",sep="")),
                    read.csv(paste(DataDir,"/viruses_pathway_viljakainen.blastx.outfmt6_clean.csv",sep="")))
BLASTx_all$sseqid <- gsub("[|]","",gsub(".*gb[|]","",BLASTx_all$sseqid))
BLASTx_all <- BLASTx_all[-match("YP_009315906.1",BLASTx_all$sseqid),] # Remove irreleveant King virus (i.e. mammal virus)

expression_matrix <- read.table(paste(DataDir,"/genes.TMM.EXPR.matrix",sep=""))
expression_matrix <- expression_matrix[-which(rowSums(expression_matrix)==0),]

# Use library sizes of reads mapped to the Argentine ant genome to normalise viral counts -> host-library-size-standardised TMM-normalised TPMs
mapped_reads_lib_size <- data.frame("sample" = rownames(x_all$samples), "alias" = colnames(expression_matrix), "lib.size" = x_all$samples$lib.size)
for(n in seq(1:ncol(expression_matrix))){
  expression_matrix[,n] <- (expression_matrix[,n]/mapped_reads_lib_size[n,3])*(mean(mapped_reads_lib_size[,3]))
}

expression_matrix$gene_id <- rownames(expression_matrix)

virus_hits <- read.csv(paste(DataDir,"/virus_hits_pathway.csv",sep=""))

list_TOI_BLASTx <- unique(BLASTx_all[,1:2])
TOI_expr_BLASTx <- merge(list_TOI_BLASTx,expression_matrix,by.x="qseqid",by.y="gene_id")[,-1]
TOI_expr_BLASTx <- aggregate(TOI_expr_BLASTx[,-1],by=list(sseqid=TOI_expr_BLASTx$sseqid),sum)

TOI_expr_BLASTx$virus <- virus_hits$organism[match(TOI_expr_BLASTx$sseqid,virus_hits$accession)]
rownames(TOI_expr_BLASTx) <- TOI_expr_BLASTx$sseqid; TOI_expr_BLASTx <- TOI_expr_BLASTx[,-1]

t_data_num_viruses <- data.frame(scale(t(TOI_expr_BLASTx[,-ncol(TOI_expr_BLASTx)]),center = FALSE,scale=TRUE))

plot_microbe_gene_site_list = list()
plot_microbe_gene_sample_list = list()
for (microbe in unique(TOI_expr_BLASTx$virus)) {
  viral_features <- subset(virus_hits,organism == microbe)
  subset_virus <- data.frame(t_data_num_viruses[,colnames(t_data_num_viruses) %in% viral_features$accession])
  rownames(subset_virus) <- rownames(t_data_num_viruses); colnames(subset_virus) <- colnames(t_data_num_viruses)[colnames(t_data_num_viruses) %in% viral_features$accession]
    # remove viral transcripts that are not expressed in the focus samples (i.e. heads)
    subset_virus <- subset_virus[sapply(subset_virus, function(x) !any(is.na(x)))]
    if (ncol(subset_virus) == 0) next
  subset_virus$sample <- rownames(subset_virus)
  subset_virus_melted <- melt(subset_virus,id="sample")
    # Calculate viral loads from all transcripts
    sum_viral_transcripts <- aggregate(subset_virus_melted$value,by=list(subset_virus_melted$sample),sum)
    colnames(sum_viral_transcripts) <- c("sample",microbe)
    sum_viral_transcripts$region <- substr(sum_viral_transcripts$sample,1,2)
    write.csv(sum_viral_transcripts,file=paste(FigDir,"/viruses/data/sum_trans_",microbe,".csv",sep=""),row.names=FALSE)
  results_microbe_gene_level=NULL
  for (transcript in unique(subset_virus_melted$variable)){
    vir_label <- paste(microbe," - ", transcript, sep="")
    subset_virus_transcript <- subset(subset_virus_melted, variable == transcript)
    if(is.na(subset_virus_transcript$value)==TRUE) next
    subset_virus_transcript$region <- substr(subset_virus_transcript$sample,1,2)
    subset_virus_transcript$sample <- factor(subset_virus_transcript$sample, levels = c("AR_A","AR_B","AR_C","AR_D","AR_LS_A","AR_LSOA_A",
                                                                                        "CA_A","CA_B","CA_C","CA_D","CA_E","CA_LS_A","CA_LSOA_A",
                                                                                        "EU_A","EU_B","EU_C",
                                                                                        "AU_A","AU_B","AU_C","AU_D","AU_LS_A","AU_LSOA_A","AU_LSOA_B",
                                                                                        "NZ_A","NZ_B","NZ_C","NZ_D","NZ_E","NZ_LS_A","NZ_LS_B"),ordered=TRUE)
    print(transcript)
    test_virus_transcript <- kruskalmc(value~region,data=subset_virus_transcript,cont="two-tailed")
    print(test_virus_transcript) 
    significance <- ifelse(test_virus_transcript$dif.com$difference == "TRUE", "p < 0.05","NS")
    ##
    mean_values_microbe_gene <- aggregate(subset_virus_transcript$value,by=list(subset_virus_transcript$region),function(x) c(mean = mean(x), sd = sd(x), se = sd(x)/sqrt(length(x))))
    microbe_gene_result <- data.frame("Transcript" = rep(transcript,nrow(mean_values_microbe_gene)),
                                      "Microbe" = microbe,
                                      "region" = mean_values_microbe_gene$Group.1,
                                      "value" = mean_values_microbe_gene$x[,1],
                                      "sd" = mean_values_microbe_gene$x[,2],
                                      "se" = mean_values_microbe_gene$x[,3])
    microbe_gene_result$region <- factor(microbe_gene_result$region, levels = c("AR","CA","EU","AU","NZ"))
    
    lim_y=round_any(max(abs(microbe_gene_result$value))+max(microbe_gene_result$se),1,f=ceiling)
    
    plot_gene_site <- ggplot(data = microbe_gene_result, aes(x = region, y = value)) +
      geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("cornflowerblue","forestgreen","gold2","darkorange","red1")) +
      geom_errorbar(aes(ymax = microbe_gene_result$value + microbe_gene_result$se, ymin = microbe_gene_result$value - microbe_gene_result$se),
                    position = position_dodge(0.90), width = 0.25) +
      theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) + ylim(0,lim_y+0.6*lim_y) +
      labs(title=vir_label,x="",y="Transcript expression\n(scaled TMM-normalised TPMs)") +
      theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
            axis.text.x = element_text(size = 14, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.title.y=element_text(size = 16, colour = "black",vjust=1),
            axis.title.x=element_text(size = 16, colour = "black")) +
      annotate("text",
               x = c(2, 3, 4,5),
               y = rep(lim_y+0.1*lim_y,4),
               label = c(significance[2],significance[3],significance[1],significance[4]),family="",fontface = 4,size=4,
               family = "", fontface = 3, size=4)
    plot_microbe_gene_site_list[[transcript]] <- plot_gene_site
    
    ##
    plot_gene_sample <- ggplot(data = subset_virus_transcript, aes(x = sample, y = value)) +
      geom_bar(stat = "identity", position = position_dodge(0.90),fill=c(rep("cornflowerblue",6),rep("forestgreen",7),rep("gold2",3),rep("darkorange",7),rep("red1",7))) +
      theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
      labs(title=vir_label,x="",y="Transcript expression\n(scaled TMM-normalised TPMs)") +
      theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
            axis.text.x = element_text(size = 14, colour = "black",angle=90),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.title.y=element_text(size = 16, colour = "black",vjust=1),
            axis.title.x=element_text(size = 16, colour = "black"))
    plot_microbe_gene_sample_list[[transcript]] <- plot_gene_sample
    
    grid.arrange(plot_gene_site,plot_gene_sample,ncol=1,nrow=2)    
    dev.copy2pdf(file=paste(FigDir,"/viruses/",microbe, " - ",transcript,".pdf",sep=""),
                 width=8, 
                 height=12)
    
    results_microbe_gene_level=rbind(results_microbe_gene_level,microbe_gene_result)
  }
  write.csv(results_microbe_gene_level,file=paste(FigDir,"/viruses/data",microbe,".csv",sep=""),row.names=FALSE)
}

# Viral loads from all viral transcripts (per virus)

data_overall_viral_loads = NULL
for (microbe_file in unique(list.files(paste(FigDir,"/viruses/data",sep=""),pattern="sum_trans"))) {
  viral_transcripts <- read.csv(paste(FigDir,"/viruses/data","/",microbe_file,sep=""))
  plot_virus_sample <- ggplot(data = viral_transcripts, aes(x = viral_transcripts$sample, y = viral_transcripts[,2])) +
    geom_segment(aes(y=0, yend=viral_transcripts[,2], x=viral_transcripts$sample, xend=viral_transcripts$sample), size=10) +
        theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
    labs(title=colnames(viral_transcripts)[2],x="",y="Transcript expression\n(scaled TMM-normalised TPMs)") +
    theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14, colour = "black",angle=90,vjust=0.5),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black"))
  print(plot_virus_sample)
  dev.copy2pdf(file=paste(FigDir,"/viruses/sum_trans_",colnames(viral_transcripts)[2],".pdf",sep=""), width=8, height=6)
  # Calculate overall viral loads
  viral_transcripts$virus <- colnames(viral_transcripts)[2]; colnames(viral_transcripts)[2] <- "value"
  data_overall_viral_loads <- rbind(data_overall_viral_loads,viral_transcripts)
}

# Figure 6: Overall viral loads

# Synthesis plot: viral loads and diversity
data_overall_viral_loads <- na.omit(data_overall_viral_loads)
data_overall_viral_loads$sample <- factor(data_overall_viral_loads$sample[c(1:6,14:20,21:23,7:13,24:30)])
plot_virus_load_div <- ggplot(data = data_overall_viral_loads, aes(x = data_overall_viral_loads$sample, y = data_overall_viral_loads$value,fill=data_overall_viral_loads$virus)) +
  geom_bar(stat = "identity") +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  labs(title="Overall viral loads",x="",y="Overall viral load\n(scaled TMM-normalised TPMs)") +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="right", legend.text=element_text(size=12),#legend.title=element_text("Viruses"),
        axis.text.x = element_text(size = 14, colour = "black",angle=90,vjust=0.5),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
print(plot_virus_load_div)
dev.copy2pdf(file=paste(FigDir,"/F6_overall_viral_loads_div.pdf",sep=""),
             width=10, 
             height=6)

### Viral diversity, presence/absence, species richness
data_overall_viral_loads$presence <- ifelse(data_overall_viral_loads$value > 0, 1,0)
viral_presence_clean=data.frame("sample"=unique(data_overall_viral_loads$sample))
for (microbe in unique(data_overall_viral_loads$virus)) {
  subset_microbe <- subset(data_overall_viral_loads, data_overall_viral_loads$virus == microbe)
  data_microbe <- data.frame("microbe" = subset_microbe$presence); colnames(data_microbe) <- microbe
  viral_presence_clean <- cbind(viral_presence_clean,data_microbe)
}
viral_presence_clean <- cbind("region"=subset_microbe[,"region"],viral_presence_clean)

# Venn diagram
viral_presence_venn <- t(aggregate(viral_presence_clean[,3:length(colnames(viral_presence_clean))],by=list(viral_presence_clean$region),max)[,-1])
colnames(viral_presence_venn) <- unique(viral_presence_clean$region)
vennDiagram(viral_presence_venn,cex=c(2,2,1),main=paste("Viral diversity, virus present when TPM > 0",sep=""))
write.csv(viral_presence_venn,file=paste(FigDir,"/viruses/data/venn_diag_viruses",".csv",sep=""),row.names = TRUE)
dev.copy2pdf(file=paste(FigDir,"/F6_venn_diag_viruses.pdf",sep=""), width=8, height=8)

# Accumulation curves
viruses_AR <- subset(viral_presence_clean, region == "AR")
viruses_CA <- subset(viral_presence_clean, region == "CA")
viruses_EU <- subset(viral_presence_clean, region == "EU")
viruses_AU <- subset(viral_presence_clean, region == "AU")
viruses_NZ <- subset(viral_presence_clean, region == "NZ")

fish_nz <- fisher.alpha(viruses_NZ[,3:ncol(viruses_NZ)],MARGIN=1)

accum_AR <- specaccum(viruses_AR[,3:ncol(viruses_AR)]); accum_AR_fmt <- data.frame("region" = rep("AR",7), "site" = seq(1:7), "richness" = as.numeric(c(accum_AR$richness,"NA")), "sd" = as.numeric(c(accum_AR$sd,"NA")))
accum_CA <- specaccum(viruses_CA[,3:ncol(viruses_CA)]); accum_CA_fmt <- data.frame("region" = rep("CA",7), "site" = seq(1:7), "richness" = as.numeric(c(accum_CA$richness)), "sd" = as.numeric(c(accum_CA$sd)))
accum_EU <- specaccum(viruses_EU[,3:ncol(viruses_EU)]); accum_EU_fmt <- data.frame("region" = rep("EU",7), "site" = seq(1:7), "richness" = as.numeric(c(accum_EU$richness,rep("NA",4))), "sd" = as.numeric(c(accum_EU$sd,rep("NA",4))))
accum_AU <- specaccum(viruses_AU[,3:ncol(viruses_AU)]); accum_AU_fmt <- data.frame("region" = rep("AU",7), "site" = seq(1:7), "richness" = as.numeric(c(accum_AU$richness)), "sd" = as.numeric(c(accum_AU$sd)))
accum_NZ <- specaccum(viruses_NZ[,3:ncol(viruses_NZ)]); accum_NZ_fmt <- data.frame("region" = rep("NZ",7), "site" = seq(1:7), "richness" = accum_NZ$richness, "sd" = accum_NZ$sd)

accum_all <- rbind(accum_AR_fmt,accum_CA_fmt,accum_EU_fmt,accum_AU_fmt,accum_NZ_fmt); accum_all$sd[accum_all$sd == 0] <- NA

plot_accum <- ggplot(data = accum_all, aes(x = site, y = richness, fill = region)) +
  geom_line(colour=c(rep("cornflowerblue",7),rep("darkorange",7),rep("forestgreen",7),rep("gold2",7),rep("red1",7))) +
  geom_errorbar(aes(ymax = richness+sd, ymin = richness-sd), width = 0.25, position = position_dodge(0.25)) +
  geom_dl(aes(label = region), method = list(dl.trans(x = x + 0.5),"last.points"), cex = 1) +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18)) +
  labs(title="Virus accumulation curve",x="Libraries",y="Virus richness") +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="", legend.text=element_text(size=14),
        axis.text.x = element_text(size = 14, colour = "black",angle=90,vjust=0.5),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
plot_accum
dev.copy2pdf(file=paste(FigDir,"/F6_accumulation_curves.pdf",sep=""), width=8, height=8)

# Table 1: viral diversity
regional_virus_presence <- aggregate(viral_presence_clean[,3:ncol(viral_presence_clean)],by=list(viral_presence_clean$region),max)
t(regional_virus_presence)
regional_virus_presence$Richness <- apply(regional_virus_presence[,2:ncol(regional_virus_presence)],1,sum)


### WGCNA

###### SECTION 1: PREPARE THE DATA

options(stringsAsFactors = FALSE);
allowWGCNAThreads()

#Read in the data set
femData = cpm_filt_tmm
colnames(femData) <- gsub("pathway_","",colnames(femData))

traitData = read.csv(paste(DataDir, "/WGCNA_sample_description.csv",sep=""));
traitData$sample <- traitData$short_name

# filter out samples that are not wanted, and remove genes not expressed in more than 90% of the samples
femData <- femData[apply(femData == 0,1,sum) <= (0.9*ncol(femData)),-ncol(femData)]

# transpose the data after a variance-stabilising log2 transformation
datExpr0 = as.data.frame(t(log(femData+1)));

# Quality checks and filter
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK # if FALSE, that means there are some data point to remove

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Check sample clustering and potentially remove outliers
sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 140, col = "red")
# Determine luster under the line
clust = cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
allTraits = traitData[match(Samples, traitData$sample),];

allTraits$range_num <- ifelse(allTraits$range == "introduced", 1, 0)
datTraits <- allTraits$range_num

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_sample_dendro.pdf",sep=""), width=8, height=6)

###### SECTION 2b: STEP-by-STEP ANALYSIS

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_sft_plots.pdf",sep=""), width=8, height=5)

# Set soft thresholding power
softPower = 3;
adjacency = adjacency(datExpr, power = softPower, type = "signed hybrid");

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_cut_tree.pdf",sep=""), width=8, height=6)
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_cluster_dendro.pdf",sep=""), width=8, height=6)



# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

###### SECTION 3: Relating modules to external traits

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 8, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,# xLabels= "Range",
               xLabels = c("range"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_trait_heatmap.pdf",sep=""), width=6, height=8)

# Define variable range_num containing the range_num column of datTrait
range_num = as.data.frame(datTraits)
names(range_num) = "range_num"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

rownames(geneModuleMembership) <- all_genes_info[match(rownames(geneModuleMembership),all_genes_info$gene_id),1]
head(geneModuleMembership)


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, range_num, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(range_num), sep="");
names(GSPvalue) = paste("p.GS.", names(range_num), sep="");

rownames(geneTraitSignificance) <- all_genes_info[match(rownames(geneTraitSignificance),all_genes_info$gene_id),1]
head(geneTraitSignificance)

###

module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for range_num",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, pch = 16, col = module)
rect(0.7, 0.8, 1, 1, col=NULL, border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_module_of_interest.pdf",sep=""), width=10, height=8)

gene_module_membership <- data.frame(rownames(geneModuleMembership[moduleGenes,]), geneModuleMembership[moduleGenes,column])
gene_module_sign <- data.frame("gene_id"=rownames(geneTraitSignificance[moduleGenes,,drop=F]), "significance"=geneTraitSignificance[moduleGenes,])
write.table(gene_module_membership,paste(FigDir,"/WGCNA/red_new_module_for_GO.txt",sep=""),row.names = FALSE, col.names = FALSE)

top_genes_mod <- all_genes_info[match(gene_module_sign[which(abs(gene_module_sign$significance) > 0.8),"gene_id"],all_genes_info$ref_gene_name),c("Gene...Description","ref_gene_name")]
write.csv(top_genes_mod,paste(FigDir,"/WGCNA/top_genes_mod.csv",sep=""))


###

for (modcolor in unique(moduleColors)) {
print(paste(modcolor,": ",length(names(datExpr)[moduleColors==modcolor]),sep=""))
  list_genes_mod <- names(datExpr)[moduleColors==modcolor]
  list_acc_mod <- unique_gtf_df_transcripts[match(list_genes_mod,unique_gtf_df_transcripts$gene_id),]
assign(paste("mod_",modcolor,sep=""), list_acc_mod)
  }

for (modcolor in unique(moduleColors)) {
  moi_ids <- names(datExpr)[moduleColors==modcolor]
  module_of_interest <- all_genes_info[match(moi_ids,all_genes_info$gene_id),]
  write.csv(module_of_interest,paste(FigDir,"/WGCNA/genes_in_",modcolor,"_module.csv",sep=""))}

genes_of_interest <- na.omit(ref_immune[match(mod_red$ref_gene_name,factor(paste("LOC",ref_immune$LOC,sep=""))),])

### Prepare file S1

DGE_genes <- de.common_names[,3:4]

red_module <- data.frame("gene_name"=mod_red[,2], "gene_descr"= all_genes_info[match(mod_red[,2],all_genes_info$ref_gene_name),3])

all_significant_genes <- rbind(DGE_genes,red_module)
duplicated_genes <- all_significant_genes[which(duplicated(all_significant_genes$gene_name)),]
dup_gene_list <- duplicated_genes[,1]

DGE_genes_only <- subset(DGE_genes,!(DGE_genes$gene_name %in% dup_gene_list)); DGE_genes_only$Analysis <- "DGE"
red_module_genes_only <- subset(red_module,!(red_module$gene_name %in% dup_gene_list)); red_module_genes_only$Analysis <- "WGCNA"
duplicated_genes$Analysis <- "Both"

file_s1 <- rbind(duplicated_genes,DGE_genes_only,red_module_genes_only)
write.csv(file_s1,file=paste(DataDir,"/File_S1.csv",sep=""))

###### Fig. S4: MODULE MEMBERSHIP AND INTRAMODULAR CONNECTIVITY

ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)
rownames(Alldegrees1) <- all_genes_info[match(rownames(Alldegrees1),all_genes_info$gene_id),1]
head(Alldegrees1)

# Calsulate gene significance
range = datTraits; names(range) = "Range"
GS1=as.numeric(cor(range,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)

colorlevels=unique(moduleColors)

i=6
par(mfrow=c(1,1))
whichmodule=colorlevels[[i]];
restrict1 = (moduleColors==whichmodule);
Degrees_red <- Alldegrees1[restrict1,]
GeneSignificance_red <- GeneSignificance[restrict1]
verboseScatterplot(Degrees_red$kWithin,
                   GeneSignificance_red, col="grey",
                   main="Module correlated with range",
                   xlab = "Connectivity", ylab = "Gene Significance", abline = FALSE)
top_genes_of_interest <- paste("LOC",genes_of_interest$LOC,sep="")
GOI <- cbind(Degrees_red[match(top_genes_of_interest, rownames(Degrees_red)),],"GeneSign"=GeneSignificance_red[match(top_genes_of_interest, rownames(Degrees_red))])
points(GOI$GeneSign~GOI$kWithin,pch=16,col="red")
with(GOI,text(GOI$GeneSign~GOI$kWithin, labels = row.names(GOI), pos = 1,cex=0.75))
dev.copy2pdf(file=paste(FigDir,"/WGCNA/signed_connectivity_significance_red_labelled.pdf",sep=""), width=8, height=6)