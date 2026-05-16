####Fig2 A-G

### Fig2 A-B

library(rstatix)
library(ggh4x)
library(ggpubr)

combat_data = read.csv("SCFA_count.csv",row.names = 1)
group = read.csv("SCFA_metadata.csv")
otu = t(combat_data) %>% as.data.frame() %>%  tibble::rownames_to_column("SampleID")
dat = dplyr::inner_join(group,otu,by = c("Sample_ID"="SampleID")) %>% mutate(Mouth=paste0("Mouth_",gsub(".*_","",Group)))
dat$Mouth <- factor(dat$Mouth,levels = c(paste0("Mouth_",c(1,3,6,9))))


lower_quantile <- quantile(dat$AA, 0.25)
upper_quantile <- quantile(dat$AA, 0.75)
dat$Group2 <- factor(dat$Group2, levels = c("WT", "RTT"))
p = dat %>% ggplot(aes(Group2,AA))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = lower_quantile, 
           fill = "cornflowerblue", alpha = .3, color = NA)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = upper_quantile, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=Mouth),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  geom_hline(yintercept = upper_quantile,linetype=2)+
  geom_hline(yintercept = lower_quantile,linetype=2)+
  facet_nested_wrap(. ~ Mouth,ncol = 4,
                    strip = strip_nested(background_x =
                    elem_list_rect(fill=c("#3B9AB2","#7294D4","#9E8CC6","#E6A0C4")))) +
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#9E8CC6","#E6A0C4"))+
  labs(x=NULL)+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"))+ 
  stat_compare_means(method = "wilcox.test",comparisons = list(c("RTT","WT")),
                     label.x = 1.3, label.y = max(dat$AA)+0.02,label = "p.signif")
ggsave("fig2_A.pdf",p,width = 6,height=4)


lower_quantile <- quantile(dat$PA, 0.25)
upper_quantile <- quantile(dat$PA, 0.75)
dat$Group2 <- factor(dat$Group2, levels = c("WT", "RTT"))
p = dat %>% ggplot(aes(Group2,PA))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = lower_quantile, 
           fill = "cornflowerblue", alpha = .3, color = NA)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = upper_quantile, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=Mouth),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  geom_hline(yintercept = upper_quantile,linetype=2)+
  geom_hline(yintercept = lower_quantile,linetype=2)+
  facet_nested_wrap(. ~ Mouth,ncol = 4,
                    strip = strip_nested(background_x =
                                           elem_list_rect(fill=c("#3B9AB2","#7294D4","#9E8CC6","#E6A0C4")))) +
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#9E8CC6","#E6A0C4"))+
  labs(x=NULL)+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"))+ 
  stat_compare_means(method = "wilcox.test",comparisons = list(c("RTT","WT")),
                     label.x = 1.3, label.y = max(dat$PA)+0.02,label = "p.signif")
ggsave("fig2_B.pdf",p,width = 6,height=4)



### Fig2 C-E

library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(biomaRt)
library(pheatmap)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(ggvenn) 
library(WGCNA)

# WGCNA -------------------------------------------------------------------
metadata = read.csv("SCFA_metadata.csv",header = T,row.names = 2)
taxonomy = read.table("taxonomy.txt",header = T,row.names = 1)
otu = read.table("otutab.txt",header = T,row.names = 1)
SCFA = read.csv("SCFA_count.csv",header = T,row.names = 1)

temp = taxonomy[taxonomy$Genus%in%"Unassigned",]
otu= otu[-which(row.names(otu)%in%row.names(temp)),]
otu= log10(otu+1)
range(otu)
range(SCFA)
SCFA = SCFA[,colnames(otu)]
identical(colnames(otu),colnames(SCFA))
otu = rbind(otu,SCFA)
datExpr0 <- t(scale(otu))
datExpr0 <- data.frame(datExpr0)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

datTraits = metadata
datTraits <- datTraits[match(rownames(datExpr0),rownames(datTraits)),]
identical(rownames(datTraits),rownames(datExpr0))

datTraits1 <- datTraits[,2] 
datTraits1 <- as.data.frame(datTraits1)
rownames(datTraits1) <- rownames(datTraits)

pheno <- binarizeCategoricalColumns(datTraits1, dropFirstLevelVsAll = F,minCount=0)
colnames(pheno) <- c('RTT1','RTT3','RTT6','RTT9','WT1','WT3','WT6','WT9')
row.names(pheno) <- rownames(datTraits1)

cor = WGCNA::cor 
mergingThresh=0.15;softPower=16;
oneStep_net = blockwiseModules(datExpr0,
                               corType="pearson",
                               networkType = "signed",TOMType = "signed",
                               power=softPower,
                               minModuleSize=20,
                               mergeCutHeight=mergingThresh,
                               numericLabels=TRUE,
                               saveTOMs=F,pamRespectsDendro=FALSE,
                               nThreads=40)
table(oneStep_net$colors)
cor = stats::cor
oneStep_modulecolor <- labels2colors(oneStep_net$colors)

plotDendroAndColors(oneStep_net$dendrograms[[1]],
                    oneStep_modulecolor[oneStep_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

geneCor=matrix(NA,nrow=8,ncol=ncol(datExpr0)) 
WT1=as.numeric(pheno$WT1)
WT3=as.numeric(pheno$WT3)
WT6=as.numeric(pheno$WT6)
WT9=as.numeric(pheno$WT9)
RTT1=as.numeric(pheno$RTT1)
RTT3=as.numeric(pheno$RTT3)
RTT6=as.numeric(pheno$RTT6)
RTT9=as.numeric(pheno$RTT9)
for(i in 1:ncol(geneCor)) {
  expr=as.numeric(datExpr0[,i])
  WT1_r=bicor(expr,WT1,use="pairwise.complete.obs")
  WT3_r=bicor(expr,WT3,use="pairwise.complete.obs")
  WT6_r=bicor(expr,WT6,use="pairwise.complete.obs")
  WT9_r=bicor(expr,WT9,use="pairwise.complete.obs")
  RTT1_r=bicor(expr, RTT1,use="pairwise.complete.obs")
  RTT3_r=bicor(expr, RTT3,use="pairwise.complete.obs")
  RTT6_r=bicor(expr, RTT6,use="pairwise.complete.obs")
  RTT9_r=bicor(expr, RTT9,use="pairwise.complete.obs")
  geneCor[,i]=c(WT1_r, WT3_r, WT6_r, WT9_r, RTT1_r,RTT3_r,RTT6_r,RTT9_r)
  cat('Done for gene...',i,'\n')
}
for(i in 1:nrow(geneCor)){
  geneCor[i,] =numbers2colors(as.numeric(geneCor[i,]),signed=TRUE,centered=TRUE,
                              colorRampPalette(c("#216FB1","white","#BE2834"))(50),
                              lim=c(-0.5,0.5))
}
rownames(geneCor)=unique(datTraits$Group)

pdf("fig2_C.pdf", width = 10, height = 6)
plotDendroAndColors(oneStep_net$dendrograms[[1]], 
                      cbind(oneStep_modulecolor,geneCor[1,], geneCor[2,],
                            geneCor[3,],geneCor[4,],geneCor[5,],geneCor[6,], geneCor[7,],geneCor[8,]),
                      groupLabels=c("Modules",'WT1','WT3','WT6','WT9','RTT1','RTT3','RTT6','RTT9'),
                      dendroLabels = FALSE, hang = 0.01,
                      addGuide = TRUE, guideHang = 0.01)
dev.off()


library(linkET)
MEs0 = moduleEigengenes(datExpr0, oneStep_modulecolor)$eigengenes
ME_merge <- orderMEs(MEs0)

moduleTraitCor <- cor(
  ME_merge,
  pheno,
  use = "pairwise.complete.obs",
  method = "spearman"
)

moduleTraitPvalue <- matrix(
  NA,
  nrow = ncol(ME_merge),
  ncol = ncol(pheno)
)

rownames(moduleTraitPvalue) <- colnames(ME_merge)
colnames(moduleTraitPvalue) <- colnames(pheno)

for (i in seq_len(ncol(ME_merge))) {
  for (j in seq_len(ncol(pheno))) {
    ok <- complete.cases(ME_merge[, i], pheno[, j])
    moduleTraitPvalue[i, j] <- cor.test(
      ME_merge[ok, i],
      pheno[ok, j],
      method = "spearman",
      exact = FALSE
    )$p.value
  }
}


ME_gene = rep(NA,length(names(ME_merge)))
for(i in 1:length(names(ME_merge))){
  modules <- names(ME_merge)[i]
  modules <- gsub("ME", "", modules)
  sites <- which(names(table(oneStep_modulecolor))==modules)
  ME_gene[i] = paste0(names(ME_merge)[i],"(",table(oneStep_modulecolor)[[sites]],")",sep="")
}

rownames(moduleTraitCor) <- ME_gene
rownames(moduleTraitPvalue) <- ME_gene

mantel <- linkET::mantel_test(pheno, ME_merge,
          spec_select = list(RTT1 = 1,RTT3 = 2,RTT6 = 3,RTT9 = 4,WT1 = 5,WT3 = 6,WT6 = 7,WT9 = 8)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, 0.8, Inf),
                  labels = c('[0, 0.4]', '(0.4, 0.6]', '(0.6, 0.8]', '(0.8, 1]')),
         pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                  labels = c('***', '**', '*', 'ns')),
         linetype = cut(r, breaks = c(-Inf,0,Inf),labels = c('Neqative Relation','Positive Relation')))
mantel = as.data.frame(mantel)
mantel[,3]=c(moduleTraitCor[,1],moduleTraitCor[,2],moduleTraitCor[,3],moduleTraitCor[,4],moduleTraitCor[,5],moduleTraitCor[,6],moduleTraitCor[,7],moduleTraitCor[,8])
mantel[,5]=cut(mantel[,3], breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),labels = c('[0, 0.2]', '(0.2, 0.4]', '(0.4, 0.6]', '(0.6, 0.8]', '(0.8, 1]'))

mantel[,4]=c(moduleTraitPvalue[,1],moduleTraitPvalue[,2],moduleTraitPvalue[,3],moduleTraitPvalue[,4],moduleTraitPvalue[,5],moduleTraitPvalue[,6],moduleTraitPvalue[,7],moduleTraitPvalue[,8])
mantel[,6]=cut(mantel[,4], breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),labels = c('***', '**', '*', 'ns'))
mantel[,7]=cut(mantel[,3], breaks = c(-Inf,0,Inf),labels = c('Neqative Relation','Positive Relation'))

p=linkET::qcorrplot(correlate(ME_merge, method = "spearman"), type = "lower", diag = FALSE)+
  set_corrplot_style(colours =colorRampPalette(c("#99CC00","#EEEEEE","#FF597B"))(100) )+##颜色
  geom_square() +
  geom_couple(aes(colour = pd, size = rd,linetype = linetype), data = mantel, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.2, 0.4,0.6,0.8)) +
  scale_colour_manual(values = color_pal(4)) +
  scale_linetype_manual(values = c("dotted","solid")) +
  guides(size = guide_legend(title = "Spearman's rho",override.aes = list(colour = "#000000"), order = 2),
         colour = guide_legend(title = "Spearman's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's rho", order = 3))
ggsave("fig2_D.pdf",p,width = 6,height=4)


library(WeightedTreemaps)
library(tibble)
module_dataframe = read.csv("module_dataframe_gene.csv")
module_dataframe = subset(module_dataframe,module_color%in%"blue")
ASVs <- taxonomy[module_dataframe$gene_id,] %>% rownames_to_column("OTUID")
ASVs$Abundance = log2(apply(otu[module_dataframe$gene_id,],1,sum))
ASVs = ASVs[!is.na(ASVs$Kingdom), ]

data <- voronoiTreemap(
  data = ASVs,
  levels = c("Phylum", "Genus"),
  filter = 0.001,
  cell_size = "Abundance",
  shape = "circle",
  positioning = "clustered_by_area",
  seed = 10)

pdf("fig2_E.pdf",width = 6,height = 5.5)
drawTreemap(data,
            color_level=1,
            color_type ="both",
            color_palette =c("#7E78D6", "#64A6B1", "#00506B","#33E6A1"),
            label_size = 2,
            label_color ="white",
            border_color="white",
            border_size = 6,
            title ="Blue module",
            title_size = 2,
            title_color ="black",
            legend = T,
            legend_position ="right",
            legend_size = 0.2,
            width =0.9,height = 0.8,
            layout =c(1,1),
            position =c(1,1),
            add =T)
dev.off()



### Fig2 F-G

library(Seurat)
metabolite.cpd2ko = read.table("hsa.metabolite.cpd2ko.240703.tsv",sep = "\t",header = T)
metabolism = read.csv("metabolism_month9.csv")
metabolism_pathways = read.csv("pathway_results.csv",row.names = 1) %>% tibble::rownames_to_column("Pathway") %>% subset(Impact>0.1) %>% arrange(Impact)
metabolism_map = read.csv("pathway_conpound.csv")

metabolisms = metabolism %>% dplyr::filter(RTT.9_vs_CON.9_VIP>1 & RTT.9_vs_CON.9_P.value<0.05 & abs(RTT.9_vs_CON.9_Log2FC)>1)
metabolisms = metabolisms %>%mutate(type = case_when(RTT.9_vs_CON.9_Log2FC > 0 ~ "up",RTT.9_vs_CON.9_Log2FC < 0 ~ "down"))

mycolor <- c("#96C3D8", "#F5B375","#67A59B", "#A5D38F", "#8D75AF", "#F19294", "#E45D61", "#BDA7CB","#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")
colors2 <- c(mycolor,colorRampPalette((pal_npg("nrc")(9)))(15))
res_inter <- as.data.frame(table(metabolisms$Class.I))
colnames(res_inter) = c("super class","value")
p <- ggpubr::ggbarplot(res_inter, x="super class", y="value", fill = "super class", color = "super class", lab.size = 3,
                  sort.by.groups=FALSE, sort.val = "desc", palette = colors2,
                  label = T, xlab = "Class I", ylab = "number") + ggpubr::rotate_x_text(90)+NoLegend()
ggsave("fig2_G.pdf", p, width = 10, height = 5)


res_heatmap = metabolisms[, 17:34]
row.names(res_heatmap) = metabolisms$Compounds
group = data.frame(
  row.names = colnames(res_heatmap),
  group = c(rep("RTT", 9), rep("WT", 9))
)
res_heatmap = res_heatmap[, c(10:18, 1:9)]
group = group[colnames(res_heatmap), , drop = FALSE]
row_anno_colors = list(
  group = c(RTT = '#FF597B', WT = '#D0E6A5')
)

pdf("fig2_F.pdf", width = 6, height = 8)
pheatmap::pheatmap(
  res_heatmap,
  annotation_col = group,
  color = colorRampPalette(c("#99CC00", "#EEEEEE", "#FF597B"))(50),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_colors = row_anno_colors,
  border_color = "#EEEEEE",
  scale = "row",
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 8
)
dev.off()


