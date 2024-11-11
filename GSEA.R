cat("\014"); rm(list = ls()); options(warn = -1); options(digits=3) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(org.Hs.eg.db);library(clusterProfiler);library(enrichplot);
library(ggplot2);library(DESeq2);library(openxlsx);library(biomaRt);library(progress)
library(enrichplot);library(latex2exp);library(ggVolcano);library(ggrepel);library(GseaVis)
# library(org.Mm.eg.db)
pvalueFilter=0.05
qvalueFilter=0.05

#--------------------------DESeq2 differential expression analysis------------
data = read.csv("TCGA-count120.csv",row.names = 1)
# load("TCGA-UCS_mrna_expr_counts.rdata");data = mrna_expr_counts
exprSet <- data

group <- read.csv("group.csv",col.names = T);colnames(group) = "label"
group <- as.factor(group$label);str(group)

# Create database
Data <- data.frame(row.names = colnames(exprSet), 
                   group = group)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)  # Standardized data
tmp <- results(dds2,contrast = c("group","pCR","non_pCR")) 
DEG_DESeq2 <- as.data.frame(tmp);DEG_DESeq2 <- na.omit(DEG_DESeq2)

# DESeq2 screens differential genes
FC <- 2;Padj <- 0.05
DEG_DESeq2$Significant <- "normal"
up <- intersect(which(DEG_DESeq2$log2FoldChange > log2(FC) ),
                which(DEG_DESeq2$padj < Padj))

down <- intersect(which(DEG_DESeq2$log2FoldChange < (-log2(FC))),
                  which(DEG_DESeq2$padj < Padj))
DEG_DESeq2$Significant[up] <- "up"
DEG_DESeq2$Significant[down] <- "down"
table(DEG_DESeq2$Significant)

#----------------------- Volcano mapping ------------------------------- 

DEG_DESeq2$color <- ifelse(-log10(DEG_DESeq2$padj) < 2, "#bfc0c1",
                        ifelse(DEG_DESeq2$log2FoldChange>0, "#e88182", "#6489b2") )

ggplot(DEG_DESeq2)+
  geom_point(aes(log2FoldChange, -log10(padj),
                 shape = Significant, size = -log10(padj)),
             color = DEG_DESeq2$color,
             alpha = 0.7
  )+
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_hline(yintercept = 2, linetype = "longdash")+
  scale_size_continuous(range = c(0.2, 3))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = TeX("$Log_2 \\textit{FC}$"),
       y = TeX("$-Log_{10} \\textit{FDR} $"))

ggsave("volcano plot.pdf", height = 5, width = 6)

# Volcano map add tags
label = rownames(DEG_DESeq2)
DEG_DESeq3 = cbind(DEG_DESeq2,label)
DEG_DESeq3$label <- rep("", nrow(DEG_DESeq3))
DEG_DESeq3$label[order(DEG_DESeq3$padj)[1:20]] <- rownames(DEG_DESeq3)[order(DEG_DESeq3$padj)[1:20]]

ggplot(DEG_DESeq3, aes(log2FoldChange, -log10(padj)))+
  geom_point(aes(shape = Significant, size = -log10(padj)), 
             color = DEG_DESeq3$color,
             alpha = 0.7)+
  geom_vline(xintercept = 0, linetype = "longdash") + 
  geom_hline(yintercept = 2, linetype = "longdash")+
  geom_text_repel(aes(label = label), size = 2, color = DEG_DESeq3$color,
                  max.overlaps = 100, key_glyph = draw_key_point)+
  scale_size_continuous(range = c(0.2, 3))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = TeX("$Log_2 \\textit{FC}$"), 
       y = TeX("$-Log_{10} \\textit{FDR} $"))

ggsave("volcano plot2.pdf", height = 4.2, width = 4.2)

# Output differential gene
up_gene = DEG_DESeq2[which(DEG_DESeq2$Significant == "up"),]
down_gene = DEG_DESeq2[which(DEG_DESeq2$Significant == "down"),]

write.csv(rownames(up_gene),"up_gene.csv",row.names = F,col.names = F)
write.csv(rownames(down_gene),"down_gene.csv",row.names = F,col.names = F)

#----------------------- GSEA analysis :GO ontology enrichment----------------------
data <- DEG_DESeq2;dim(data)
data_sort <- data %>%   #  Make the gene_list set by sorting logFC columns in descending order
  arrange(desc(log2FoldChange))
gene_list <- data_sort$log2FoldChange
names(gene_list) <- rownames(data_sort)

res <- gseGO(
  gene_list, 
  ont = "BP",  
  OrgDb = org.Hs.eg.db, 
  keyType = "SYMBOL", 
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
)

# Draw bubble diagram
pdf("GO_1_dotplot.pdf", wi=9, h=7)
dotplot(res)
dotplot(
  res, 
  showCategory=10, 
  split=".sign") + facet_grid(.~.sign)
dev.off()

# Cycle to make a typical enrichment map
dir.create("GSEAPlot");setwd(paste0(getwd(),"/GSEAPlot"))
N = length(res@result[["ID"]]);pb <- progress_bar$new(total=N) 

for (i in 1:length(res@result[["ID"]])) {
  pb$tick()
  pdf_file <- paste0(res$Description[i], "-GSEA.pdf")
  pdf(pdf_file, wi=9, h=7)
  a <- gseaplot2(res, title = res$Description[i], geneSetID = i);print(a) 
  dev.off()
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pairwise_termsim(res)

# Cycle to make a typical enrichment map
pdf("GSEA.pdf", wi=9, h=7)
gseaplot2(res, title = res$Description[1], geneSetID = 1);
dev.off()

# Network diagram
pdf("GO_3_CNEplot.pdf", wi=9, h=7)
cnetplot(res, categorySize="pvalue", foldChange=gene_list)
dev.off()

# Ridge map drawing
pdf("GO_4_ridgeplot.pdf", wi=9, h=7)
ridgeplot(res) + labs(x = "enrichment distribution")
dev.off()

# NES score as X axis
pdf("GO_6_GSEAdotplot2.pdf", wi=9, h=7)
dotplotGsea(data = res,topn = 10,
            order.by = 'NES')
dev.off()

# free scale Adds a grid
pdf("GO_8_lollipopplot2.pdf", wi=16, h=7)
dotplotGsea(data = res,topn = 10,
            order.by = 'NES',
            add.seg = T,
            scales = 'free')
dev.off()

# adjust the pajust and pval thresholds, and choose the term that suits you to show
p1 <- dotplotGsea(data = res,topn = 5,
                  pajust = 0.1)$plot

p2 <- dotplotGsea(data = res,topn = 5,
                  pval = 0.05)$plot

# combine
pdf("GO_9_lollipopplot.pdf", wi=9, h=8)
cowplot::plot_grid(p1,p2,nrow = 2,align = 'hv')
dev.off()

showNum = 5

#-----------------------GSEA analysis :KEGG pathway enrichment-----------------------
# Convert symbol to ENTREZID
DEG_DESeq2$ENTREZID <- mapIds(org.Hs.eg.db, keys = rownames(data), column = "ENTREZID", keytype = "SYMBOL")

data3 <- DEG_DESeq2;dim(data3)
data_sort3 <- data3 %>%  
  arrange(desc(log2FoldChange))
gene_list <- data_sort3$log2FoldChange
names(gene_list) <- data_sort3$ENTREZID

res <- gseKEGG(
  gene_list, 
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# 绘制气泡图
pdf("KEGG_1_dotplot.pdf", wi=9, h=7)
dotplot(res) 
dotplot(
  res, 
  showCategory=10, 
  split=".sign") + facet_grid(.~.sign)
dev.off()


pdf("KEGG_2_gseaplot.pdf", wi=9, h=7)
gseaplot2(res, title = res$Description[1], geneSetID = 1) 
dev.off()

# gseaNb(object = res,
#        geneSetID = 'GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS')

pairwise_termsim(res)

# Network diagram
pdf("KEGG_3_CNEplot.pdf", wi=9, h=7)
cnetplot(res, categorySize="pvalue", foldChange=gene_list)
dev.off()

# Ridge map drawing
pdf("KEGG_4_ridgeplot.pdf", wi=9, h=7)
ridgeplot(res) + labs(x = "enrichment distribution")
dev.off()

# Draw a lollipop chart
library(clusterProfiler);library(GseaVis)
pdf("KEGG_5_GSEAdotplot1.pdf", wi=9, h=7)
dotplotGsea(data = res,topn = 10)
dev.off()

# NES score as X axis
pdf("KEGG_6_GSEAdotplot2.pdf", wi=9, h=7)
dotplotGsea(data = res,topn = 10,
            order.by = 'NES')
dev.off()

# Add lines and draw lollipops
pdf("KEGG_7_lollipopplot.pdf", wi=9, h=7)
dotplotGsea(data = res,topn = 10,
            order.by = 'NES',
            add.seg = T)
dev.off()

# free scale Adds a grid
pdf("KEGG_8_lollipopplot2.pdf", wi=16, h=7)
dotplotGsea(data = res,topn = 10,
            order.by = 'NES',
            add.seg = T,
            scales = 'free')
dev.off()

# adjust the pajust and pval thresholds, and choose the term that suits you to show
p1 <- dotplotGsea(data = res,topn = 5,
                  pajust = 0.1)$plot

p2 <- dotplotGsea(data = res,topn = 5,
                  pval = 0.05)$plot

# combine
pdf("KEGG_9_lollipopplot.pdf", wi=9, h=8)
cowplot::plot_grid(p1,p2,nrow = 2,align = 'hv')
dev.off()
