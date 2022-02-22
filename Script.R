#DESeq2 Necesita datos sin normalizar. Después de descargar los archivos de TCGA usando el pipeline de Leukemias se puede usar la función DESeqDataSetFromMatrix para 
crear un objeto dds

dds <- DESeqDataSetFromMatrix(countData = round(rnas1[,factors$Group=="MM_III"|factors$Group=="NormalBM"]),
                              colData = factors[factors$Group=="MM_III"|factors$Group=="NormalBM",],
                              design = ~ Group)
dds <- DESeq(dds)
#Se aclara a DESeq quién es el grupo control o de referencia
dds$Group <- relevel(dds$Group, ref = "NormalBM")
dds <- DESeq(dds)
res <-  results(dds)
res
summary(res)
res05 <- results(dds, alpha=0.05)
resLFC <- lfcShrink(dds, coef="Group_MM_III_vs_NormalBM", type="apeglm")
write.table(resLFC, file = "resLFC_Group_MM_III_vs_NormalBM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]q

resLFC[abs(resLFC$log2FoldChange) >= 2 & resLFC$padj <= 0.05, ]

#### Volcano plot ####

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
resLFC <- read.table("resLFC_Group_MM_III_vs_MM_II.tsv", header=TRUE)

rownames(resLFC) <- gsub("-", ".", rownames(resLFC))

volcano <- EnhancedVolcano(resLFC,
                           lab= rownames(resLFC),
                           x = 'log2FoldChange',
                           y = 'padj',
                           pCutoff = 10^-5,
                           FCcutoff = 1.0,
                           pointSize = 3.0,
                           labSize = 6.0,
                           labCol = 'black',
                           labFace = 'bold',
                           boxedLabels = TRUE,
                           parseLabels = TRUE,
                           col = c('black', 'pink', 'purple', 'red3'),
                           colAlpha = 4/5,
                           legendPosition = 'bottom',
                           legendLabSize = 14,
                           legendIconSize = 4.0,
                           drawConnectors = TRUE,
                           widthConnectors = 1.0,
                           colConnectors = 'black') + coord_flip()

