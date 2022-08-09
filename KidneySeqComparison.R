library(ggplot2)
library(org.Mm.eg.db)
library(limma)
library(edgeR)

##### Kidney FGF21 LP Data (Fang et al. PMID 34151592)

fgfData <- read.csv("kidneyFgfCounts.csv")
geneVec <- fgfData[,1]
rownames(fgfData) <- fgfData[,1]
fgfData <- fgfData[-1]

fgfProtGroup <- factor(c(rep("Ctrl", 6), rep("LP", 6),
                      rep("Ctrl", 6), rep("LP", 6)))

fgfGtGroup <- factor(c(rep("KO", 12), rep("WT", 12)),
                  levels = c("WT", "KO"))

fgfProtGtGroup <- interaction(fgfGtGroup, fgfProtGroup)

xF <- DGEList(counts = fgfData, genes = geneVec)
#x <- DGEList(counts = rawData, genes = geneVec)
xF$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(xF),
                         keytype = "ENSEMBL", column = "SYMBOL")
xF$genes$entrez <- mapIds(org.Mm.eg.db, rownames(xF),
                         keytype = "ENSEMBL", column = "ENTREZID")

fgfKeep <- rowSums( cpm(xF) >1) >= 1
table(fgfKeep)
xF <- xF[fgfKeep, , keep.lib.sizes=FALSE]
xF <- calcNormFactors(xF)

##### Kidney bulk seq

rawDataK <- read.csv("kidneyCounts.csv")
geneVecK <- rawDataK[,1]
rownames(rawDataK) <- rawDataK[,1]
rawDataK <- rawDataK[-1]

rawDataSelK <- rawDataK[, -grep("HI", colnames(rawDataK))]
rawDataSelK <- rawDataSelK[, -grep("CI", colnames(rawDataK))]

colnames(rawDataK)

txGroupK <- factor(c(rep("HP_LC", 5), rep("HP_IGF", 3), rep("HP_HC", 6),
                    rep("HC_IGF", 5), rep("LP_HC", 5), rep("LP_LC", 6)),
                  levels = c("HP_LC", "HP_HC", "LP_LC", "LP_HC",
                             "HP_IGF", "HC_IGF"))

carbGroupK <- factor(c(rep("LC", 5), rep("HC", 11), rep("LC", 6)),
                    levels = c("LC", "HC"))
protGroupK <- factor(c(rep("HP", 11), rep("LP", 11)),
                    levels = c("HP", "LP"))
protGroupK2 <- factor(c(rep("High", 5), rep("Med", 6),
                       rep("Low", 5), rep("Med", 6)),
                     levels = c("High", "Med", "Low"))
dietGroupK <- interaction(carbGroupK, protGroupK)

xK <- DGEList(counts = rawDataSelK, genes = geneVecK)
#x <- DGEList(counts = rawData, genes = geneVec)
xK$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(xK),
                         keytype = "ENSEMBL", column = "SYMBOL")
xK$genes$entrez <- mapIds(org.Mm.eg.db, rownames(xK),
                         keytype = "ENSEMBL", column = "ENTREZID")

keepK <- rowSums( cpm(xK) >1) >= 1
table(keep)
xK <- xK[keepK, , keep.lib.sizes=FALSE]
xK <- calcNormFactors(xK)

##### Liver seq data (MacArthur et al Cell Reports)

rawDataL <- read.csv("prCounts.csv")
geneVecL <- rawDataL[,1]
rownames(rawDataL) <- rawDataL$X
rawDataL <- rawDataL[,-1]
rawDataL <- rawDataL[, grep("PER", colnames(rawDataL))]

xL <- DGEList(counts = rawDataL, genes = geneVecL)
xL$genes$symbol <- mapIds(org.Mm.eg.db, rownames(xL),
                         keytype = "ENSEMBL", column = "SYMBOL")
xL$genes$entrez <- mapIds(org.Mm.eg.db, rownames(xL),
                         keytype = "ENSEMBL", column = "ENTREZID")

txGroupL <- factor(c(rep("pct00", 4), rep("pct10", 4), rep("pct14", 4),
                    rep("pct18", 4), rep("pct02", 4), rep("pct06", 4)),
                  levels = c("pct18", "pct14", "pct10", "pct06", "pct02", "pct00"))

keepL <- rowSums( cpm(xL) >1) >= 2
table(keepL)
xL <- xL[keepL, , keep.lib.sizes=FALSE]
xL <- calcNormFactors(xL)

save(rawDataK, rawDataL, fgfData,
     xF, xK, xL, file = "compSeq.RData")


