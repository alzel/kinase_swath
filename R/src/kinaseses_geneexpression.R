data.raw <- read.delim("~/projects/kinase_swath/data/2015-09-27/DataS1.txt", header=F)
load("./R/objects/protein_annotations._load_.RData")


toExclude = grep(x=as.character(t(data.raw[2,])), pattern="p_value", perl=T)
data.raw = data.raw[-2,-toExclude]

colNames = gsub(pattern="-del", replacement="", t(data.raw[1,]))
colNames = gsub(pattern=" |vs.|wt", replacement="", colNames)
colNames = gsub(pattern="[+-]", replacement="_", colNames)

data.raw = data.raw[-1,]
colnames(data.raw) = colNames
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

data = data.raw
data[, -c(1:2)] <- sapply(data[, -c(1:2)], as.numeric.factor)

rownames(data) = data$geneSymbol

orf2name = unique(data.frame(ORF = protein_annotations$SystName, 
                             gene_name = protein_annotations$sgdName))
data.selected = data[grep(pattern="HXk2|HXT", ignore.case=T, x=data$geneSymbol),]
rownames(data.selected) = data.selected$geneSymbol
pairs(t(data.selected[,-c(1,2)]))
cor(t(data.selected[,-c(1,2)]))
pairs(t(data.selected)

load("./R/objects/metabolitesTCA_metadata._clean_.RData")

abs(data.selected[,-c(1,2)]) > 0.15

metabolitesTCA_metadata$gene_name = orf2name$gene_name[match(metabolitesTCA_metadata$ORF, orf2name$ORF)]
data.selected.tca = data.selected[,c(1,2,which(toupper(colnames(data.selected)) %in% metabolitesTCA_metadata$gene_name))]
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)   {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- cor(x, y)
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * r)
 }

pairs(t(data.selected.tca[,-c(1,2)]), upper.panel=panel.cor)
cor(t(data.selected.tca[,-c(1,2)]))


example(pairs)
