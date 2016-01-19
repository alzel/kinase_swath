#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "peptides_correlations"

## ---- data_load ----
load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations._load_.RData")

## ---- checking for correlations of peptides per protein ---- 
peptides.dataset = tbl_df(peptides.peak_sums.trimmed)
protein_annot.selected = droplevels(filter(protein_annotations, specificity == "unique"))

peptide2orfs = dplyr::select(protein_annot.selected, strippedSequence, SystName)
names(peptide2orfs) = c("EG.StrippedSequence", "ORF")

peptides.dataset.f = droplevels(peptides.dataset[peptides.dataset$EG.StrippedSequence %in% protein_annot.selected$strippedSequence,])

peak_sums = peptides.dataset.f
peak_sums.f = tbl_df(droplevels(merge(peak_sums, peptide2orfs, by="EG.StrippedSequence")))
peak_sums.f = peak_sums.f %>% distinct(R.Label, EG.StrippedSequence) %>% arrange(R.Label, ORF, EG.StrippedSequence)

peptides.wide = dcast(data=peak_sums.f, formula=ORF+EG.StrippedSequence~R.Label+batch.exp.n+batch, value.var="T_signal")
peptides.wide.selected = dcast(data=peak_sums.f, formula=ORF+EG.StrippedSequence+batch.exp.n~R.Label+batch, value.var="T_signal")

#peptides.long %>% group_by(ORF, variable, value, R.Label, batch.exp.n, batch_date, batch) %>% dplyr::summarise(sum(value))

batch.lengths = peak_sums.f %>% group_by(batch.exp.n) %>% summarize(counts = length(unique(R.Label))) %>% arrange(-counts)

#selecting top by size of batch
peptides.wide.selected = droplevels(peptides.wide.selected[peptides.wide.selected$batch.exp.n %in% batch.lengths$batch.exp.n[1:5],]) 

peptides.sel.cor = ddply(peptides.wide.selected, .(ORF, batch.exp.n), .parallel=T,
                        .fun = function(z) {
                          v <<- z
                          rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                          if (nrow(z) == 1) {
                            return(data.frame(id = z$EG.StrippedSequence, variable=row.names(z), value=NA))
                          }
                          
                          rownames(z) = paste(z$EG.StrippedSequence, sep=".")  
                          tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                          tmp.matrix[!upper.tri(tmp.matrix)] = NA
                          
                          tmp.df = data.frame(tmp.matrix)
                          tmp.df$id = rownames(tmp.df)
                          tmp.long = melt(tmp.df, id.vars=c("id"))
                          tmp.long = tmp.long[complete.cases(tmp.long),]
                          return(tmp.long)
                        })

peptides.sel.cor[sapply(peptides.sel.cor, is.character)] <- lapply(peptides.sel.cor[sapply(peptides.sel.cor, is.character)],  as.factor)

p = ggplot(peptides.sel.cor, aes(x=value)) +
           geom_density() +
           facet_wrap(~batch.exp.n)

plots.list = lappend(plots.list, p)


p = ggplot(peptides.sel.cor, aes(x=value)) +
          geom_histogram(colour = "black", fill = "white", binwidth = 0.05) + 
          #ggtitle(paste("Grouped by", paste(attr(toPlot,which="vars"), collapse=".")))
          xlab("Correlation between pairs of peptides per protein") +
          geom_vline(xintercept=0, linetype=2) +
          theme(axis.title=element_text(size=20),
                axis.text=element_text(size=14),
                aspect.ratio = 1)
        
plots.list = lappend(plots.list, p)

peptides.correlations = ddply(peptides.wide, .(ORF), .parallel=T,
                              .fun = function(z) {
                                v <<- z                      
                                rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                if (nrow(z) == 1) {
                                  return(data.frame(id = z$EG.StrippedSequence, variable=row.names(z), value=NA))
                                }
                                rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                tmp.matrix = cor(t(z[,-c(1:2)]), method="spearman", use="pairwise.complete.obs")
                                tmp.matrix[!upper.tri(tmp.matrix)] = NA
                                
                                tmp.df = data.frame(tmp.matrix)
                                tmp.df$id = rownames(tmp.df)
                                tmp.long = melt(tmp.df, id.vars=c("id"))
                                tmp.long = tmp.long[complete.cases(tmp.long),]
                                return(tmp.long)
                              })


file_name = paste("correlation_histogram", fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
png(file_path, width=210, height=297, units="mm", res=150)

par(pty="s")
hist(peptides.correlations$value, breaks=50, xlab="Spearman's correaltion, rho", main="", cex.axis=1.2, cex.lab=2)
legend("topleft", bg=NULL, bty="n", 
       legend=paste("Mean correlation ", round(mean(peptides.correlations$value, na.rm=T),2)))
p = recordPlot()
plots.list = lappend(plots.list, p)
dev.off()


peptides.correlations.subset = droplevels(peptides.correlations[(peptides.correlations$value > 0.5 | peptides.correlations$value < 0.5) & !is.na(peptides.correlations$value), ])
tmp1 = tbl_df(melt(peptides.wide[match(peptides.correlations.subset$id, peptides.wide$EG.StrippedSequence),], id.vars=c("ORF", "EG.StrippedSequence")))
tmp2 = tbl_df(melt(peptides.wide[match(peptides.correlations.subset$variable, peptides.wide$EG.StrippedSequence),], id.vars=c("ORF", "EG.StrippedSequence")))

peptides.correlations.long = tbl_df(cbind(tmp1,tmp2[,c("EG.StrippedSequence","value")]))
names(peptides.correlations.long)[2] = "EG.StrippedSequence1"
names(peptides.correlations.long)[5] = "EG.StrippedSequence2"
names(peptides.correlations.long)[4] = "value1"
names(peptides.correlations.long)[6] = "value2"

peptides.correlations.long$peptides = factor(paste(peptides.correlations.long$EG.StrippedSequence1, peptides.correlations.long$EG.StrippedSequence2, sep="."))
peptides.correlations.long = peptides.correlations.long %>% 
    group_by(ORF, EG.StrippedSequence1, EG.StrippedSequence2) %>% 
    arrange()

peptides.correlations.long = peptides.correlations.long %>% 
  extract(variable, into=c("R.Label", "batch.exp.n", "batch"), regex="(.*?)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")



# selecting top lowest correlated peptides
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(value))[1:20],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2)) + 
  geom_point() +
  ggtitle("top 20 low correlation peptides") + 
  facet_wrap(~peptides, scales="free") + theme(aspect.ratio=1) +
  geom_smooth(method="lm", fill=NA)
p

file_name = "top20_lowCorrelation.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, width=11.69+0.1*11.69, height=8.27+0.1*8.27)
plots.list = lappend(plots.list, p)


# selecting top highest correlated peptides
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(-value))[1:20],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2)) + 
  geom_point() +
  ggtitle("top 20 high correlation peptides") + 
  facet_wrap(~peptides, scales="free") +
  geom_smooth(method="lm", fill=NA)

file_name = "top20_highCorrelation.png"
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, width=11.69+0.1*11.69, height=8.27+0.1*8.27)
plots.list = lappend(plots.list, p)

# example of peptide correlations
#best 
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(-value))[1],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2)) + 
  geom_point() +
  ggtitle(toPlot$ORF[1]) + 
  xlab(paste("Peptide ", toPlot$EG.StrippedSequence1[1], ", log(intensity)", sep="")) +
  ylab(paste("Peptide ", toPlot$EG.StrippedSequence2[1], ", log(intensity)", sep="")) +
  geom_smooth(method="lm", fill=NA) +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=14),
        aspect.ratio = 1)

plots.list = lappend(plots.list, p)

#worst
toSelect = droplevels(unique(peptides.correlations.subset[with(peptides.correlations.subset, order(value))[1],c("id","variable")])) 
toSelect$peptides = with(toSelect, paste(id, variable, sep="."))

toPlot = droplevels(peptides.correlations.long[peptides.correlations.long$peptides %in% toSelect$peptides,])
p = ggplot(toPlot, aes(x=value1, y=value2)) + 
  geom_point() +
  ggtitle(toPlot$ORF[1]) + 
  xlab(paste("Peptide ", toPlot$EG.StrippedSequence1[1], ", log(intensity)", sep="")) +
  ylab(paste("Peptide ", toPlot$EG.StrippedSequence2[1], ", log(intensity)", sep="")) +
  geom_smooth(method="lm", fill=NA) +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=14),
        aspect.ratio = 1)

plots.list = lappend(plots.list, p)





#ranking peptides and additional stats

peak_sums.f[sapply(peak_sums.f, is.character)] <- lapply(peak_sums.f[sapply(peak_sums.f, is.character)],  as.factor)

peptides.ranks = peak_sums.f %>% group_by(R.Label, ORF) %>% 
                                    mutate(rev.rank = rank(-T_signal)) %>% # peptide rank per ORF
                                 group_by(R.Label) %>%
                                    mutate(chunk = cut(T_signal, 5, labels=c(1:5))) %>% # devide sample signals into chunks 
                                 group_by(R.Label, ORF) %>%
                                    mutate(peptide.count = n()) # how many peptides per protein


peptides.ranks$peptide.count.f = factor(with(peptides.ranks, ifelse(peptide.count > 10, "gt10", peptide.count)))


peptides.ranks.summary =  peptides.ranks %>% group_by(ORF, EG.StrippedSequence) %>%  
                                                summarize(mean.rank = mean(rev.rank),
                                                          CV.rank = sd(rev.rank)/mean(rev.rank),
                                                          mean.chunk= mean(chunk),
                                                          peptide.count = unique(peptide.count),
                                                          peptide.count.f = unique(peptide.count.f)) 
peptides.ranks.summary$peptide.count.f = factor(peptides.ranks.summary$peptide.count.f, 
  levels = with(peptides.ranks.summary, c(sort(as.numeric(levels(peptide.count.f)[-11])), levels(peptide.count.f)[11])))


# variation of peptide ranks among samples
p = ggplot(peptides.ranks.summary, aes(x=peptide.count.f, y=CV.rank)) +
           stat_boxplot(data = peptides.ranks.summary, aes(x=peptide.count.f, y=CV.rank), geom ='errorbar')+
           geom_boxplot() +
           ylab("Coefficient of variation of peptide ranks across samples") +
           xlab("Number of peptides per protein") +
           theme(axis.title=element_text(size=20),
                 axis.text=element_text(size=14),
                 aspect.ratio = 1)

plots.list = lappend(plots.list, p)


#checking if there is a strong effect if batches per peptide pair 
#apparently there is a strong effect
batch.cor.var = peptides.sel.cor %>% group_by(id,variable) %>% summarize(CV = sd(abs(value), na.rm=T)/mean(abs(value), na.rm=T))
#hist(batch.cor.var$CV, xlim=c(0,2))

str(peptides.wide)

peptides.cor.stats = ddply(peptides.wide, .(ORF), 
                            .fun = function(z) {
                              v <<- z                      
                              z = v
                              rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                              if (nrow(z) == 1) {
                                return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                              }
                              rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                              tmp.matrix = cor(t(z[,-c(1:2)]), method="spearman", use="pairwise.complete.obs")
                              
                              cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                              cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                              cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                              
                              ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                              return(ret)
                            })

p = ggplot(melt(peptides.cor.stats, id.vars = c("ORF", "EG.StrippedSequence")), aes(x=value, color=variable)) +
           geom_density() +
           xlab("Expected (average,median,shorth) Spearman's correlation coefficient between peptides") +
           theme(aspect.ratio = 1)

plots.list = lappend(plots.list, p)

peptides.ranks.summary = merge(peptides.ranks.summary, peptides.cor.stats[,c("ORF", "EG.StrippedSequence", "cor.means")], by = c("ORF", "EG.StrippedSequence"))
peptides.ranks.summary = peptides.ranks.summary %>% group_by(ORF) %>% mutate(topN = rank(mean.rank))


## ---- correlations of peptides considering their ranks ---- 
peptides.wide.topN = merge(peptides.wide, peptides.ranks.summary[,c("ORF", "EG.StrippedSequence", "topN")], by=c("ORF", "EG.StrippedSequence"))
peptides.wide.topN = peptides.wide.topN[,c(colnames(peptides.wide.topN)[ncol(peptides.wide.topN)], colnames(peptides.wide.topN)[-ncol(peptides.wide.topN)])]

peptides.cor.stats.top5 = ddply(filter(peptides.wide.topN, topN <=5), .(ORF), 
                           .fun = function(z) {
                             v <<- z                      
                             z = v
                             rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                             if (nrow(z) == 1) {
                               return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                             }
                             rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                             tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                             
                             cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                             cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                             cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                             
                             ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                             return(ret)
                           })


peptides.cor.stats.gt_top5 = ddply(filter(peptides.wide.topN, topN > 5), .(ORF), 
                                .fun = function(z) {
                                  v <<- z                      
                                  z = v
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  if (nrow(z) == 1) {
                                    return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                                  }
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                                  
                                  cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                                  cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                                  cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                                  
                                  ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                                  return(ret)
                                })

peptides.cor.stats.top3 = ddply(filter(peptides.wide.topN, topN <=3), .(ORF), 
                                .fun = function(z) {
                                  v <<- z                      
                                  z = v
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  if (nrow(z) == 1) {
                                    return(data.frame(EG.StrippedSequence = z$EG.StrippedSequence, cor.shorths=NA, cor.medians=NA, cor.means=NA))
                                  }
                                  rownames(z) = paste(z$EG.StrippedSequence, sep=".")
                                  tmp.matrix = cor(t(z[,-c(1:3)]), method="spearman", use="pairwise.complete.obs")
                                  
                                  cor.shorths = apply(tmp.matrix,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                                  cor.medians = apply(tmp.matrix,1, FUN=function(x) {median(x[x != 1 & !is.na(x)], na.rm=T)})
                                  cor.means   = apply(tmp.matrix,1, FUN=function(x) {mean(x[x != 1 & !is.na(x)], na.rm=T)})
                                  
                                  ret = data.frame(EG.StrippedSequence=names(cor.medians), cor.shorths, cor.medians, cor.means)
                                  return(ret)
                                })


peptides.ranks.summary = merge(peptides.ranks.summary, peptides.cor.stats[,c("ORF", "EG.StrippedSequence", "cor.means")], by = c("ORF", "EG.StrippedSequence"))

peptides.cor.stats.top3$top = 3
peptides.cor.stats.top5$top = 5
peptides.cor.stats.gt_top5$top = "gt_5"
tmp = peptides.cor.stats
tmp$top = "all"

peptides.cor.stats.top = rbind(peptides.cor.stats.top3, peptides.cor.stats.top5, peptides.cor.stats.gt_top5, tmp)
peptides.cor.stats.top.long = melt(peptides.cor.stats.top, id.vars=c("ORF", "EG.StrippedSequence", "top"))

#correlations tops density plots
toPlot = filter(peptides.cor.stats.top.long, variable=="cor.means")
toPlot.stats = toPlot %>% group_by(top) %>% summarize(mean.cor = mean(value, na.rm=T),
                                                      median.cor = median(value, na.rm=T),
                                                      shorth.cor = genefilter::shorth(value, na.rm=T)) 

p = ggplot(toPlot, aes(x=value, fill=top)) +
           geom_density(alpha=0.5) +
           geom_vline(data=toPlot.stats, aes(xintercept=shorth.cor, color=top), linetype=2, size=1) +
           xlab("Spearman's correlation, rho") +
           theme(axis.title=element_text(size=20),
                  axis.text=element_text(size=14),
                  aspect.ratio = 1)

plots.list = lappend(plots.list, p)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 

file_name = "peptides.cor.stats.top.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peptides.cor.stats.top,file=file_path)  

