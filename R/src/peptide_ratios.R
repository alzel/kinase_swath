#!/usr/bin/env Rscript
rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "peptide_ratios"

## ---- data_load ----
load("./R/objects/peptides.peak_sums.trimmed.RData")
load("./R/objects/protein_annotations.RData")
load("./R/objects/sample_exp.map.RData")
load("./R/objects/experiment.map.RData")

peptides.peak_sums.trimmed
peptides.dataset = tbl_df(peptides.peak_sums.trimmed)
protein_annot.selected = droplevels(filter(protein_annotations, specificity == "unique"))

peptide2orfs = dplyr::select(protein_annot.selected, strippedSequence, SystName)
names(peptide2orfs) = c("EG.StrippedSequence", "ORF")

peptides.dataset.f = droplevels(peptides.dataset[peptides.dataset$EG.StrippedSequence %in% protein_annot.selected$strippedSequence,])

peak_sums = peptides.dataset.f
peak_sums.f = tbl_df(droplevels(merge(peak_sums, peptide2orfs, by="EG.StrippedSequence")))

max(peak_sums.f$peptide.prevalence)
peak_sums.f = peak_sums.f %>% group_by(EG.StrippedSequence) %>% mutate(peptide.prevalence = n())
peak_sums.f = peak_sums.f %>% distinct(R.Label, EG.StrippedSequence) %>% arrange( ORF,R.Label, EG.StrippedSequence, -signal)

tmp.counts = peak_sums.f %>% group_by(ORF, R.Label) %>% summarize(pep.count = n())
unique(tmp.counts[which(tmp.counts$pep.count == 1)[],"ORF"])

stopifnot(max(table(tmp.counts$ORF)) == min(table(tmp.counts$ORF)))



peptides.ranks = peak_sums.f %>% group_by(R.Label, ORF) %>% 
  mutate(rev.rank = rank(-signal)) %>% # peptide rank per ORF
  group_by(R.Label) %>%
  mutate(chunk = cut(signal, 5, labels=c(1:5))) %>% # devide sample signals into chunks 
  group_by(R.Label, ORF) %>%
  mutate(peptide.count = n()) # how many peptides per protein

library(doSNOW)
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)

peak_sums.ratios =  ddply(peptides.ranks, .(ORF, R.Label ), .parallel=T,
      .fun = function(z) {
        v <<- z
        z = v
        z = droplevels(z)
        if (length(z$ORF) > 1) {
          #tmp.names  = t(combn(as.character(z$EG.StrippedSequence), 2))
          ind = 1:length(z$EG.StrippedSequence)
          tmp.ind  = t(combn(ind, 2))
          tmp.names = matrix(z$EG.StrippedSequence[as.vector(tmp.ind)], ncol=2, byrow=F) 
          tmp.signal = matrix(z$signal[as.vector(tmp.ind)], ncol=2, byrow=F) 
          tmp.rank = matrix(z$rev.rank[as.vector(tmp.ind)], ncol=2, byrow=F) 
                    
          #tmp.signal = t(combn(z$signal, 2))
          ratio = tmp.signal[,2]/tmp.signal[,1]
          mean.rank = rowMeans(tmp.rank)
          ret = data.frame(tmp.names, ratio, mean.rank, z$peptide.count[1])
          names(ret) = c("peptide1", "peptide2", "ratio", "mean.rank", "peptide.count")
          return(ret)
          
        } else {
          ret = data.frame(z$EG.StrippedSequence,z$EG.StrippedSequence, 1, 1, z$peptide.count[1])
          names(ret) = c("peptide1", "peptide2", "ratio", "mean.rank", "peptide.count")
          return(ret)
        }
      })
stopCluster(cl)




peak_sums.ratios = tbl_df(peak_sums.ratios)
peak_sums.ratios$pair = paste(peak_sums.ratios$peptide1, peak_sums.ratios$peptide2, sep=".")

peak_sums.ratios$KO = factor(sample_exp.map$ORF[match(peak_sums.ratios$R.Label, sample_exp.map$SampleName)])
peak_sums.ratios$KO[is.na(peak_sums.ratios$KO)] = "none"

file_name = "peak_sums.ratios.RData"
file_path = paste(output_dir, file_name, sep="/")
save(peak_sums.ratios, file=file_path)


load("./R/objects/peak_sums.ratios.RData")
peak_sums.ratios.f = tbl_df(peak_sums.ratios[grep(x=peak_sums.ratios$R.Label, pattern="mix", ignore.case=T, invert=T),])


#tmp = peak_sums.ratios.f %>% distinct(ORF, R.Label) %>% group_by(ORF) %>% summarise(insamples = n())

peak_sums.ratios.f.CV = peak_sums.ratios.f %>% group_by(ORF, pair) %>% dplyr::summarize( CV = sd(ratio, na.rm=T)/mean(ratio, na.rm=T),
                                                                                        rCV = mad(ratio, na.rm=T)/median(ratio, na.rm=T),
                                                                                        median.rank = median(mean.rank, na.rm=T))

cor(peak_sums.ratios.f.CV$rCV, peak_sums.ratios.f.CV$rank)


peak_sums.ratios.f.CV$pair = factor(peak_sums.ratios.f.CV$pair)
peak_sums.ratios.f.CV.long = melt(data.frame(peak_sums.ratios.f.CV), id.vars=c("ORF", "pair"))


p1 = ggplot(filter(peak_sums.ratios.f.CV.long, variable != "median.rank"), aes(x=value)) +
      geom_histogram(fill="white", color="black") +
      facet_wrap(~variable, scales="free") +
      theme(aspect.ratio = 1)

pep.count = 6
selectedORF = as.character(tmp.counts$ORF[which(tmp.counts$pep.count == pep.count)[1]])

p2 = ggplot(filter(peak_sums.ratios, ORF == selectedORF ), aes(x=pair, y=ratio)) +
      geom_boxplot() +
      ggtitle(paste(selectedORF,"with", pep.count, "peptides", sep=" ")) +
      theme(axis.text.x = element_text(angle = 15, hjust = 1))

p3 = ggplot(peak_sums.ratios.f.CV, aes(x=peak_sums.ratios.f.CV$median.rank, y=peak_sums.ratios.f.CV$rCV)) +
      geom_point() + 
      theme(aspect.ratio=1)

g = arrangeGrob(p2, arrangeGrob(p1,p3, ncol=2))
plots.list = lappend(plots.list, g)


# -- Peptide ratios from Picotti et al ----

load("./R/objects/peptides.LIT.RData")
peptides.LIT.long = melt(peptides.LIT, id.vars = c("ORF", "EG.StrippedSequence"))
peptides.LIT.long = peptides.LIT.long %>% arrange(ORF, variable, EG.StrippedSequence, -value)
                    
peptides.LIT.ratios =  ddply(peptides.LIT.long, .( ORF,  variable),
                          .fun = function(z) {
                            v <<- z
                            z = v
                            z = droplevels(z)
                            if (length(z$ORF) > 1) {
                              #tmp.names  = t(combn(as.character(z$EG.StrippedSequence), 2))
                              ind = 1:length(z$EG.StrippedSequence)
                              tmp.ind  = t(combn(ind, 2))
                              tmp.names = matrix(z$EG.StrippedSequence[as.vector(tmp.ind)], ncol=2, byrow=F) 
                              tmp.signal = matrix(z$value[as.vector(tmp.ind)], ncol=2, byrow=F) 
                            
                              #tmp.signal = t(combn(z$signal, 2))
                              ratio = tmp.signal[,2]/tmp.signal[,1]
                              mean.rank = NA
                              ret = data.frame(tmp.names, ratio, mean.rank, length(z$EG.StrippedSequence))
                              names(ret) = c("peptide1", "peptide2", "ratio", "mean.rank", "peptide.count")
                              return(ret)
                              
                            } else {
                              ret = data.frame(z$EG.StrippedSequence,z$EG.StrippedSequence, 1, 1, length(z$EG.StrippedSequence))
                              names(ret) = c("peptide1", "peptide2", "ratio", "mean.rank", "peptide.count")
                              return(ret)
                            }
                          })

peptides.LIT.ratios = tbl_df(peptides.LIT.ratios)
peptides.LIT.ratios$pair = paste(peptides.LIT.ratios$peptide1, peptides.LIT.ratios$peptide2, sep=".")

peptides.LIT.ratios.CV = peptides.LIT.ratios %>% 
                         group_by(ORF, pair) %>% dplyr::summarize(CV = sd(ratio, na.rm=T)/mean(ratio, na.rm=T),
                                                                  rCV = mad(ratio, na.rm=T)/median(ratio, na.rm=T),
                                                                  median.rank = NA)

peptides.LIT.ratios.CV$pair = factor(peptides.LIT.ratios.CV$pair)
peptides.LIT.ratios.CV.long = melt(data.frame(peptides.LIT.ratios.CV), id.vars=c("ORF", "pair"))

p1 = ggplot(filter(peptides.LIT.ratios.CV.long, variable != "median.rank"), aes(x=value)) +
  geom_histogram(fill="white", color="black") +
  facet_wrap(~variable, scales="free") +
  theme(aspect.ratio = 1)

pep.count = 10
selectedORF = as.character(peptides.LIT.ratios$ORF[which(peptides.LIT.ratios$peptide.count == pep.count)[1]])


p2 = ggplot(filter(peptides.LIT.ratios, ORF == selectedORF ), aes(x=pair, y=log(ratio,2))) +
  geom_boxplot() +
  ggtitle(paste(selectedORF,"with", pep.count, "peptides", sep=" ")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

g = arrangeGrob(p2, p1)

plots.list = lappend(plots.list, g)

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 


