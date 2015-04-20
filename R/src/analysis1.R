rm(list=ls())
source("./R/functions.R")
source("./R/boot.R")

plots.list = list()
fun_name = "analysis1"

load("./R/objects/proteins.matrix.f.deseq.combat.RData")
load("./R/objects/proteins.matrix.f.sva.adjusted.RData")

load("./R/objects/sample_exp.map.RData")
load("./R/objects/peptides.matrix.f.RData")


proteins.matrix = proteins.matrix.f.deseq.combat


tmp.wide = as.data.frame(proteins.matrix)
tmp.wide$ORF = rownames(proteins.matrix)
proteins.long = melt(tmp.wide, id.vars="ORF")
names(proteins.long) = c("ORF", "variable", "value")

proteins.long = proteins.long %>% extract(variable, 
                                          into=c("R.Label", "batch_date", "batch.exp.n", "batch"), 
                                          regex="(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")

col_names <- names(proteins.long)[names(proteins.long) != "value"]
proteins.long[,col_names] <- lapply(proteins.long[,col_names] , factor)
proteins.long = tbl_df(proteins.long)
proteins.long$KO = factor(sample_exp.map$ORF[match(proteins.long$R.Label, sample_exp.map$SampleName)])

#proteins.long   = proteins.deseq.combat.long

pattern.p = "(.*?)_([0-9]+_[0-9]+_[0-9]+|[A-Za-z]?|[A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$"
matches = stringr::str_match_all(pattern=pattern.p, colnames(proteins.matrix))

stopifnot(sum(lapply(matches,length)!=0) == ncol(proteins.matrix))
pheno = data.frame(matrix(unlist(matches), ncol=length(matches[[1]]), byrow=T))
colnames(pheno) = c("name", "R.Label", "batch_date", "batch.exp.n", "batch" )
rownames(pheno) = colnames(proteins.matrix)
pheno$ORF = droplevels(sample_exp.map$ORF[match(pheno$R.Label, sample_exp.map$SampleName)])
pheno$ORF[pheno$R.Label == "KL_Try_027_c"] = "WT"
pheno$batch.exp.n[pheno$R.Label == "KL_Try_027_c"] = 5

pheno = droplevels(pheno[pheno$ORF != "none",]) #removing none-phenotype samples
proteins.matrix.f = proteins.matrix[,match(pheno$name, colnames(proteins.matrix))]


X = model.matrix(~pheno$ORF + 0)
colnames(X) = levels(pheno$ORF)

reference = "WT"
matrix = proteins.matrix.f

lm.fit_model = lmFit(matrix, X)
ph = unique(as.character(pheno$ORF))
contrasts = paste0( ph[ph !=reference] ,"-", reference)  

mc = makeContrasts(contrasts=contrasts, levels=X)    
c.fit = contrasts.fit(lm.fit_model, mc)
eb = eBayes(c.fit)
  
eb$coefficients
folds = rowFolds(data=exp(matrix), groups=pheno$ORF, reference=reference)
folds = log(folds, 2)

#folds_tmp = melt(as.matrix(folds), id.vars="row.names")

#merging results
folds_tmp = melt(eb$coefficients, id.vars="row.names")
#folds_tmp$contrasts = factor(paste(folds_tmp$contrasts, "-", reference, sep=""))
pvals_tmp = melt(eb$p.value, id.vars="row.names")


names(folds_tmp) = c("ORF", "contrasts", "logFC")
names(pvals_tmp) = c("ORF", "contrasts", "p.value")

folds_tmp$contrasts = factor(folds_tmp$contrasts)
pvals_tmp$contrasts = factor(pvals_tmp$contrasts)


proteins.FC = merge(folds_tmp, pvals_tmp, all=T,
                    by=c("ORF", "contrasts"))

##multiple testing correction
proteins.FC$p.value_BH = p.adjust(proteins.FC$p.value, method="BH")
proteins.FC$p.value_bonferroni = p.adjust(proteins.FC$p.value, method="bonferroni")

proteins.FC$KO = sub(x = proteins.FC$contrasts, pattern=paste("(.*?)-", reference, sep=""), replacement="\\1")
proteins.FC$reference = reference

file_name = "proteins.FC.RData"
file_path = paste(output_dir, file_name, sep="/")
save(proteins.FC, file=file_path)


lb = -4
ub = 4
toPlot = proteins.FC
toPlot[toPlot$logFC < 0,]$logFC = ifelse(toPlot[toPlot$logFC < 0,]$logFC < lb, lb, toPlot[toPlot$logFC < 0,]$logFC)
toPlot[toPlot$logFC > 0,]$logFC = ifelse(toPlot[toPlot$logFC > 0,]$logFC > ub, ub, toPlot[toPlot$logFC > 0,]$logFC)

FC_thr = 1
pval_thr = 0.05
proteins.FC.stats = data.frame( ratio_sign = round(sum(proteins.FC$logFC < -FC_thr & proteins.FC$p.value_BH < pval_thr)/sum(proteins.FC$logFC > FC_thr & proteins.FC$p.value_BH < pval_thr),2),
                                ratio =  round(sum(proteins.FC$logFC < 0)/sum(proteins.FC$logFC > 0),2),
                                n_prot = length(unique(proteins.FC$ORF)),
                                n_sign = round(sum(proteins.FC$p.value_BH<0.05)/length(unique(proteins.FC$contrasts)),2),
                                x_min = -4,
                                y_max =max(-log10(proteins.FC$p.value_BH)))

toPlot$sign = ifelse(abs(toPlot$logFC) >= FC_thr & toPlot$p.value_BH < pval_thr, 1,0)


p1 = ggplot(toPlot, aes(y=-log10(p.value_BH), x=logFC)) +
         geom_point(aes(color = sign), alpha=0.5) + 
         xlim(c(lb,ub)) +
         geom_hline(y=-log(pval_thr,10),linetype=3) +
         geom_vline(x=c(FC_thr,-FC_thr),linetype=3) +
         geom_text(data=proteins.FC.stats, aes(x=x_min, y=y_max, hjust=0, vjust=1,
                                             label=paste(paste0("#_mean_sign = ", n_sign),
                                                         paste0("up/down_ratio = ",ratio),
                                                         paste0("ratio_sign = ",   ratio_sign),
                                                         paste0("#_prot = ", n_prot), sep="\n"))) +
        xlab(paste("Log2(fold-change)")) +
        theme(aspect.ratio = 1, legend.position = "none", text=element_text(size=18))

p2 = ggplot(toPlot, aes(x=logFC)) +
            geom_histogram(colour = "white", fill = "black", binwidth = 0.15) +
            xlab(paste("Log2(fold-change)")) +
            theme(aspect.ratio = 1, text=element_text(size=18))

g = arrangeGrob(p1,p2 , ncol=2, main=paste("Reference:", reference))
plots.list = lappend(plots.list, g)


proteins.FC.f = proteins.FC[proteins.FC$KO %in% droplevels(unique(sample_exp.map$ORF[grep("Kinase", sample_exp.map$Type)])),]

p_thr = 0.05
fc_thr = 1
top = 150
tmp.f = table(droplevels(filter(proteins.FC.f, p.value_BH < p_thr, (logFC < fc_thr | logFC > fc_thr))$ORF))
selected = names(sort(-tmp.f)[1:top])

proteins.FC.f.selected = droplevels(proteins.FC.f[proteins.FC.f$ORF %in% selected,])
proteins.FC.wide = droplevels(dcast(proteins.FC.f.selected, formula=ORF~KO, value.var="logFC"))


# cor_dist = as.dist(1 - abs(cor(proteins.FC.wide[,-1])))
# toClust = t(proteins.FC.wide[,-1])
# colnames(toClust) = proteins.FC.wide$ORF
# cl = NbClust(data = t(proteins.FC.wide[,-1]), diss=cor_dist, distance=NULL, method="ward.D2", index=c("ch"))


cl_res = ConsensusClusterPlus(as.matrix(proteins.FC.wide[,-1]), maxK=6,
                              reps=1000, pItem=0.8, pFeature=1,finalLinkage="complete",
                              clusterAlg="hc",distance="spearman",seed=123 )

# cl_res2 = ConsensusClusterPlus(as.matrix(t(proteins.FC.wide[,-1])), maxK=20,
#                               reps=1000, pItem=0.8, pFeature=1,title=title, finalLinkage="complete",
#                               clusterAlg="hc",distance="euclidean",seed=123 )



cl_V = 3
#cl_H = 2

annotation = data.frame(batch_date = sample_exp.map$date[match(colnames(proteins.FC.wide)[-1], sample_exp.map$ORF)])
rownames(annotation) = colnames(proteins.FC.wide)[-1]
annotation$clusters = factor(cl_res[[cl_V]]$consensusClass)

col_breaks = c(-4,-3,-2,-1,1,2,3,4)
#col_breaks = c(0)
file_name = paste("heatmap_proteins.p_thr.fc_thr.top", p_thr, fc_thr, top, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")

aheatmap(proteins.FC.wide[,-1], annCol=annotation, 
         #annRow=factor(cl_res2[[cl_H]]$consensusClass),
         Colv=cl_res[[cl_V]]$consensusTree,
         #Rowv=cl_res2[[cl_H]]$consensusTree,
         color = paste("-RdBu", length(col_breaks)-1, sep=":"),
         breaks = col_breaks,
         filename=file_path)


proteins.FC.clusters = droplevels(merge(proteins.FC, data.frame(KO=rownames(annotation), cl = annotation$clusters), by="KO"))


##-- cluster enrichementes ----


load("./R/objects/pathway2orf.RData")
load("./R/objects/pathways.RData")
load("./R/objects/GO_slim.function.RData")
load("./R/objects/GO_slim.process.RData")
load("./R/objects/GO_slim.compartment.RData")

enrichments = list()
for (i in unique(proteins.FC.clusters$cl)) {
  for (j in c("up","down")) {
    tmp.selected = c()
    cluster = i
#     cluster = 1
#     j = "down"
    if (j == "up") {
      tmp.selected = table(droplevels(filter(proteins.FC.clusters, p.value_BH < p_thr, (logFC > fc_thr), cl == cluster)$ORF))  
    } else {
      tmp.selected = table(droplevels(filter(proteins.FC.clusters, p.value_BH < p_thr, (logFC < -fc_thr), cl == cluster)$ORF))  
    }
    
    important_size = round(max(table(droplevels(filter(proteins.FC.clusters, cl == cluster)$ORF)))/3)
    
    orf_thr = names(tmp.selected[tmp.selected >= important_size])
    orf_universe = unique(as.character(droplevels(filter(proteins.FC.clusters, cl == cluster)$ORF)))
    
    pathway.counts = pathway_enrichments(orf_thr=orf_thr, orf_universe=orf_universe, pathway2orf=GO_slim.process)
    pathway.counts$cl = factor(cluster)
    pathway.counts$direction = factor(j)
    enrichments[[i]] = pathway.counts
  }
}

enrichments.clusters = do.call(rbind.data.frame, enrichments)

View(enrichments.clusters[enrichments.clusters$p.value<0.05,])

##---- protein variations ----

proteins.long.ko_mean = proteins.long %>% 
                        group_by(KO, ORF) %>%
                        summarise(signal = mean(value))

proteins.long.ko_mean = proteins.long.ko_mean[grep("none", proteins.long.ko_mean$KO, ignore.case=T, invert=T),]
     
proteins.ko_mean = dcast(proteins.long.ko_mean, ORF~KO, value.var="signal")

proteins.matrix.ko_mean = as.matrix(proteins.ko_mean[,-1])
rownames(proteins.matrix.ko_mean) = proteins.ko_mean$ORF

proteins.cv.variation = apply(proteins.matrix.ko_mean, 1, 
                      FUN=function(x) {
                        sqrt(exp(sd(x)^2)-1)        
                      } )


#proteins variation
toPlot=data.frame(ORF=names(proteins.cv.variation), 
                  cv.variation = proteins.cv.variation) 
p = ggplot(data=toPlot, aes(x=cv.variation)) +
      geom_histogram(binwidth=0.05, fill="white", colour="black") +
      xlab("Coefficient of variation of protein intensity across all KO") +
      theme(axis.title=element_text(size=20),
            axis.text=element_text(size=14),
            aspect.ratio = 1)
plots.list = lappend(plots.list, p)

upper = 0.95
lower = 0.05
signal_upper = names(proteins.cv.variation[proteins.cv.variation > quantile(proteins.cv.variation,upper)])
signal_lower = names(proteins.cv.variation[proteins.cv.variation < quantile(proteins.cv.variation,lower)])
universe  = names(proteins.cv.variation)


View(pathway_enrichments(orf_thr=signal_lower, orf_universe=universe, pathway2orf=GO_slim.compartment))
View(pathway_enrichments(orf_thr=signal_lower, orf_universe=universe, pathway2orf=pathway2orf))



proteins.z_scores = apply(proteins.matrix.ko_mean, 1, 
                     FUN=function(x) {
                       unlog_x = exp(x)
                       sd.unlog_x = sd(unlog_x,na.rm=T)
                       mean.unlog_x = mean(unlog_x)
                       
                       sd.log = sqrt(log(1+sd.unlog_x^2/mean.unlog_x^2))
                       mean.log = log(mean.unlog_x)-0.5*sd.log^2
                       z.log = (x - mean.log)/sd.log
                       return(z.log)
                     } )



#write.table(unique(colnames(proteins.z_scores)), file="proteins.txt", quote=F, row.names=F, col.names=F)
mutants.z_sums = apply(proteins.z_scores, 1,
                      FUN=function(x) {
                        sum(abs(x))
                      })

hist(mutants.z_sums, breaks=30, xlim=c(400,2000))
points(mutants.z_sums[names(mutants.z_sums)=="WT"], 2, col="red", pch=20)
points(mutants.z_sums[names(mutants.z_sums)=="YPL115C"], 2, col="red", pch=20)
points(mutants.z_sums[names(mutants.z_sums)=="YPL185W"], 2, col="red", pch=20)

text(mutants.z_sums[names(mutants.z_sums)=="WT"], 3, labels="WT")

## ---- comprarison with PPI network ----
load("./R/objects/BIOGRID.ORGANISM.RData")

BIOGRID.ORGANISM = tbl_df(BIOGRID.ORGANISM)
BIOGRID.ORGANISM$Pubmed.ID = factor(BIOGRID.ORGANISM$Pubmed.ID)
interactions = BIOGRID.ORGANISM %>% select(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B,Experimental.System, Experimental.System.Type, Pubmed.ID)

interactions.p = filter(interactions, Experimental.System.Type == "physical") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
interactions.g = filter(interactions, Experimental.System.Type == "genetic") %>% distinct(Systematic.Name.Interactor.A, Systematic.Name.Interactor.B)
G.phys = graph.data.frame(select(interactions.p, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)
G.gene = graph.data.frame(select(interactions.g, Systematic.Name.Interactor.A, Systematic.Name.Interactor.B), directed=F)

proteins.FC.sums = proteins.FC.clusters %>% group_by(KO) %>% summarize(perturbation = sum(abs(logFC)))

tmp = data.frame(KO = names(mutants.z_sums), 
           perturb = proteins.FC.sums$perturbation[match(names(mutants.z_sums), proteins.FC.sums$KO)],
           z_sums = mutants.z_sums, 
           degree.phys = degree(G.phys)[match(names(mutants.z_sums), names(degree(G.phys)))],
           degree.gene = degree(G.gene)[match(names(mutants.z_sums), names(degree(G.gene)))])


proteins.z_scores.long = melt(proteins.z_scores, id.vars=row.names)
names(proteins.z_scores.long) = c("KO", "ORF", "value")
all_enrichments = ddply(proteins.z_scores.long, .(KO), 
      .fun = function(z) {
        #z = proteins.z_scores.long[proteins.z_scores.long$KO == "WT",]
        signal = z$ORF[abs(z$value) >= 1.98]
        universe =  z$ORF
        
        tmp.p = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.process)
        tmp.p$type = "process"
        tmp.f = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.function)
        tmp.f$type = "function"
        tmp.c = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=GO_slim.compartment)
        tmp.c$type = "compartment"
        tmp.path = pathway_enrichments(orf_thr=signal, orf_universe=universe, pathway2orf=pathway2orf)
        tmp.path$type = "kegg_pathways"
        
        return(rbind(tmp.p, tmp.c, tmp.c, tmp.path))
        
      })

tmp$p.adf_all = tmp$p.value
annotation$clusters = NULL

file_name = paste(fun_name, "report.pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
save_plots(plots.list, filename=file_path, type="l") 




# #enrichments in clusters
# proteins.FC
# #modules
# module2orf.universe = droplevels(module2orf[module2orf$ORF %in% orf_universe,])
# module2orf.thr      = droplevels(module2orf[module2orf$ORF %in% orf_thr,])
# 
# module2orf.thr$ORF = factor(module2orf.thr$ORF)
# module2orf.universe$ORF = factor(module2orf.universe$ORF)
# 
# module.signal = module2orf.thr %>% group_by(md) %>% dplyr::summarise(count = n())
# module.universe = module2orf.universe %>% group_by(md) %>% dplyr::summarise(count = n())
# 
# 
# module.merged = merge(module.signal, module.universe, by="md", suffixes=c(".signal", ".universe"))
# 
# total.universe = length(unique(module2orf.universe$ORF))
# total.signal   = length(unique(module2orf.thr$ORF))
# 
# module.counts = module.merged %>% group_by(md) %>%
#   mutate(notSignal.inPath = count.universe - count.signal,
#          Signal.notinPath = total.signal - count.signal,
#          notSignal.notinPath = total.universe - notSignal.inPath )
# 
# module.counts = module.merged %>% group_by(md) %>%
#   mutate(notSignal.inPath = count.universe - count.signal,
#          Signal.notinPath = total.signal - count.signal,
#          notSignal.notinPath = total.universe - notSignal.inPath )
# 
# 
# module.counts$description = modules$description[match(module.counts$md, modules$md)]
# 
# module.counts =  module.counts %>% group_by(md) %>%             
#   mutate(p.value = 1 - phyper(q=count.signal-1, m=total.signal, n=total.universe-total.signal,  k=count.universe, lower.tail=T))
# 
# module.counts$p.adj = p.adjust(module.counts$p.value, method="BH")
# names(module.counts)[1] = "md|pathway"
# module.counts$cl = cluster


