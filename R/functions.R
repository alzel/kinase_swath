
source("./R/boot.R")
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


cor_all_within = function (matrix, groups) {
  within = data.frame()
  for (l in levels(groups)) {
    to_compare = which(groups==l)
    i = 1
    j = 1
    while (i <= length(to_compare)){
      j = i + 1
      while (j <= length(to_compare)) {
        ind1 = to_compare[i]
        ind2 = to_compare[j]
        tmp = data.frame(cor = cor(matrix[,ind1], matrix[,ind2], use="pairwise.complete"), 
                         spearman = cor(matrix[,ind1], matrix[,ind2], use="pairwise.complete", method = "spearman") ,  
                         group1=l, group2=NA, perspective="within")
        within = rbind(within, tmp)        
        j = j + 1
      }
      i = i + 1
    }
  }
  return(within)
}

plot_all_within = function (matrix, groups) {
  within = data.frame()
  for (l in levels(groups)) {
    to_compare = which(groups==l)
    i = 1
    j = 1
    while (i <= length(to_compare)){
      j = i + 1
      while (j <= length(to_compare)) {
        ind1 = to_compare[i]
        ind2 = to_compare[j]
        plot(matrix[,ind1], matrix[,ind2])
        tmp = data.frame(cor = cor(matrix[,ind1], matrix[,ind2],use="pairwise.complete"), 
                         spearman = cor(matrix[,k], matrix[,l],use="pairwise.complete", method="spearman") ,  
                         group1=l, group2=NA, perspective="within")
        within = rbind(within, tmp)        
        j = j + 1
      }
      i = i + 1
    }
  }
  return(within)
}

cor_all_between = function (matrix, groups) {
  between = data.frame()
  l_groups = levels(groups)
  
  i = 1 #index group 1
  j = 1 #index group 2
  while (i <= length(l_groups)){
    to_compare1 = which(groups==l_groups[i])
    j = i + 1
    while (j <= length(l_groups)) {
      to_compare2 = which(groups==l_groups[j])
      for (k in to_compare1) {
        for (l in to_compare2) {
          tmp = data.frame(cor = cor(matrix[,k], matrix[,l],use="pairwise.complete"),
                           spearman = cor(matrix[,k], matrix[,l],use="pairwise.complete", method="spearman"),  
                           group1=groups[k], group2=groups[l], perspective="between")
          between = rbind(between, tmp)
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  return(between)
}

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}


cor_peptides = function(data) {
  #computes peptides correlations among proteins
  #peptide_cor = data.table()
  peptide_cor = list()
  for (protein in levels(data$ORF)) {
    
    
    to_compare = which(data$ORF %in% protein)
    
    if (length(to_compare) == 1) {
      x = as.numeric(data[to_compare[1],c(-1,-2)])
      tmp = data.frame( protein = protein,
                        cor=NA, cor_log2 = NA, spearman = NA,
                        peptide1=data[to_compare[1],2], 
                        peptide2=data[to_compare[1],2], 
                        pair = paste(data[to_compare[1],2], data[to_compare[1],2], sep = "."),
                        nr_peptides = length(to_compare),
                        prevalence_p1 = 1 - sum(is.na(x))/length(x),
                        prevalence_p2 = 1 - sum(is.na(x))/length(x))
            
      peptide_cor[[length(peptide_cor) + 1]] = tmp
      
      #l = list(peptide_cor, as.data.table(tmp))
      #peptide_cor = rbindlist(l)
    } else {
      i = 1
      j = 1
      while (i <= length(to_compare)){
        j = i + 1
        while (j <= length(to_compare)) {
          ind1 = to_compare[i]
          ind2 = to_compare[j]
          x = as.numeric(data[ind1,c(-1,-2)])
          y = as.numeric(data[ind2,c(-1,-2)])
          
          cor = cor(x,y, use="pairwise.complete")
          cor_log2 = cor(log2(x),log2(y), use="pairwise.complete")
          spearman = cor(x, y, use="pairwise.complete", method="spearman")
          tmp = data.frame(protein = protein,
                           cor=cor, cor_log2 = cor_log2, spearman = spearman, 
                           peptide1=data[ind1,2], 
                           peptide2=data[ind2,2], 
                           pair = paste(data[ind1,2],data[ind2,2], sep = "."),
                           nr_peptides = length(to_compare),
                           prevalence_p1 = 1 - sum(is.na(x))/length(x),
                           prevalence_p2 = 1 - sum(is.na(y))/length(y))
#           l = list(peptide_cor, as.data.table(tmp))
#           peptide_cor = rbindlist(l)
          peptide_cor[[length(peptide_cor) + 1]] = tmp
          j = j + 1
        }
        i = i + 1
      }
    }
  }
  
  return(do.call(rbind.data.frame(peptide_cor)))
}

plot_allMA = function (matrix, groups) {
  ##TODO: estimate limits for plot(0) data is out of range 
  
  matrix = proteins.matrix.f[,-1]
  groups
  plot(0, type="n",axes=T, ylab="A, log2(replicate1/replicate2)", 
       xlab="M, (log2(replicate1)+log2(replicate2))/2 ",  cex.lab=1.33)
  for (l in levels(groups)) {
    to_compare = which(groups==l)
    i = 1
    j = 1
    while (i <= length(to_compare)){
      j = i + 1
      while (j <= length(to_compare)) {
        ind1 = to_compare[i]
        ind2 = to_compare[j]
        
        M=log2(matrix[,ind1])-log2(matrix[,ind2])
        A=(log2(matrix[,ind1])+log2(matrix[,ind2]))/2
        points(A,M, col=rgb(0,0,0,0.25))
        fit<-loess(M~A)
        ord <- order(A)
        lines(A[ord],fit$fitted[ord], col=rainbow(length(levels(groups)))[which(levels(groups)==l)],lwd=2)
        j = j + 1
      }
      i = i + 1
    }
  }
  abline(h=0)
  legend("topright", pch=rep(15,length(levels(groups))), legend=as.vector(levels(groups)),
         col=rainbow(length(levels(groups)))[1:length(levels(groups))],
         cex=0.66)
  
}



# cluster_peptides = function(peptides.matrix.all.filtered, cor_thr) {
#   
#   tmp = ddply(peptides.matrix.all.filtered[complete.cases(peptides.matrix.all.filtered),], .(protein.Accession, type), 
#              .fun = function(x) {
#                tmp_data = x[complete.cases(x),]
#                if(nrow(tmp_data) == 1) {
#                  return(data.frame(protein.Accession = tmp_data$protein.Accession,
#                                    type = tmp_data$type,
#                                    peptide.seq = tmp_data$peptide.seq,
#                                    clusters = c(1),
#                                    selected_cluster = c("Y") ))
#                } else {
#                  #tmp_data = peptides.matrix.all[peptides.matrix.all$type == "quantiles" & peptides.matrix.all$protein.Accession=="P00815",]
#                  tmp_matrix = t(tmp_data[,-c(1,2,3)])
#                  #tmp_matrix
#                  colnames(tmp_matrix) = as.character(tmp_data$peptide.seq)
#                  d = as.dist(1 - cor(tmp_matrix, method="spearman"))
#                  h = hclust(d, method="complete")
#                  plot(h)
#                  clusters = cutree(h, h = (1 - cor_thr))
#                  peptide.seq = names(clusters)
#                  cuts.t = table(clusters)
#                  cuts.t
#                  
#                  selected_cluster = ifelse(clusters == names(cuts.t[match(max(cuts.t), cuts.t)]), "Y", "N") #selecting best cluster
#                  selected_cluster
#                  names(clusters) = NULL
#                  names(selected_cluster) = NULL
#                  return(data.frame(protein.Accession = tmp_data$protein.Accession,
#                                    type = tmp_data$type,
#                                    peptide.seq = peptide.seq,
#                                    clusters = clusters,
#                                    selected_cluster = selected_cluster)) 
#                }
#              })
#     return (tmp)
# }

cluster_peptides = function(peptides.matrix.all.filtered, cor_thr) {
  
  #tmp = ddply(peptides.matrix.all.filtered[complete.cases(peptides.matrix.all.filtered),], .(ORF), 
  tmp = ddply(peptides.matrix.all.filtered, .(ORF), 
              .fun = function(x) {
                tmp_data <<- x
                
                if(nrow(tmp_data) == 1) {
                  return(data.frame(ORF = tmp_data$ORF,
                                    EG.StrippedSequence = tmp_data$EG.StrippedSequence,
                                    clusters = c(1),
                                    selected_cluster = c("Y") ))
                } else {
                  
                  tmp_matrix <<- tmp_data[,-c(1,2)]
                  tmp_matrix = t(tmp_data[,-c(1,2)])
                  colnames(tmp_matrix) = as.character(tmp_data$EG.StrippedSequence)
                  
                  tmp_cor = cor(tmp_matrix, method="spearman", use="pairwise.complete.obs")
                  
                  shorths = apply(tmp_cor,1, FUN=function(x) {genefilter::shorth(x[x != 1 & !is.na(x)], na.rm=T,tie.limit=0.5)})
                  
                  selected_cluster = ifelse(round(shorths,2) >= cor_thr,  "Y", "N") #selecting best cluster
                  names(shorths) = NULL
                  names(selected_cluster) = NULL
                  
                  return(data.frame(ORF = tmp_data$ORF,
                                    EG.StrippedSequence = rownames(tmp_cor),
                                    clusters = shorths,
                                    selected_cluster = selected_cluster)) 
                }
              })
  return (tmp)
}


Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(p.val)
}

save_plots = function(plots.list, filename, type="p") {
  width=8.27+0.1*8.27
  height=11.69+0.1*11.69
  
  if (type == "l") {
    width = 11.69+0.15*11.69
    height = 8.27+0.15*8.27
  }
  
  pdf(filename, width=width, height=height)
  for (i in 1:length(plots.list)) {
    #print(plots.list[[i]])
    plot(plots.list[[i]])
  }
  dev.off()
}

plot_pcacor = function(x,y,data, xlab, ylab) {
  plot(x, y, col=rainbow(length(levels(data$pheno)))[as.numeric(data$pheno)], pch=16, xlab=xlab, ylab=ylab)
  tmp = cor.test(x, y)
  coef = tmp$estimate
  pval = tmp$p.value
  legend("topright", legend=c(paste(paste("R^2 = ", round(coef^2,2)),paste("P = ", format(pval,digits = 2, scientific=T)), sep="\n")),
         bty="n")
  text(x, y, labels=data$pheno, cex = 0.66)
}


geosd <- function(x, na.rm = FALSE, ...)
{
  exp(sd(log(x, ...), na.rm = na.rm, ...))
}

geomean <- function(x, na.rm = FALSE, trim = 0, ...)
{
  exp(mean(log(x, ...), na.rm = na.rm, trim = trim, ...))
}


rowFolds = function(data, groups, reference) {
  
#   data = metabolites.matrix.f
#   reference = "WT"
#   groups = phenotypes.f$ORF
  
  #   reference = "PSEUDO"
  #   groups = factor(phenotypes$name)
  #   
  
  
  ref_ind = which(groups == reference)
  
  if (length(ref_ind) == 0 ) {
    stop(paste("No such reference:", reference))
  }
  
  stopifnot(length(unique(groups)) >=2)  
    
  tmp_mat=data.frame(zero = rep(0,nrow(data)))
  ref_row = rowMeans(data[,ref_ind],na.rm=T)
  
  for (i in levels(groups)) {
    if (i == reference) {
      next
    }
    
    condition_index = which(groups == i)
    ret=c()
    if (length(condition_index) == 1 ){
      ret = data[,condition_index]/ref_row
    } else {
      ret = rowMeans(data[,condition_index],na.rm=T)/ref_row  
    }
    tmp_mat = cbind(tmp_mat,ret)  
  }
  tmp_mat$zero = NULL
  colnames(tmp_mat) = levels(groups)[levels(groups) != reference]
  
  return(tmp_mat)  
}



sva.patched=function (dat, mod, mod0 = NULL, n.sv = NULL, method = c("irw",  
                                                                     "two-step"), vfilter = NULL, B = 5, numSVmethod = "be")  
{ 
  if (is.null(n.sv)) { 
    n.sv = num.sv(dat, mod, method = numSVmethod, vfilter = vfilter) 
  } 
  if (!is.null(vfilter)) { 
    if (vfilter < 100 | vfilter > dim(dat)[1]) { 
      stop(paste("The number of genes used in the analysis must be between 100 and",  
                 dim(dat)[1], "\n")) 
    } 
    tmpv = rowVars(dat) 
    ind = which(rank(-tmpv) < vfilter) 
    dat = dat[ind, ] 
  } 
  if (n.sv > 0) { 
    cat(paste("Number of significant surrogate variables is: ",  
              n.sv, "\n")) 
    method <- match.arg(method) 
    if (method == "two-step") { 
      return(twostepsva.build(dat = dat, mod = mod, n.sv = n.sv)) 
    } 
    if (method == "irw") { 
      # return(irwsva.build(dat = dat, mod = mod, mod0 = mod0, n.sv = n.sv, B = B)) 
      return(irwsva.build.patched(dat = dat, mod = mod, mod0 = mod0, n.sv = n.sv, B = B)) 
    } 
  } 
  else { 
    print("No significant surrogate variables") 
    return(list(sv = 0, pprob.gam = 0, pprob.b = 0, n.sv = 0)) 
  } 
} 

irwsva.build.patched<- 
  function (dat, mod, mod0 = NULL, n.sv, B = 5)  
  { 
    n <- ncol(dat) 
    m <- nrow(dat) 
    if (is.null(mod0)) { 
      mod0 <- mod[, 1] 
    } 
    Id <- diag(n) 
    resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%  
                        t(mod)) 
    uu <- eigen(t(resid) %*% resid) 
    vv <- uu$vectors 
    ndf <- n - dim(mod)[2] 
    pprob <- rep(1, m) 
    one <- rep(1, n) 
    Id <- diag(n) 
    df1 <- dim(mod)[2] + n.sv 
    df0 <- dim(mod0)[2] + n.sv 
    rm(resid) 
    cat(paste("Iteration (out of", B, "):")) 
    for (i in 1:B) { 
      mod.b <- cbind(mod, uu$vectors[, 1:n.sv]) 
      mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv]) 
      ptmp <- f.pvalue(dat, mod.b, mod0.b) 
      pprob.b <- (1 - sva:::edge.lfdr(ptmp)) 
      mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv]) 
      mod0.gam <- cbind(mod0) 
      ptmp <- f.pvalue(dat, mod.gam, mod0.gam) 
      pprob.gam <- (1 - sva:::edge.lfdr(ptmp)) 
      pprob <- pprob.gam * (1 - pprob.b) 
      dats <- dat * pprob 
      dats <- dats - rowMeans(dats) 
      uu <- eigen(t(dats) %*% dats) 
      cat(paste(i, " ")) 
    } 
    # sv = fast.svd(dats, tol = 0)$v[, 1:n.sv]  
    # Patch code 
    if(any(dats!=0)) sv = fast.svd(dats, tol = 0)$v[, 1:n.sv] 
    else sv=svd(dats)$v[, 1:n.sv] 
    
    retval <- list(sv = sv, pprob.gam = pprob.gam, pprob.b = pprob.b,  
                   n.sv = n.sv) 
    return(retval) 
}


pathway_enrichments = function(orf_thr, orf_universe, pathway2orf ) {
  
  pathway2orf.universe = droplevels(pathway2orf[pathway2orf$ORF %in% orf_universe,])
  pathway2orf.thr      = droplevels(pathway2orf[pathway2orf$ORF %in% orf_thr,])
  
  pathway2orf.thr$ORF = factor(pathway2orf.thr$ORF)
  pathway2orf.universe$ORF = factor(pathway2orf.universe$ORF)
  
  pathway.signal = pathway2orf.thr %>% group_by(pathway) %>% dplyr::summarise(count = n())
  pathway.universe = pathway2orf.universe %>% group_by(pathway) %>% dplyr::summarise(count = n())
  
  pathway.merged = merge(pathway.signal, pathway.universe, by="pathway", suffixes=c(".signal", ".universe"))
  pathway.merged[is.na(pathway.merged)] = 0
  
  total.universe = length(unique(pathway2orf.universe$ORF))
  total.signal   = length(unique(pathway2orf.thr$ORF))
  
  pathway.counts = pathway.merged %>% group_by(pathway) %>%
    mutate(notSignal.inPath = count.universe - count.signal,
           Signal.notinPath = total.signal - count.signal,
           notSignal.notinPath = total.universe - notSignal.inPath )
  
  pathway.counts = pathway.merged %>% group_by(pathway) %>%
    mutate(notSignal.inPath = count.universe - count.signal,
           Signal.notinPath = total.signal - count.signal,
           notSignal.notinPath = total.universe - notSignal.inPath )
  
  
  pathway.counts$description = pathway2orf$description[match(pathway.counts$pathway, pathway2orf$path)]
  
  pathway.counts =  pathway.counts %>% group_by(pathway) %>%             
    mutate(p.value = 1 - phyper(q=count.signal-1, m=total.signal, n=total.universe-total.signal,  k=count.universe, lower.tail=T))
  
  pathway.counts$p.adj = p.adjust(pathway.counts$p.value, method="BH")
  return(pathway.counts)
  
}

my_diagnostics = function(fit, filename = NULL, influence=T) {
  
  stopifnot(class(fit) == "lm")
  plots.list = list()
  
  myfortdata  = fortify(fit)
  fit.summary = summary(fit)
  colnames(myfortdata)[1] = "metabolite"
  outliers = outlierTest(fit)
  
  cooks_thrI = 4/(nrow(myfortdata) - length(fit$coefficients)-2)  
  
  myfortdata$isOutlier = ifelse(rownames(myfortdata) %in% names(outliers$bonf.p < 0.05), TRUE, FALSE)
  myfortdata$isInfluence = ifelse(myfortdata$.cooksd > cooks_thrI, TRUE, FALSE)
  myfortdata$label = rownames(myfortdata)
  myfortdata$rows=1:nrow(myfortdata)
  myfortdata$.y_star = stats::predict(object=fit)
  
  beta.changes = data.frame(((fit$coefficients + lm.influence(fit)$coef) - fit$coefficients )/fit$coefficients)
  beta.changes$sample = rownames(beta.changes)
  beta.changes.long = melt(beta.changes[,-1], id.vars="sample")
  
  X <- fit$model[,-1]
  y <- fit$model[,1]
  
  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
  
#  set.seed(123) 
#   CVs = c()
#   for (i in 1:1000) {
#     cv.results <- crossval(x=X,y=y,theta.fit=theta.fit,theta.predict=theta.predict, ngroup=10)
#     cv.r.squared = cor(y,cv.results$cv.fit)**2
#     CVs = c(CVs, cv.r.squared)
#   }
    
  p1 = ggplot(data = myfortdata, aes(x = .fitted, y = .resid)) +
    geom_hline(yintercept = 0, colour = "firebrick3") +
    geom_point(size=2) +
    geom_smooth(method="loess",se = FALSE, span=2/3) +
    ylab("Residuals") +
    xlab("Fitted values")
  plots.list = lappend(plots.list, p1)
  
  p2 = ggplot(data = myfortdata, aes(sample = .stdresid)) +
    stat_qq() +
    geom_abline(colour = "firebrick3") +
    ylab("Standartized residuals")
  plots.list = lappend(plots.list, p2)
  
  p3 = ggplot(data = myfortdata, aes(x = .fitted, y = sqrt(abs(.stdresid)))) +
    geom_point(size=2) +
    geom_smooth(method="loess",se = FALSE)+
    xlab("Fitted values")
  plots.list = lappend(plots.list, p3)
  
  p4 = ggplot(data = myfortdata, aes(x = .hat, y=.stdresid)) +
    geom_point(aes(size = .cooksd)) +
    geom_hline(yintercept=c(-2,0,2)) +
    geom_vline(xintercept=c(2*sum(myfortdata$.hat)/nrow(myfortdata),
                            3*sum(myfortdata$.hat)/nrow(myfortdata)))+
    scale_x_continuous("Leverage") +
    scale_y_continuous("Standardized residuals") +
    scale_size_area("Cook’s distance", max_size=6)
  
  if (influence) {
    p4 = p4 + geom_text(data=filter(myfortdata, isOutlier ==T),
                        hjust=0, vjust=0,
                        aes(x = .hat, y=.stdresid, label=label)) +
              geom_text(data=filter(myfortdata, isInfluence ==T),
                        hjust=0, vjust=1, col="red",
                        aes(x = .hat, y=.stdresid, label=label))
  }
  plots.list = lappend(plots.list, p4)
  
  p5 <- ggplot(myfortdata, aes(.hat, .cooksd)) +
        geom_point(size=3) +
        geom_smooth(method="loess",se=FALSE, span=2/3) +
        scale_x_continuous("Leverage") +
        scale_y_continuous("Cook's distance")
  
  if (influence) {
    p5 = p5 + geom_hline(yintercept=cooks_thrI) +
              geom_text(data=filter(myfortdata, isInfluence ==T),
                hjust=0, vjust=0, col="red",
                aes(x = .hat, y=.cooksd, label=label))
  }
    
  plots.list = lappend(plots.list, p5)
  
  p6 <- ggplot(myfortdata, aes(y=metabolite, x=.fitted)) +
                geom_point() +
                geom_smooth(method="lm", se=F) +
                scale_x_continuous(paste("Predicted",i)) +
                scale_y_continuous(paste("Measured",i)) +
                geom_text(x=min(myfortdata$.fitted)[1], 
                          y=max(myfortdata[,1])[1],
                          hjust=0, vjust=1, 
                          label=paste("R^2 = ",round(fit.summary$r.squared,2)), parse=T)
  plots.list = lappend(plots.list, p6)
  
  p7 = ggplot(beta.changes.long, aes(x=variable,y=value,2)) + 
              geom_boxplot() +
              geom_text(data=filter(beta.changes.long, abs(value) > 0.5), 
                        aes(x=variable, y=value, label=sample))
  plots.list = lappend(plots.list, p7)
  filter(myfortdata, isInfluence ==T)
  
  p8 <-  ggplot(myfortdata, aes(rows, .cooksd, ymin=0, ymax=.cooksd)) +
                geom_point() + geom_linerange() +
                scale_x_continuous("Observation Number") +
                scale_y_continuous("Cook's distance")
  if (sum(myfortdata$isInfluence ==T) > 0 ) {
    
   p8 <- p8 + geom_text(data=filter(myfortdata, isInfluence ==T),
                        hjust=0, vjust=0, col="red",
                        aes(x = rows, y=.cooksd, label=label))                
    
  }
    
  plots.list = lappend(plots.list, p8)
  
#   p9 = ggplot(data=CVs, aes(x=CVs)) + geom_density()
#   plots.list = lappend(plots.list, p9)
  
  return(plots.list)
}

theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())


splinesfit <- function(i,x,y, subset=NULL) {
  return(AIC(lm(y~bs(x,df=i, subset))))
} 




splines.normalization<-function(x,time,subset=NULL, best = NULL, progress=TRUE, adjust="median",...){
  #subset = logical specifying which subset of the data to be used for fitting
  #adjust = function used to adjust post normalized data statistic to that of the pre normalized
  #best = vector of best df for spline fits  
  if (progress == TRUE){ pb <- txtProgressBar(min = 0, max = ncol(x), style = 3)} else {pb<-NULL}
  res<-do.call("cbind",lapply(1:ncol(x),function(i){
    tmp.y = x[,i]
    df.time = time[subset]
    df.tmp.y = tmp.y[subset]
    DF.bs = ifelse(is.null(best), as.integer(optimize(splinesfit, df.time , df.tmp.y, interval = c(1,15))$minimum),
                   best[i])
    fit<- lm(tmp.y ~ bs(time, df = DF.bs), subset=subset, ...)
    pred<-predict(fit,data.frame(time=time))
    if (progress == TRUE){setTxtProgressBar(pb, i)}
    return(tmp.y-pred) # residuals for train and test
  }))
  
  if (progress == TRUE){close(pb)}
  if(!is.null(adjust)){
    scale1<-apply(x,2,adjust,na.rm=TRUE)
    
    tmp<-sweep(res,2,apply(res,2,min,na.rm=TRUE),"-") # get rid of negative values
    mins<-apply(x,2,min,na.rm=TRUE)
    tmp<-sweep(tmp,2,mins,"+")
    scale2<-apply(tmp,2,adjust,na.rm=TRUE)
    adjustment<-scale1/scale2
    res<-sweep(tmp,2,adjustment,"*")
  }
  colnames(res) = colnames(x)
  return(res)
}


cor_heatmap = function (clean.data) {
  #clean.data - matrix or data.frame of numbers 
  correlations = cor(clean.data)
  hc = hclust(as.dist(correlations))
  correlations.long = melt(correlations)
  
  ## set color representation for specific values of the data distribution
  quantile_range <- quantile(correlations, probs = seq(0, 1, 0.2))
  
  ## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)
  
  ## prepare label text (use two adjacent values for range text)
  label_text <- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
  
  ## discretize matrix; this is the most important step, where for each value we find category of predefined ranges (modify probs argument of quantile to detail the colors)
  mod_mat <- matrix(findInterval(correlations, quantile_range, all.inside = TRUE), nrow = nrow(correlations))
  colnames(mod_mat) = colnames(correlations)
  rownames(mod_mat) = rownames(correlations)
  
  toPlot = merge(correlations.long, melt(mod_mat), by=c("X1", "X2"))
  toPlot$X1 = factor(toPlot$X1, levels=rownames(correlations)[hc$order])
  toPlot$X2 = factor(toPlot$X2, levels=rownames(correlations)[hc$order])
  
  ## output the graphics
  pheat = ggplot(toPlot, aes(x = X1, y = X2, fill = factor(value.y))) +
    geom_tile(color = "black") +
    geom_text(aes(label=round(value.x,2)),size=5) +
    scale_fill_manual(values = color_palette, name = "", labels = label_text) +
    theme_grey() +
    theme_change
  return(pheat)
}


repeatedCV = function(fit, repeats = 100) {
  
  #cross-valitation of all data
  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 
  
  X <- fit$model[,-1]
  y <- fit$model[,1]
  
  CVs = c()
  for (tmp.i in 1:repeats) {
    cv.results <- crossval(x=X, y=y, 
                           theta.fit=theta.fit, 
                           theta.predict=theta.predict, ngroup=10)
    cv.r.squared = cor(y,cv.results$cv.fit)**2
    CVs = c(CVs, cv.r.squared)
  }
  
  CVs = data.frame(CVs)  
  return(CVs)
}

getFC_thr = function(proteins.matrix = proteins.matrix.combat, pval_thr = 0.01) {
  
#   load("./R/objects/proteins.matrix.combat.RData")
#   proteins.matrix = proteins.matrix.combat
#   pval_thr = 0.01
#   
  #checking WT samples to define FC
  
  exp_metadata$aquisition_date.str = as.POSIXct(strftime(exp_metadata$aquisition_date, format="%Y-%m-%d %H:%M:%S"))
  cl = pam(exp_metadata$aquisition_date.str, 7)
  exp_metadata$batch_kmeans = cl$clustering
  
  pheno_wt = as.data.frame(exp_metadata[match(colnames(proteins.matrix), exp_metadata$sample_name),])
  pheno_wt = pheno_wt[pheno_wt$type == "Standard Mix",]
  pheno_wt = pheno_wt[pheno_wt$batch_kmeans %in% names(table(pheno_wt$batch_kmeans))[table(pheno_wt$batch_kmeans)  >= 3],]
  
  #plot(exp_metadata$aquisition_date, exp_metadata$batch_kmeans)
  
  #points(pheno_wt$aquisition_date, pheno_wt$batch_kmeans, col="red")
  
  proteins.matrix.WT = proteins.matrix[,match(pheno_wt$sample_name, colnames(proteins.matrix))]
  
  #   s = prcomp(t(proteins.matrix.WT))
  #   plot(s$x[,c(1,2)], col=pheno_wt$batch_kmeans)
  
  pheno_wt$factor = factor(paste(pheno_wt$ORF, pheno_wt$batch_kmeans, sep="."))
  
  X = model.matrix(~factor + 0, data=pheno_wt)
  colnames(X) = levels(pheno_wt$factor)
  
  tbl.tmp = table(pheno_wt$factor)
  reference = names(tbl.tmp)[tbl.tmp == max(tbl.tmp)][1]
  
  matrix = proteins.matrix.WT
  
  lm.fit_model = lmFit(matrix, X)
  ph = unique(as.character(pheno_wt$factor))
  contrasts = paste0( ph[ph !=reference] ,"-", reference)  
  
  mc = makeContrasts(contrasts=contrasts, levels=X)    
  c.fit = contrasts.fit(lm.fit_model, mc)
  eb = eBayes(c.fit)
  
  folds = rowFolds(data=exp(matrix), groups=pheno_wt$factor, reference=reference)
  folds = log(folds, 2)
  
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
  
  data = abs(proteins.FC$logFC[proteins.FC$p.value_BH < pval_thr])
  
  file_name = paste(fun_name, "getFC_thr", "pdf", sep=".")
  file_path = paste(figures_dir, file_name, sep="/")
  
  pdf(file_path, paper="a4")
  par(pty="s")
  hist(data, breaks=50, main="Expected fold changes in Standatd Mix")
  fc_thr = median(data)
  
  abline(v=fc_thr, lty=2)
  legend("topleft", bg=NULL, bty="n", 
         legend=paste("Median FC=", round(fc_thr,2)))
  dev.off()
  
  return(abs(fc_thr))
}


library(Matrix)
jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}


sample_removal <- function(fragments.tmp, fdr_thr1) {
  
  fragments.tmp[,-1] <- ifelse(fragments.tmp[,-1] < fdr_thr1, 1, 0) 
  fragments.tmp[is.na(fragments.tmp)] = 0
  fragments.tmp.matrix <- as.matrix(fragments.tmp[,-1])
  rownames(fragments.tmp.matrix) <- fragments.tmp$value
  
  sample_names <- colnames(fragments.tmp.matrix)
  tmp.matrix = as.matrix(jaccard(t(fragments.tmp.matrix)))
  diag(tmp.matrix) = NA
  
  tmp.list = list()
  
  rank_order = order(rank(rowMeans(tmp.matrix, na.rm = T)))
  for (i in 1:(nrow(tmp.matrix)-2)) {
    idx_remove <- rank_order[1:i]
    fragments.tmp.matrix.f <- fragments.tmp.matrix[, -idx_remove]
    fragments.tmp.matrix.f[fragments.tmp.matrix.f == 0] <- NA
    
    #z <<- data.frame(filename = sample_names[i], N = nrow(na.omit(fragments.tmp.matrix.f)))
    fragments.tmp.matrix.f
    tmp.list[[i]] <- data.frame(filename = sample_names[i], N = nrow(na.omit(fragments.tmp.matrix.f)), i=i)
  }
  #adding if dataset left untouched
  idx <- length(tmp.list)+1
  fragments.tmp.matrix.f = fragments.tmp.matrix
  fragments.tmp.matrix.f[fragments.tmp.matrix.f == 0] <- NA
  tmp.list[[idx]] <- data.frame(filename = "none", N = nrow(na.omit(fragments.tmp.matrix.f)), i=0)
  entity_removal <- do.call(rbind.data.frame, tmp.list)
  
  entity_removal$z_score = (entity_removal$N - mean(entity_removal$N))/sd(entity_removal$N)
  return(entity_removal)
}



untransform <- function(trans, data.transformed) {
  
  stopifnot(class(trans) == "preProcess")
  ignored = c()
  if (!is.null(trans$method)) {
    methods <- unlist(as.list(my_models$trans$call$method)[-1])  
  } else {
    methods <- names(trans$method)
    ignored <- trans$method$ignore
  }
  
  new_data <- data.transformed[, !colnames(data.transformed) %in% ignored]
  
  if (any(methods == "pca")) {
    new_data = as.matrix(new_data) %*% t(trans$rotation)
  }
  
  if (any(methods == "scale")) {
    new_data <- scale(new_data, center = FALSE , scale=1/trans$std)
  }
  
  if (any(methods == "center")) {
    new_data <- scale(new_data, center = -1 * trans$mean, scale=FALSE)
  }
  
  if (any(methods == "BoxCox")) {
    new_data <- sapply(names(trans$bc), 
                       FUN = function(i) {
                         inverse.BoxCoxTrans(object = trans$bc[[i]], newdata = new_data[,i])
                       })  
  } 
  
  new_data <- as.data.frame(new_data)
  return (new_data)
}

inverse.BoxCoxTrans <- function(object, newdata) {
  if(!is.vector(newdata) || !is.numeric(newdata)) stop("newdata should be a numeric vector")
  if(is.na(object$lambda)) return(newdata) 
  
  lambda <- object$lambda
  if(object$lambda < object$fudge & object$lambda > -object$fudge)
    lambda <- 0
  else if(object$lambda < 1 + object$fudge & object$lambda > 1 - object$fudge) {
    #lambda <- 1
    warning(paste("No transformation applied, lambda within the fudge", "lambda:", lambda, "tolerance:", object$fudge))
    return(newdata)
  }
  if(lambda == 0) exp(newdata) else (lambda*newdata + 1)^(1/lambda) 
}
