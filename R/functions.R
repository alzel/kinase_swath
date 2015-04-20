
source("./R/boot.R")

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
    width = 11.69+0.1*11.69
    height = 8.27+0.1*8.27
  }
  
  pdf(filename, width=width, height=height)
  for (i in 1:length(plots.list)) {
    print(plots.list[[i]])
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
  
  stopifnot(length(levels(groups)) >=2)  
    
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


