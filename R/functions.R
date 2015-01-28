
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
  peptide_cor = data.table()
  for (protein in levels(data$protein.Accession)) {
    
    to_compare = which(data$protein.Accession %in% protein)
    
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
      l = list(peptide_cor, as.data.table(tmp))
      peptide_cor = rbindlist(l)
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
          l = list(peptide_cor, as.data.table(tmp))
          peptide_cor = rbindlist(l)
          j = j + 1
        }
        i = i + 1
      }
    }
  }
  return(as.data.frame(peptide_cor))
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



cluster_peptides = function(peptides.matrix.all.filtered, cor_thr) {
  
  tmp = ddply(peptides.matrix.all.filtered[complete.cases(peptides.matrix.all.filtered),], .(protein.Accession), 
             .fun = function(x) {
               tmp_data = x[complete.cases(x),]
               if(nrow(tmp_data) == 1) {
                 return(data.frame(protein.Accession = tmp_data$protein.Accession,
                                   peptide.seq = tmp_data$peptide.seq,
                                   clusters = c(1),
                                   selected_cluster = c("Y") ))
               } else {
                 #tmp_data = peptides.matrix.all[peptides.matrix.all$protein.Accession=="P00815",]
                 tmp_matrix = t(tmp_data[,-c(1,2)])
                 
                 colnames(tmp_matrix) = as.character(tmp_data$peptide.seq)
                 d = as.dist(1 - cor(tmp_matrix, method="spearman"))
                 h = hclust(d, method="complete")
                 clusters = cutree(h, h = (1 - cor_thr))
                 peptide.seq = names(clusters)
                 cuts.t = table(clusters)
                 
                 selected_cluster = ifelse(clusters == names(cuts.t[match(max(cuts.t), cuts.t)]), "Y", "N") #selecting best cluster
                 names(clusters) = NULL
                 names(selected_cluster) = NULL
                 return(data.frame(protein.Accession = tmp_data$protein.Accession,
                                   peptide.seq = peptide.seq,
                                   clusters = clusters,
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


