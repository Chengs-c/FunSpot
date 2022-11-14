#' @param geno the expression matrix of scRNA data, a matrix  ngene * ncell
#' @param ST.matrix the expression matrix of ST data, a matrix ngene * nspot
#' @param cell_type_factor the annotation of scRNA data
#' @param min_count the minmum weight contribution
#' @param ct_mode the method that generate Q, ct_mode=c("median","mean")
#' @param nkonts number of knots
#' @param norder the order of the basis function
#' @param nfpca the numer of the principle component
#' @param basis_fun basis function
#' @param marker_gene_order gene arrangement
#' @param cluster_marker the marker gene information
#' @param normalize whether to normalize the data and how
#' @param Upper_Bound the upper bound of variance estimator of gene_error
#' @param Lower_Bound the lower bound of variance estimator of gene_error
#' @param weight_fpca the weight of functional pca



fun_spotlight2=function(geno,
                        ST.matrix,
                        cell.type.factor,
                        min_count=0.09,
                        ct_mode="median",
                        nkonts = 250,
                        norder = 3,
                        nfpca = 50,
                        basis_fun="bspline",
                        marker_gene_order="0",
                        cluster_marker=NULL,
                        normalize="0",
                        Upper_Bound=1,
                        Lower_Bound=1,
                        weight_gene='0',
                        weight_fpca='0'){








  #order sort
  if(marker_gene_order!=0){
    if(marker_gene_order=="2"){
      gene_mean=apply(geno, 1, median)
      gene_mean_sort=sort(gene_mean)
      geno=geno[names(gene_mean_sort),]
      ST.matrix = ST.matrix[rownames(geno), ]
    }
    if(marker_gene_order=="3"){
      cell.type.all.factor=unique(cell.type.factor)
      gene_sort=c()
      cluster_marker1=cluster_marker[!duplicated(cluster_marker$gene),]
      for (i in 1:length(cell.type.all.factor)) {
        gene_cluster_marker_name= unique(cluster_marker1$gene[cluster_marker1$cluster==cell.type.all.factor[i]])
        gene_cluster_marker_geno=geno[gene_cluster_marker_name,]
        gene_cluster_marker_mean=apply(gene_cluster_marker_geno, 1, median)
        if(i%%2==0){
          decrease=TRUE
        }else{
          decrease=TRUE
        }
        gene_cluster_marker_mean_sort=sort(gene_cluster_marker_mean,decreasing = decrease)
        gene_sort=c(gene_sort,names(gene_cluster_marker_mean_sort))
      }
      geno=geno[gene_sort,]
      ST.matrix = ST.matrix[rownames(geno), ]
    }
    if(marker_gene_order=="4"){
      cell.type.all.factor=unique(cell.type.factor)
      gene_sort=c()
      cluster_marker1=cluster_marker[!duplicated(cluster_marker$gene),]
      for (i in 1:length(cell.type.all.factor)) {
        gene_cluster_marker_name= cluster_marker1$gene[cluster_marker1$cluster==cell.type.all.factor[i]]
        gene_cluster_marker_geno=geno[gene_cluster_marker_name,cell.type.factor==cell.type.all.factor[i]]
        gene_cluster_marker_mean=apply(gene_cluster_marker_geno, 1, median)
        if(i%%2==0){
          decrease=TRUE
        }else{
          decrease=FALSE
        }
        gene_cluster_marker_mean_sort=sort(gene_cluster_marker_mean,decreasing = decrease)
        gene_sort=c(gene_sort,names(gene_cluster_marker_mean_sort))
      }
      geno=geno[gene_sort,]
      ST.matrix = ST.matrix[rownames(geno), ]
    }
    if(marker_gene_order=="5"){
      cell.type.all.factor=unique(cell.type.factor)
      gene_sort=c()
      gene_sort1=c()
      gene_sort2=c()
      cluster_marker1=cluster_marker[!duplicated(cluster_marker$gene),]
      for (i in 1:length(cell.type.all.factor)) {
        gene_cluster_marker_name= cluster_marker1$gene[cluster_marker1$cluster==cell.type.all.factor[i]]
        gene_cluster_marker_geno=geno[gene_cluster_marker_name,cell.type.factor==cell.type.all.factor[i]]
        gene_cluster_marker_mean=apply(gene_cluster_marker_geno, 1, mean)
        gene_cluster_marker_mean1=gene_cluster_marker_mean[gene_cluster_marker_mean>0]
        gene_cluster_marker_mean2=gene_cluster_marker_mean[gene_cluster_marker_mean=0]
        if(i%%2==0){
          decrease=TRUE
        }else{
          decrease=FALSE
        }
        gene_cluster_marker_mean_sort1=sort(gene_cluster_marker_mean1,decreasing = decrease)
        gene_cluster_marker_mean_sort2=sort(gene_cluster_marker_mean2,decreasing = decrease)
        gene_sort1=c(gene_sort1,names(gene_cluster_marker_mean_sort1))
        gene_sort2=c(gene_sort2,names(gene_cluster_marker_mean_sort2))
      }
      gene_sort=c(gene_sort1,gene_sort2)
      geno=geno[gene_sort,]
      ST.matrix = ST.matrix[rownames(geno), ]
    }}


  if(normalize!="0"){
    geno=geno[rowSums(geno)>0,]
    ST.matrix=ST.matrix[rownames(geno),]
    ST.matrix=ST.matrix[rowSums(ST.matrix)>0,]
    geno=geno[rownames(ST.matrix),]
  }


  start_time=Sys.time()
  ## base_function
  p = nrow(geno)
  t = seq.int(0,1,1/(p-1))
  nbasis = nkonts + norder - 2
  rangeval = c(0,1)
  if(basis_fun=="bspline"){
    genobasis <- fda::create.bspline.basis(rangeval= rangeval, nbasis=nbasis,  norder=norder)
  }
  if(basis_fun=="fourier"){
    genobasis = fda::create.fourier.basis(rangeval= rangeval, nbasis=nbasis)
  }
  if(basis_fun=="exponential"){
    genobasis=fda::create.exponential.basis(rangeval = rangeval,nbasis=nbasis)
  }
  if(basis_fun=="power"){
    genobasis=fda::create.power.basis(rangeval = rangeval,nbasis=nbasis)
  }
  if(basis_fun=="polygonal"){
    genobasis=fda::create.polygonal.basis(rangeval = rangeval)
  }
  if(basis_fun=="monomial"){
    genobasis=fda::create.monomial.basis(rangeval = rangeval,nbasis = nbasis)
  }
  Phi <- fda::eval.basis(t, genobasis)

  ## cell.type
  cell.type.factor.all=unique(cell.type.factor)
  table(cell.type.factor)

  if(normalize=="2"){
    geno=t(geno)
    geno=(geno/rowSums(geno))
    geno=t(geno)
    ST.matrix=t(ST.matrix)
    ST.matrix=ST.matrix/rowSums(ST.matrix)
    ST.matrix=t(ST.matrix)
  }
  if(normalize=="6"){
    geno=(geno-apply(geno,1,mean))/sqrt(apply(geno,1,var))
    ST.matrix=(ST.matrix-apply(ST.matrix,1,mean))/sqrt(apply(ST.matrix,1,var))
  }
  ##fpca
  if(nfpca!=0){
    if(weight_gene=='0'){
    #variance_geno_feature=apply(geno, 1, var)
    #W_vec=1/variance_geno_feature
    #W_vec[which(W_vec>Upper_Bound)]=Upper_Bound
    #W_vec[which(W_vec<Lower_Bound)]=Lower_Bound
    #sigma_e=diag(variance_geno_feature)
      cell_types=unique(cell.type.factor)
      var_pre=matrix(rep(0,length(cell_types)*nrow(geno)),nrow = nrow(geno))

      for (i in 1:nrow(var_pre)) {
        for (j in 1:ncol(var_pre)) {
          var_pre[i,j]=var(geno[i,which(cell.type.factor==cell_types[j])])
        }
      }
      W_vec=apply(var_pre,1,mean)
      W_vec[which(W_vec==0)]=1
      W_vec=1/W_vec
      #above_bound=10
      #down_bound=0.01
      W_vec[which(W_vec>Upper_Bound)]=Upper_Bound
      W_vec[which(W_vec<Lower_Bound)]=Lower_Bound
      W=diag(W_vec)
    #W=diag(W_vec)
    }else{
      W=diag(weight_gene)
    }
    #W=ginv(sigma_e)

    new_basia_coef=fpca_single_cell2(geno=geno,
                                     t=t,
                                     nbasis = nbasis,
                                     genobasis = genobasis,
                                     cell.type.factor = cell.type.factor,
                                     nfpca = nfpca,
                                     W=W)
    Phi=new_basia_coef
  }
  #Phi ngene*nbasis

  ## transform
  geno.max = max(geno)
  geno = geno / geno.max
  #Phi的维度为特征数量（基因数量）*基函数数�?
  #geno的维度为样本�?*特征数量（基因数量）
  ST.matrix.max = max(ST.matrix)
  ST.matrix= ST.matrix/ST.matrix.max


  if(normalize!=0){
    if(normalize=="1"){
    geno=(geno/rowSums(geno))
    ST.matrix=ST.matrix/rowSums(ST.matrix)
  }
    if(normalize=="4"){
      geno=(geno-apply(geno,1,mean))/sqrt(apply(geno,1,var))
      ST.matrix=(ST.matrix-apply(ST.matrix,1,mean))/sqrt(apply(ST.matrix,1,var))
    }
    if(normalize=="5"){
      geno=t(geno)
      geno=(geno-apply(geno,1,mean))/sqrt(apply(geno,1,var))
      geno=t(geno)
      ST.matrix=t(ST.matrix)
      ST.matrix=(ST.matrix-apply(ST.matrix,1,mean))/sqrt(apply(ST.matrix,1,var))
      ST.matrix=t(ST.matrix)
    }
    if(normalize=="3"){
      geno=t(geno)
      geno=(geno/rowSums(geno))
      geno=t(geno)
      ST.matrix=t(ST.matrix)
      ST.matrix=ST.matrix/rowSums(ST.matrix)
      ST.matrix=t(ST.matrix)
    }}


  ## funtional data
  U=MASS::ginv(t(Phi) %*% W %*% Phi) %*% t(Phi) %*% W %*% geno
  ST.U=MASS::ginv(t(Phi) %*% W %*% Phi) %*% t(Phi) %*% W %*% ST.matrix

  U=t(U)
  ST.U=t(ST.U)

  ## deconvolution data
  ct_topic_matrix=ct_topic_fun(U,cell.type.factor,ct_mode)
  if(weight_fpca=='0'){
  decon_mtrx=mixture_deconvolution_fun(ST.U,
                                       ct_topic_matrix,
                                       cell.type.factor,
                                       min_count=min_count)
  }else{
    decon_mtrx=mixture_deconvolution_fun2(ST.U,
                                         ct_topic_matrix,
                                         cell.type.factor,
                                         min_count=min_count,
                                         p_value_vector_weight=weight_fpca)
                                       }
  total_t <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  print(sprintf("the total time was %smins",total_t))
  return(decon_mtrx)
}
