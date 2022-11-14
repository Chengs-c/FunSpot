#' @param ST.U the funcitonal data matrix of ST data
#' @param ct_topic_fun the topic matrix
#' @param cell.type.factor the annotation of scRNA data
#' @param min_count the minmum weight contribution
#' @param p_value_vector_weight the weight of fpc

mixture_deconvolution_fun2=function(ST.U,ct_topic_fun,cell.type.factor,min_count,p_value_vector_weight){

  ST.U=t(ST.U)
  cell.type.factor.all=unique(cell.type.factor)
  decon_mtrx=matrix(nrow = length(cell.type.factor.all),
                    ncol = ncol(ST.U),
                    0)
  rownames(decon_mtrx)=colnames(ct_topic_fun)

  ct_topic_fun1=diag(sqrt(p_value_vector_weight)) %*% ct_topic_fun
  ST.U1=diag(sqrt(p_value_vector_weight))%*%ST.U
  #print("Deconvoluting spot to")
  #total=ncol(profile_mtrx)
  #pb=txtProgressBar(min=0,max=total,style = 3)

  for (i in 1:ncol(ST.U)) {
    ##nnls?õ?ϸ??????????
    nnls_pred=nnls::nnls(A=ct_topic_fun1,b=ST.U1[,i])
    weights=nnls_pred$x

    #?õ?????
    comp=weights/sum(weights)

    #ȥ????С
    comp[comp<min_count]=0
    weights[comp<min_count]=0

    comp_prop=comp/sum(comp)
    comp_prop[is.na(comp_prop)]=0

    decon_mtrx[,i]=comp_prop
    #setTxtProgressBar(pb,i)

  }

  #close(pb)

  return(decon_mtrx)
}
