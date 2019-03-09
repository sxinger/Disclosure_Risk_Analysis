#####use sdcMicro to assess re-identification risk####
eval_ReID_risk<-function(dat,ns=10,rsp=0.1,csp=0.8,wt,verb=T,keep_est=T){
  #--global info
  n<-nrow(dat)
  vars<-colnames(dat)
  wt<-wt/rsp
  
  #--initialization container for risk scores
  # risk1_vec<-c()
  # risk2_vec<-c()
  # risk3_vec<-c()
  risk4_vec<-c()
  risk5_vec<-c()
  risk6_vec<-c()
  # risk7_vec<-c()
  
  #--loop over resamples
  for(i in seq_len(ns)){
    start_i<-Sys.time()
    
    #--sample data
    if(rsp<1){
      subset<-dat %>% sample_n(round(n*rsp)) 
    }else{
      subset<-dat
    }
    
    if(csp<1){
      vars_sel<-sample(vars,round(csp*length(vars)))
      vars_sel<-vars[vars %in% vars_sel] #put it in corrected order
    }else{
      vars_sel<-vars
    }
    
    subset1<-subset %>%
      dplyr::select(vars_sel) %>%
      mutate(weights=wt) %>%
      group_by_(.dots=vars_sel) %>%
      dplyr::summarize(counts=n(),
                       weights=sum(weights)) %>%
      ungroup
    
    #zero-frequency patterns
    for(col in length(vars_sel):(length(vars_sel)-pmin(round(length(vars_sel)*0.2),length(vars_sel)-1))){
      subset2<-subset %>%
        dplyr::select(vars_sel) %>%
        dplyr::mutate(temp=abs(get(vars_sel[col])-1)) %>%
        dplyr::select(-one_of(vars_sel[col])) %>%
        dplyr::rename_at(vars(starts_with("temp")),
                         funs(str_replace(.,"temp",vars_sel[col]))) %>%
        anti_join(subset %>%
                    dplyr::select(vars_sel),by=vars_sel) %>%
        dplyr::mutate(counts=0,
                      weights=0)
      
      subset1 %<>%
        bind_rows(subset2)
    }
    
    #https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/modRisk.R#L198
    subset1 %<>%
      mutate(EC=counts/weights) %>%
      mutate(EC=ifelse(is.nan(EC),0,EC)) %>%
      mutate(offset=log(EC+0.1))
    
    subset1_h2o<-as.h2o(as.matrix(subset1))
    
    if(verb){
      cat("...data aggregation done.\n")
    }
    
    #--inidivial risk
    subset2<-subset1 %>% dplyr::select(-counts,-EC,-offset)
    subset2<-as.data.frame(subset2) #tibble object doesn't work in measure_risk()
    rk_ms<-measure_risk(subset2,
                        keyVars=vars_sel,
                        w="weights")
    print(rk_ms) # report
    Fk<-rk_ms$Res[,3]
    indiv_rk2<-rk_ms$Res[,1]
    
    rk_bd<-2*(median(indiv_rk2,na.rm=T)+2*stats::mad(indiv_rk2,na.rm=T)) 
    rm(rk_ms);gc()
    cat("...inidividual risk measured.\n")
    
    #--model based global risk (possion)
    predb_idx<-which(!(colnames(subset1) %in% c("counts","EC","weights","offset")))
    target_idx<-which(colnames(subset1) == "counts")
    # x<-as.matrix(subset1 %>% dplyr::select(-weights,-counts))
    # y<-as.matrix(subset1 %>% dplyr::select(counts))
    # EC<-subset1$counts/subset1$weights  #offset term
    # offset<-log(EC + 0.1)
    # form<-as.formula(paste(c("counts", c(as.character(formulaM)), "+ offset(EC)"),
    #                        collapse=""))
    # pmod<-cv.glmnet(x,y,offset=offset,family="poisson")
    # fit_cnt<-predict(pmod,newx=x,newoffset=EC,s="lambda.min",type="response")
    pmod<-h2o.glm(x=predb_idx,
                  y=target_idx,
                  offset_column="offset",
                  training_frame=subset1_h2o,
                  family="poisson",
                  # link="log",
                  solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                  alpha=1,                       #penalized regression
                  nfolds=3,
                  lambda_search=TRUE,
                  early_stopping = TRUE,
                  ignore_const_cols= TRUE,
                  remove_collinear_columns=TRUE
                  # keep_cross_validation_predictions =T
    )
    h2o_fit<-predict(pmod,newdata=subset1_h2o)
    fit_cnt<-subset1 %>% 
      dplyr::select(-weights,-offset) %>%
      unite("pattern_str",vars_sel,sep="") %>%
      dplyr::select(pattern_str,EC,counts) %>%
      dplyr::rename(sample_freq=counts) %>%
      mutate(est_freq1=Fk,
             est_freq2=as.data.frame(h2o_fit)$predict)
    
    if(verb){
      cat("...log-linear model for global risk done.\n")
    }
    h2o.removeAll()
    
    #--marketer attacker model
    # mmod<-microaggregation(subset,
    #                        method="mdav",
    #                        aggr=20) #for better speed
    # # mmod<-mafast(subset,variables=vars,aggr=10) #faster performance--not always result in meaningful clusters
    # if(verb){
    #   cat("...nearest neightbor model for global risk done.\n")
    # }
    
    #--calculate risks
    # https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/modRisk.R#L248
    # indiv_rk3<-rescale(fit_cnt$est_freq2,c(1,max(fit_cnt$est_freq1)))
    indiv_rk3<-fit_cnt$est_freq2
    r1<-risk1(indiv_rk3, subset1$EC)/n 
    r2<-risk2(indiv_rk3, subset1$EC)/n 
    
    rk1<-round(sum(as.numeric(fit_cnt$sample_freq == 1))/n,4)                                             #disclosure risk model1 - sample uniqueness
    # rk2<-ifelse(nrow(fit_cnt[(fit_cnt$sample_freq==1),])==0,0,
    #             round(mean(as.numeric(fit_cnt[(fit_cnt$sample_freq==1),]$est_freq1==1)),4))               #disclosure risk model2 - population uniquness (NB)
    # rk3<-round(sum(as.numeric(fit_cnt$sample_freq == 1) * indiv_rk2)/n,4)                                 #disclosure risk model3 - matching rate (NB)
    rk4<-round(sum(as.numeric(fit_cnt$sample_freq == 1) * r1),4)                                          #disclosure risk model4 - population uniquness (Poisson)
    rk5<-round(sum(as.numeric(fit_cnt$sample_freq == 1) * r2),4)                                          #disclosure risk model5 - matching rate (Poisson)
    rk6<-round(sum((indiv_rk2 >= max(0.1,rk_bd)),na.rm=T)/n,4)                                            #disclosure risk model6 - records at risk
    # rk7<-round(mean(max(indiv_rk2),max(r2)),4)                                                            #disclosure risk model7 - worst-case senario
    # rk6<-dRisk(obj=subset,xm=mmod$mx)                                                                   #disclosure risk model6 - adjusted uniquness
    cat("...risk scores calculation done.\n")
    
    #--attach results
    # risk1_vec<-c(risk1_vec,rk1)
    # risk2_vec<-c(risk2_vec,rk2)
    # risk3_vec<-c(risk3_vec,rk3)
    risk4_vec<-c(risk4_vec,rk4)
    risk5_vec<-c(risk5_vec,rk5)
    risk6_vec<-c(risk6_vec,rk6)
    # risk7_vec<-c(risk7_vec,rk7)
    
    lapse_i<-Sys.time()-start_i
    if(verb){
      cat("finish evaluate sample",i,"in",lapse_i,units(lapse_i),".\n")
    }
  }
  
  #--calculate confidence interval
  risk_summ<-data.frame(
                        # risk1=risk1_vec,
                        # risk2=risk2_vec,
                        # risk3=risk3_vec,
                        risk4=risk4_vec,
                        risk5=risk5_vec,
                        risk6=risk6_vec
                        # risk7=risk7_vec
                        )
  
  if(keep_est & ns==1){
    risk_df<-list(est_freq=fit_cnt,
                  risk_summ=risk_summ)
  }else{
    risk_df<-list(risk_summ=risk_summ)
  }
  
  return(risk_df)
}



sel_safe_col<-function(dat,wt,col_incl=c("X0","X1","X2"),
                       ns=10,rsp=0.1,csp=0.8,
                       incl_inc=20,uni_thresh=0.02,
                       sparse_rank=data.frame()){
  col_add<-colnames(dat)[colnames(dat) %in% col_incl]
  dat_i<-dat %>% dplyr::select(col_add)
  out<-eval_ReID_risk(dat=dat_i,ns=ns,rsp=rsp,csp=csp,wt = wt)
  chk_pt<-chk_pt_old<-mean(out$risk_summ$risk4)
  
  #--iteration
  remain_col<-nrow(sparse_rank)
  col_add_try<-col_add
  incl_until<-0
  risk_inc<-c()
  rewind<-F
  while(incl_inc>0 && remain_col>0){
    while(chk_pt<=uni_thresh && remain_col>0){
      chk_pt_old<-chk_pt
      col_add<-col_add_try
      if(!rewind){
        risk_inc %<>%
          bind_rows(out$risk_summ %>% 
                      summarise_all(.,mean,na.rm=T) %>%
                      mutate(num_vars=length(unique(col_add)),
                             vars_lst=paste(c(col_incl,
                                              sparse_rank$vars[1:min(length(sparse_rank$vars),incl_until)]),collapse=","))) %>%
          unique
      }
      
      incl_until<-incl_until+incl_inc
      col_add_try<-c(col_add_try,sparse_rank$vars[(incl_until-incl_inc+1):min(length(sparse_rank$vars),incl_until)])
      dat_i<-dat %>% dplyr::select(col_add_try)
      out<-eval_ReID_risk(dat_i,ns=ns,rsp=rsp,csp=csp,wt = wt)
      chk_pt<-mean(out$risk_summ$risk4)
      
      remain_col<-nrow(sparse_rank)-incl_until
      rewind<-F
    }
    #rewind
    incl_until<-incl_until-incl_inc
    col_add_try<-col_add
    incl_inc<-floor(incl_inc/2)
    chk_pt<-chk_pt_old
    rewind<-T
  }
  risk_sel<-list(col_sel=col_add,
                 risk_inc=risk_inc)
  
  return(risk_sel)
}