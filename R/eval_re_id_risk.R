#####use sdcMicro to assess re-identification risk####
eval_ReID_risk<-function(dat,ns=10,rsp=0.6,csp=0.8, verb=T){
  #--global info
  n<-nrow(dat)
  vars<-colnames(dat)
  
  #--initialization
  risk1_vec<-c()
  risk2_vec<-c()
  risk3_vec<-c()
  risk4_vec<-c()
  risk5_vec<-c()
  risk6_vec<-c()
  
  #--loop over resamples
  for(i in seq_len(ns)){
    start_i<-Sys.time()
    
    #--sample data
    subset<-dat %>% sample_n(round(n*rsp)) 
    vars_sel<-sample(vars,round(csp*length(vars)))
    
    subset1<-subset %>%
      dplyr::select(vars_sel) %>%
      mutate(weights=round(n*rsp)) %>%
      group_by_(.dots=vars_sel) %>%
      mutate(counts=n(),
             weights=sum(weights)) %>%
      ungroup
    
    #https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/modRisk.R#L198
    subset1 %<>%
      mutate(EC=counts/weights) %>%
      mutate(offset=log(EC + 0.1))

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
    indiv_rk<-rk_ms$Res[,1]
    rk_bd<-2*(median(indiv_rk)+2*stats::mad(indiv_rk)) 
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
                  solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                  alpha=0,                       #ridge regression
                  nfolds=5,
                  lambda_search=TRUE,
                  early_stopping = TRUE,
                  ignore_const_cols= TRUE,
                  remove_collinear_columns=TRUE,
                  keep_cross_validation_predictions =T
    )
    h2o_fit<-h2o.getFrame(pmod@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
    fit_cnt<-as.data.frame(h2o_fit)$predict
    if(verb){
      cat("...procecutor attacker model done.\n")
    }
    h2o.removeAll()
    
    #--marketer attacker model
    mmod<-microaggregation(subset,
                           method="mdav",
                           aggr=20) #for better speed
    # mmod<-mafast(subset,variables=vars,aggr=10) #faster performance--not always result in meaningful clusters
    if(verb){
      cat("...marketer attacker model done.\n")
    }
    
    #--calculate risks
    # https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/modRisk.R#L248
    r1<-risk1(fit_cnt, subset1$EC)/n # 1. estimates the number of sample uniques that are population unique
    r2<-risk2(fit_cnt, subset1$EC)/n # 2. estimates the number of correct matches of sample uniques
    
    rk1<-sum(as.numeric(subset1$counts == 1) * r1) #disclosure risk model1 - uniquness
    rk2<-dRisk(obj=subset,xm=mmod$mx)              #disclosure risk model3 - uniquness adjusted for continuous variable
    rk3<-sum(as.numeric(subset1$counts == 1) * r2) #disclosure risk model2 - success rate
    rk4<-sum((indiv_rk >= max(0.1,rk_bd)))         #disclosure risk model5 - records at risk
    rk5<-sum(indiv_rk)                             #disclosure risk model4 - global risk
    rk6<-max(indiv_rk)                             #disclosure risk model6 - worst-case senario
    cat("...risk score calculation done.\n")
    
    #--attach results
    risk1_vec<-c(risk1_vec,rk1)
    risk2_vec<-c(risk2_vec,rk2)
    risk3_vec<-c(risk3_vec,rk3)
    risk4_vec<-c(risk4_vec,rk4)
    risk5_vec<-c(risk5_vec,rk5)
    risk6_vec<-c(risk6_vec,rk6)
    
    lapse_i<-Sys.time()-start_i
    if(verb){
      cat("finish evaluate sample",i,"in",lapse_i,units(lapse_i),".\n")
    }
  }
  
  #--calculate confidence interval
  risk_df<-data.frame(risk1=risk1_vec,
                      risk2=risk2_vec,
                      risk3=risk3_vec,
                      risk4=risk4_vec,
                      risk5=risk5_vec,
                      risk6=risk6_vec)

  return(risk_df)
}