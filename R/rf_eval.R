#####use sdcMicro to assess re-identification risk####
rm(list=ls())
gc()

#turn off scientific notation
options(scipen=999)

######### prepare #######
#--require libraries
source("./R/util.R")
source("./R/eval_re_id_risk.R")
require_libraries(c("tidyr",
                    "dplyr",
                    "magrittr",
                    "stringr",
                    "sdcMicro",
                    "h2o",
                    "scales"))

#--load data
i<-1
dat<-read.csv(paste0("./data/Perspective_",i,".csv"))
dat$X1917<-as.factor(dat$X1917)

#sparsity rank matrix
sparse_rank<-dat %>% 
  dplyr::select(-X1917) %>%
  summarise_all(function(x) mean(as.numeric(x!=0))) %>%
  gather(vars,sparsity) %>%
  arrange(desc(sparsity))

#partition index
rsample_idx<-data.frame(row_num=1:nrow(dat),
                        part=sample(c("T","V"),nrow(dat),c(0.8,0.2),replace=T))

#hyper-parameters
hyper_params<-list(ntrees = c(1,5,seq(10,100,10)), 
                   max_depth = seq(1,20), 
                   min_rows = c(1,5,10,20,50,100),
                   sample_rate = seq(0.3,1,0.1),
                   col_sample_rate_per_tree = seq(0.3,1,0.1)
                   )

search_criteria<-list(strategy = "RandomDiscrete", 
                      max_runtime_secs = 200, 
                      max_models = 20, 
                      stopping_metric = "AUC", 
                      stopping_tolerance = 0.0001, 
                      stopping_rounds = 5)


#--initialize h2o cluster
h2o.init(nthreads=-1)

#--start experiment
rslt_lst<-list()
start_v<-Sys.time()

train<-dat[rsample_idx$part=="T",]
test<-dat[rsample_idx$part=="V",]
# rm(dat);gc()

#scrub columns
start_vv<-Sys.time()
out<-eval_ReID_risk(dat[,c("X0","X1","X2")],
                    ns=5,
                    rsp=0.1,
                    csp=1,
                    wt=1)

#baseline with only demo
start_vv<-Sys.time()
tr_h2o<-as.h2o(train[,unique(c("X1917","X0","X1","X2"))])
ts_h2o<-as.h2o(test[,unique(c("X1917","X0","X1","X2"))])

rf_j<-h2o.grid(algorithm = "drf",
               grid_id="rf_grid",
               x=2:4,
               y=1,
               training_frame=tr_h2o,
               nfolds=3,
               hyper_params=hyper_params,
               search_criteria = search_criteria)

rf_sorted_grid<-h2o.getGrid(grid_id="rf_grid",sort_by = "AUC",decreasing = T)
rf_bst<-h2o.getModel(rf_sorted_grid@model_ids[[1]])
ts_pred<-h2o.predict(rf_bst,newdata=ts_h2o)

out_risk<-data.frame(num_var=3,
                     var_lst="X0,X1,X2",
                     out$risk_summ %>% summarise_all(.,mean,na.rm=T),
                     opt_model_spec=paste(rf_sorted_grid@summary_table[1,1:4],collapse = ","))

out_eval<-data.frame(num_var=3,
                     pred=as.data.frame(ts_pred)$p1,
                     real=as.numeric(as.character(as.vector(ts_h2o[,1]))))

rslt<-list(risk=out_risk,eval=out_eval)
rslt_lst[[1]]<-rslt

lapse_vv<-Sys.time()-start_vv
cat("...finish initial risk evaluation in",lapse_vv,units(lapse_vv),".\n")


#scrub columns
start_vv<-Sys.time()

out<-sel_safe_col(dat[,which(colnames(dat)!="X1917")],
                  wt=1,
                  col_incl=c("X0","X1","X2"),
                  ns=5,
                  rsp=0.1,
                  csp=1,
                  incl_inc=10,
                  uni_thresh=0.3,
                  sparse_rank=sparse_rank)

lapse_vv<-Sys.time()-start_vv
cat("...finish iterative risk evaluation in",lapse_vv,units(lapse_vv),".\n")

#iterative rf evaluation
num_var<-out$risk_inc$num_vars
for(j in seq_along(num_var)){
  start_vj<-Sys.time()
  
  h2o.removeAll()
  var_j<-unlist(strsplit(out$risk_inc$vars_lst,",")[[1]])
  tr_h2o<-as.h2o(train[,unique(c("X1917",var_j))])
  ts_h2o<-as.h2o(test[,unique(c("X1917",var_j))])
  
  rf_j<-h2o.grid(algorithm = "drf",
                 grid_id="rf_grid",
                 x=2:length(var_j),
                 y=1,
                 training_frame=tr_h2o,
                 nfolds=5,
                 hyper_params=hyper_params,
                 search_criteria = search_criteria)
  
  rf_sorted_grid<-h2o.getGrid(grid_id="rf_grid",sort_by = "AUC",decreasing = T)
  rf_bst<-h2o.getModel(rf_sorted_grid@model_ids[[1]])
  ts_pred<-h2o.predict(rf_bst,newdata=ts_h2o)
  
  out_risk<-data.frame(num_var=out$risk_inc$num_var[j],
                       var_lst=out$risk_inc$vars_lst[j],
                       risk1=out$risk_inc$risk1[j],
                       risk2=out$risk_inc$risk2[j],
                       risk3=out$risk_inc$risk3[j],
                       risk4=out$risk_inc$risk4[j],
                       risk5=out$risk_inc$risk5[j],
                       risk6=out$risk_inc$risk6[j],
                       risk7=out$risk_inc$risk7[j],
                       opt_model_spec=paste(rf_sorted_grid@summary_table[1,1:4],collapse = ","))
  
  out_eval<-data.frame(num_var=out$risk_inc$num_var[j],
                       pred=as.data.frame(ts_pred)$p1,
                       real=as.numeric(as.character(as.vector(ts_h2o[,1]))))
  
  rslt<-list(risk=out_risk,eval=out_eval)
  rslt_lst[[(j+1)]]<-rslt
  
  lapse_vj<-Sys.time()-start_vj
  bm_v<-c(bm_v,lapse_vj)
  bm_nm_v<-c(bm_nm_v,paste0("include_var_",j))
  cat("...finish risk evaluation for",out$risk_inc$num_var[j],"variables in",lapse_vj,units(lapse_vj),".\n")
}

lapse_v<-Sys.time()-start_v
bm_v<-c(bm_v,lapse_v)
bm_nm_v<-c(bm_nm_v,"overall")

saveRDS(rslt_lst,file="./output/deid_vs_pred_results2.rda")


#shut down h2o instance
h2o.shutdown(prompt=F)


