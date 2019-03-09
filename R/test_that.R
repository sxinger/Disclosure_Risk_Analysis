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

#--initialize h2o cluster
h2o.init(nthreads=-1)

#--test eval_ReID_risk()
out<-eval_ReID_risk(dat=dat[,1:100],
                    ns=5,
                    rsp=0.1,
                    csp=1,
                    wt=1)
out$risk_summ

#--test sel_safe_col()
#sparsity rank matrix
sparse_rank<-dat %>% 
  dplyr::select(-X1917) %>%
  summarise_all(function(x) mean(as.numeric(x!=0))) %>%
  gather(vars,sparsity) %>%
  arrange(desc(sparsity))

out<-sel_safe_col(dat=dat,
                  wt=1,
                  col_incl=c("X0","X1","X2"),
                  ns=5,
                  rsp=0.1,
                  csp=1,
                  incl_inc=10,
                  uni_thresh=0.05,
                  sparse_rank=sparse_rank)
out$col_sel
out$risk_inc
