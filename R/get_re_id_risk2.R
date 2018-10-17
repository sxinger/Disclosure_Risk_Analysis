#####use sdcMicro to assess re-identification risk####
rm(list=ls())
gc()

#--global parameters
#https://www.hcup-us.ahrq.gov/reports/statbriefs/sb225-Inpatient-US-Stays-Trends.jsp
AKI_US<-35000000*0.1*10


######### prepare #######
#--require libraries
source("./R/util.R")
source("./R/eval_re_id_risk.R")
require_libraries(c("tidyr",
                    "dplyr",
                    "magrittr",
                    "sdcMicro",
                    "h2o"))

#--initialize h2o cluster
h2o.init(nthreads=-1)

#--load data
dat<-read.csv("./data/Perspective_1/Perspective_1.csv")
# mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor

#--remove sparse rows
# sparse_bd<-0.995
# dat %<>% 
#   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
#   filter(zero_infl < sparse_bd) %>%
#   dplyr::select(-zero_infl)

## key variable selection - remove med
col_med_rm<-colnames(dat)[colnames(dat) %in% paste0("X",c(1:366,1638:1918))]
dat %<>% dplyr::select(col_med_rm) #610

out<-eval_ReID_risk(dat,ns=1,rsp=1,csp=1, verb=T)
saveRDS(out,file="./output/ReID_risk_persp1.rda")


######### perspective 2 ########
#--load data
dat<-read.csv("./data/Perspective_2/Perspective_2.csv")
# mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor

#--remove sparse rows
# sparse_bd<-0.995
# dat %<>% 
#   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
#   filter(zero_infl < sparse_bd) %>%
#   dplyr::select(-zero_infl)

out<-eval_ReID_risk(dat,ns=10,rsp=0.6,csp=0.8, verb=T)
saveRDS(out,file="./output/ReID_risk_persp2.rda")


######### perspective 3 ########
#--load data
# dat_suffix<-c(1,2,3,5,7,30)
dat_suffix<-1
for(i in seq_along(dat_suffix)){
  dat<-read.csv(paste0("./data/Perspective_3/CSV/Perspective_3_d0_p",
                       dat_suffix[i],".csv"))
  # mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor
  
  #--remove sparse rows
  # sparse_bd<-0.995
  # dat %<>% 
  #   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
  #   filter(zero_infl < sparse_bd) %>%
  #   dplyr::select(-zero_infl)

  out<-eval_ReID_risk(dat,ns=10,rsp=0.6,csp=0.8, verb=T)
  saveRDS(out,file=paste0("./output/ReID_risk_persp3_p",dat_suffix[i],".rda"))
}

######### perspective 4 ########
#--load data
# dat_suffix<-c(0,1,2,3,4)
dat_suffix<-1
for(i in seq_along(dat_suffix)){
  dat<-read.csv(paste0("./data/Perspective_4/Perspective_4_d",
                       dat_suffix[i],"_p1.csv"))
  # mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor
  
  #--remove sparse rows
  # sparse_bd<-0.995
  # dat %<>% 
  #   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
  #   filter(zero_infl < sparse_bd) %>%
  #   dplyr::select(-zero_infl)
  
  out<-eval_ReID_risk(dat,ns=10,rsp=0.6,csp=0.8, verb=T)
  saveRDS(out,file=paste0("./output/ReID_risk_persp4_d",dat_suffix[i],".rda"))
}


######### perspective 4a ########
#--load data
# dat_suffix<-c(0,1,2,3,4,5,6,7,8,9)
# for(i in seq_along(dat_suffix)){
#   dat<-read.csv(paste0("./data/Perspective_4a/",dat_suffix[i],
#                        "/l1_train",dat_suffix[i],".csv"))
#   # mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor
#   
#   #--remove sparse rows
#   # sparse_bd<-0.995
#   # dat %<>% 
#   #   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
#   #   filter(zero_infl < sparse_bd) %>%
#   #   dplyr::select(-zero_infl)
#   
#   out<-eval_ReID_risk(dat,ns=10,rsp=0.6,csp=0.8, verb=T)
#   saveRDS(out,file=paste0("./output/ReID_risk_persp4a_f",dat_suffix[i],".rda"))
# }


#--close h2o cluster
h2o.shutdown(prompt = FALSE)



######### failed attemps #############
# dat<-read.csv("./data/Perspective_1/Perspective_1.csv")
# n<-nrow(dat)
# vars<-colnames(dat)
# form<-as.formula(paste0("~",paste(vars,collapse = "+")))
# 
# dat %<>% mutate(sampling_wt=n)
# dat_sdc<-createSdcObj(test_dat,
#                       keyVars=vars,
#                       weightVar = "sampling_wt")
# #
# risk1<-modRisk(dat_sdc,
#                formulaM=form)
# Error in rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep) : 
# invalid 'times' value
#debug
# x<-test_dat
# grid<-data.table(expand.grid(lapply(1:20, function(t) unique((x[[vars[t]]])))))
#expanded grid can be too large
