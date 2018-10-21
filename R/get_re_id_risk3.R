#####use sdcMicro to assess re-identification risk####
rm(list=ls())
gc()

#turn off scientific notation
options(scipen=999)

#--global parameters
#https://www.hcup-us.ahrq.gov/reports/statbriefs/sb225-Inpatient-US-Stays-Trends.jsp
# wt<-sample(seq_len(round(pop_N/rsp)),nrow(.),replace=T)
# wt<-round(1/0.1*3)
# wt<-round(35000000*0.1*10/70000) #hospitalization counts/year * AKI incidence rate * 10 years/sample size
# wt<-1  #https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/freqCalc.r#L97

######### prepare #######
#--require libraries
source("./R/util.R")
source("./R/eval_re_id_risk.R")
require_libraries(c("tidyr",
                    "dplyr",
                    "magrittr",
                    "sdcMicro",
                    "h2o",
                    "scales"))

#--initialize h2o cluster
h2o.init(nthreads=-1)

#--load data
i<-2 #...done
dat<-read.csv(paste0("./data/Perspective_",i,".csv"))
out<-sel_safe_col(dat,
                  col_incl=c("X0","X1","X2"),
                  incl_inc=20,
                  uni_thresh=0.0001,
                  wt=round(1/0.1*3))

saveRDS(out,file=paste0("./output/ReID_risk_persp",i,".rda"))

h2o.shutdown(prompt=F)

