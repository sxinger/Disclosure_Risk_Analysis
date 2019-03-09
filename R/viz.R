#####use sdcMicro to assess re-identification risk####
rm(list=ls())
gc()

#turn off scientific notation
options(scipen=999)

######### prepare #######
#--require libraries
source("./R/util.R")
require_libraries(c("tidyr",
                    "dplyr",
                    "magrittr",
                    "ggplot2",
                    "ggrepel"))

#load data
perf_summ<-readRDS("./output/perf_summary.rda")
calib_equal_bin<-readRDS("./output/perf_calib.rda")
var_sp<-read.csv("./output/sparsity_rank.csv",stringsAsFactors = F) %>%
  left_join(data.frame(type="demo",vars=paste0("X",1:3-1),stringsAsFactors = F) %>%
              bind_rows(data.frame(type="vital",vars=paste0("X",4:8-1),stringsAsFactors = F)) %>%
              bind_rows(data.frame(type="lab",vars=paste0("X",9:22-1),stringsAsFactors = F)) %>%
              bind_rows(data.frame(type="drg",vars=paste0("X",23:337-1),stringsAsFactors = F)) %>%
              bind_rows(data.frame(type="comorb",vars=paste0("X",338:366-1),stringsAsFactors = F)) %>%
              bind_rows(data.frame(type="med",vars=paste0("X",367:1637-1),stringsAsFactors = F)) %>%
              bind_rows(data.frame(type="ccs",vars=paste0("X",1638:1917-1),stringsAsFactors = F)),
            by="vars") %>%
  left_join(data.frame(type=c("demo",
                              "vital",
                              "lab",
                              "drg",
                              "comorb",
                              "med",
                              "ccs"),
                       log10_cat_cnt=c(log10(6*4*2),
                                       log10(5*5*5*6*6),
                                       log10(3^length(paste0("X",9:22-1))),
                                       log10(2^length(paste0("X",23:337-1))),
                                       log10(2^length(paste0("X",338:366-1))),
                                       log10(2^length(paste0("X",367:1637-1))),
                                       log10(2^length(paste0("X",1638:1917-1)))),
                       stringsAsFactors = F),
            by="type") %>%
  dplyr::filter(!is.na(type)) %>%
  dplyr::mutate(type=recode(type,
                            "demo"="1.Demographic",
                            "vital"="2.Vital_Sign",
                            "lab"="3.Labs",
                            "comorb"="4.Comorbidity",
                            "ccs"="5.Diagnosis(CCS)",
                            "med"="6.Medication",
                            "drg"="7.Admission DRG"))


ggplot(var_sp,aes(x=type,y=sparsity))+
  geom_boxplot()+
  geom_text(data=var_sp %>% dplyr::select(type,sparsity,log10_cat_cnt) %>% 
                 group_by(type,log10_cat_cnt) %>% 
                 dplyr::summarize(sparsity=median(sparsity)) %>%
                 ungroup,
               aes(x=type,y=sparsity,label=paste0("10^",round(log10_cat_cnt))),nudge_y = 0.05)+
  labs(x="Data Type",y="Sparsity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#de-id risk
risk<-perf_summ %>%
  dplyr::select(num_var,risk2,risk3,risk4,risk5,risk6,risk7) %>%
  gather(risk,score,-num_var)

ggplot(risk,aes(x=num_var,y=score))+
  geom_point()+
  facet_wrap(~risk)


# overall performance vs risk
ss_risk<-perf_summ %>%
  dplyr::select(num_var,risk3,risk5,fold,overall_meas,meas_val) %>%
  filter(overall_meas %in% c("roauc","opt_sens","opt_spec",
                             "prauc1","opt_prec","opt_rec")) %>%
  mutate(overall_meas=recode(overall_meas,
                             "roauc" = "1.AUROC",
                             "opt_sens" = "2.Sensitivity",
                             "opt_spec" = "3.Specificity",
                             "prauc1" = "4.AUPRC",
                             "opt_prec" = "5.Precision",
                             "opt_rec" = "6.Recall")) %>%
  gather(risk,score,-fold,-num_var,-overall_meas,-meas_val) %>%
  group_by(num_var,risk) %>%
  dplyr::mutate(risk_score=mean(score)) %>%
  ungroup %>%
  group_by(num_var,risk,risk_score,overall_meas) %>%
  dplyr::summarize(score=mean(meas_val),
                   score_low=min(meas_val),
                   score_up=max(meas_val)) %>%
  ungroup %>% unique %>%
  dplyr::filter(num_var<230)

ggplot(ss_risk,aes(x=risk_score,y=score,label=num_var))+
  geom_point(aes(size=num_var),alpha=0.5)+
  geom_text_repel()+
  # geom_line()+
  geom_smooth(method = 'loess',formula='y ~ x')+
  # geom_errorbar(aes(ymin=roauc_low,
  #                   ymax=roauc_up))+
  labs(x="Re-id Risk Score",y="Performance Metric")+
  guides(size=FALSE)+
  facet_wrap(~overall_meas,scales="free")

