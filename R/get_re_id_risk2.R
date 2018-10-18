#####use sdcMicro to assess re-identification risk####
rm(list=ls())
gc()

#--global parameters
#https://www.hcup-us.ahrq.gov/reports/statbriefs/sb225-Inpatient-US-Stays-Trends.jsp
AKI_US<-35000000*0.1*10
#hospitalization counts/year * AKI incidence rate * 10 years


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
# dat %<>%
#   mutate_at(vars(starts_with("X")),funs(factor)) #group conversion to factor

#--remove sparse rows
# sparse_bd<-0.995
# dat %<>% 
#   mutate(zero_infl=rowSums((.==0))/length(vars)) %>%
#   filter(zero_infl < sparse_bd) %>%
#   dplyr::select(-zero_infl)

## key variable selection - demo
# col_sel<-colnames(dat)[colnames(dat) %in% paste0("X",1:3)]
# dat %<>% dplyr::select(col_sel) #3
# out<-eval_ReID_risk(dat,ns=1,rsp=1,csp=1)
# pass

data_type<-data.frame(type=c("demo",
                             "vital",
                             "lab",
                             "drg",
                             "comorb",
                             "med",
                             "ccs"),
                      cat_cnt=c(6*4*2,
                                5*5*5*6*6,
                                3^length(paste0("X",9:22-1)),
                                2^length(paste0("X",23:337-1)),
                                2^length(paste0("X",338:366-1)),
                                2^length(paste0("X",367:1637-1)),
                                2^length(paste0("X",1638:1917-1))),
                      stringsAsFactors = F) %>%
  arrange(cat_cnt)


idx_lst<-list(demo=paste0("X",1:3-1),
              vital=paste0("X",4:8-1),
              lab=paste0("X",9:22-1),
              drg=paste0("X",23:337-1),
              comorb=paste0("X",338:366-1),
              med=paste0("X",367:1637-1),
              ccs=paste0("X",1638:1917-1))

risk_ind<-c()
for(i in 1:nrow(data_type)){
  col_sel<-unique(c(idx_lst[["demo"]],
                    colnames(dat)[colnames(dat) %in% idx_lst[[data_type$type[i]]]]))
  dat_i<<-dat %>% dplyr::select(col_sel)
  out<-eval_ReID_risk(dat_i,ns=1,rsp=1,csp=1)
  
  risk_ind %<>%
    bind_rows(out$risk_summ %>% mutate(data_type_incl=data_type$type[i]))
}
saveRDS(risk_ind,file="./output/ReID_risk_var_ind.rda")



# search for minimal set of columns that still secure enough to release
risk_order<-seq_len(nrow(data_type))[order(risk_ind$risk3)]
risk_inc<-c()
col_add<-c()
chk_pt<-0
for(i in 1:nrow(data_type)){
  #--copy the last check point value
  chk_pt_last<-chk_pt
  
  #--add the whole data type, tentatively
  col_add_try<-c(col_add,colnames(dat)[colnames(dat) %in% idx_lst[[data_type$type[risk_order[i]]]]])
  dat_i<<-dat %>% dplyr::select(col_add_try)
  out<-eval_ReID_risk(dat_i,ns=1,rsp=1,csp=1)
  chk_pt<-max(unlist(out$risk_summ))
  
  #--inner loop for screening out rare features
  col_ii<-colnames(dat)[colnames(dat) %in% idx_lst[[data_type$type[risk_order[i]]]]]
  dat_ii<-dat_i %>% 
    dplyr::select(col_ii) %>% 
    group_by_(col_ii) %>%
    summarise_all(function(x) mean(x!=0)) %>%
    ungroup %>%
    gather(vars,sparsity) %>%
    arrange(desc(sparsity))
  
  chk_pt<-chk_pt_last
  remain_col<-nrow(dat_ii)
  col_add_try<-col_add
  incl_until<-1
  while(chk_pt<=0.03 | remain_col>0){
    col_add_try<-c(col_add_try,dat_ii$vars[1:incl_until])
    dat_i<<-dat %>% dplyr::select(col_add_try)
    out<-eval_ReID_risk(dat_i,ns=1,rsp=1,csp=1)
    chk_pt<-max(unlist(out$risk_summ))
    remain_col<-nrow(dat_ii)-incl_until
  }
  
  col_add<-col_add_try
  risk_inc %<>%
    bind_rows(out$risk_summ %>% 
                mutate(num_data_type=i,
                       num_vars=length(col_add)))
}

risk_sel<-list(col_sel=col_add,
               risk_inc=risk_rinc)

saveRDS(risk_sel,file="./output/ReID_risk_var_inc.rda")
