# get independent variables
get_ind_vars <- function(df, include_race) {
  unique_cols <- apply(df, 2, function(x) length(unique(x)))

  ind_vars <- setdiff(
    c(
      "neuro_post", "sex", "age_group",
      "pre_admission_cns", "pre_admission_pns",
      paste0(".fittedPC", 1:10)
    ),
    names(unique_cols)[unique_cols == 1]
  )

  if (include_race) {
    ind_vars <- c(ind_vars, "race")
  }
  ind_vars
}

# 8.15.2022: Updated function from Xuan
# tcut = the censor time
run_coxregressions <- function(df, include_race = TRUE, tcut=60, blur_abs, mask_thres, is_pediatric = NULL) {

  if(is_pediatric = FALSE){
    df <- df %>% filter(age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus"))
  } else if(is_pediatric = TRUE){
    df <- df %>% filter(age_group %in% c("00to02", "06to11", "12to17"))
  }

  ind_vars <- get_ind_vars(df, include_race)

  deceased_reg_elix <-
    run_coxregression(df, 'deceased', ind_vars, tcut=tcut, blur_abs, mask_thres)

  time_first_discharge_reg_elix <-
    run_coxregression(df, "time_to_first_discharge", ind_vars, tcut=tcut, blur_abs, mask_thres)

  severe_reg_elix <-
    run_coxregression(df, 'severe', ind_vars, tcut=tcut, blur_abs, mask_thres)

  list(
    severe_reg_elix = severe_reg_elix,
    deceased_reg_elix = deceased_reg_elix,
    time_first_discharge_reg_elix= time_first_discharge_reg_elix
  )
}


run_coxregression <-function(df, depend_var, ind_vars, tcut=60, blur_abs, mask_thres) {
  if (length(unique(df[, depend_var, drop = T])) <= 1)
    return(NULL)
  independ_vars <- paste(ind_vars, collapse = ' + ')

  print(paste('Begin evaluating', depend_var,  'outcome'))

  if (depend_var=="deceased"){
    # censor time is the shortest time between days_since_admission or tcut
    df$c=pmin(df$days_since_admission,tcut)
    # if missing, temporarily set 1000 days
    df$time_to_death[is.na(df$time_to_death)]=1000
    # time is the shortest time between time_to_death or censor tme
    df$time=apply(cbind(df$time_to_death,df$c),1,min)
    # delta==1 if time_to_death occurs before censor time
    df$delta=as.numeric(df$time_to_death<=df$c)
    # remove those who meet outcome on day of admission
    df <- df %>% filter(!(time == 0 & delta == 1))

  }else if (depend_var=='time_to_first_discharge'){
    # censor time is the shortest time between days_since_admission or tcut
    df$c=pmin(df$days_since_admission,tcut)
    # if missing, temporarily set 1000 days
    df$time_to_first_discharge[is.na(df$time_to_first_discharge)]=1000
    # if missing, temporarily set 1000 days
    df$time_to_death[is.na(df$time_to_death)]=1000
    # time is the shortest time_to_first_discharge, time_to_death, censor time
    df$time=apply(cbind(df$time_to_first_discharge,df$time_to_death,df$c),1,min)
    # which.min will identify the column integer with the minimum time (e.g. if time_to_first_discharge is the smallest, delta = 1)
    df$delta=apply(cbind(df$time_to_first_discharge,df$time_to_death,df$c),1,which.min)
    # remove those who meet outcome on day of admission
    df <- df %>% filter(!(time == 0 & delta == 1))

  } else if (depend_var=='severe'){
    # if missing, temporarily set 1000 days
    df$time_to_severe[is.na(df$time_to_severe)]=1000
    df$time_to_neuro[is.na(df$time_to_neuro)]=1000
    # remove those who meet outcome on day of admission
    df <- df %>% filter(!time_to_severe == 0)
  }

  if (depend_var!='severe'){

    # mean adjusted survival time
    tryCatch(
      {

        covariate=model.matrix(as.formula(paste("Surv(time,delta==1)", '~',
                                                independ_vars)),data=df)[,-1] #[-1] removes intercept
        data=data.frame( cbind('time'=df$time,'delta'=df$delta,covariate) )
        fit=coxph(as.formula(paste("Surv(time,delta==1)", '~',
                                   paste(colnames(data[,-(1:2)]),collapse='+'))),data=data)
        newdata=NULL
        newdata[[1]]= data.frame(cbind(VTM(c(0,0),nrow(covariate)),covariate[,-(1:2)]) )
        newdata[[2]]= data.frame(cbind(VTM(c(1,0),nrow(covariate)),covariate[,-(1:2)]) )
        newdata[[3]]= data.frame(cbind(VTM(c(0,1),nrow(covariate)),covariate[,-(1:2)]) )
        survout=NULL;surv=NULL
        for (i in 1:length(newdata)){
          colnames(newdata[[i]])=colnames(covariate)
          survout[[i]]=survfit(fit,newdata=newdata[[i]] )
          surv[[i]]=apply(survout[[i]]$surv, 1, mean)
          # t=survout[[i]]$time
          # plot(t,surv[[i]],col=i,type='l',ylim=c(0,1))
          # par(new=T)
        }
        # legend('bottomleft',c('None','Peripheral','Central'),lwd=rep(1,length(newdata)),
        #        col = 1:length(newdata))
        names(survout)=c('none','pns','cns')

        event_table_obfs <- tryCatch({

          message("generating event_tables for adjusted survival curves")
          none_surv_table <- data.frame(status = "neuro_postNone",
                                   time = survout$none$time,
                                   n.censor = survout$none$n.censor,
                                   n.risk = survout$none$n.risk,
                                   n.event = survout$none$n.event)

          cns_surv_table <- data.frame(status = "neuro_postCentral",
                                  time = survout$cns$time,
                                  n.censor = survout$cns$n.censor,
                                  n.risk = survout$cns$n.risk,
                                  n.event = survout$cns$n.event)

          pns_surv_table <- data.frame(status = "neuro_postPeripheral",
                                  time = survout$pns$time,
                                  n.censor = survout$pns$n.censor,
                                  n.risk = survout$pns$n.risk,
                                  n.event = survout$pns$n.event)

          # combine tables
          event_table_surv_adjust <- rbind(none_surv_table, cns_surv_table, pns_surv_table)

          # mask and blur for obfuscation
          message("blurring event_tables for adjusted survival curves")
          event_surv_table_obfs <- blur_it(event_table_surv_adjust, vars = c("n.risk", "n.event", "n.censor"), blur_abs, mask_thres)

          # save survival model results to a list
          output=list('fit'=fit,
                      'surv_avg'=surv,
                      'survtable'=event_surv_table_obfs)

          # remove patient level data
          if (!is.null(output)) {
            output$fit$linear.predictors <- NULL
            output$fit$residuals <- NULL
            output$fit$n <- NULL
            output$fit$y <- NULL
            output$fit$nevent <- NULL
            output$fit$terms <- NULL

            output$survout$none$n <- NULL
            output$survout$none$surv <- NULL
            output$survout$none$cumhaz <- NULL
            output$survout$none$std.err <- NULL
            output$survout$none$lower <- NULL
            output$survout$none$upper <- NULL

            output$survout$cns$n <- NULL
            output$survout$cns$surv <- NULL
            output$survout$cns$cumhaz <- NULL
            output$survout$cns$std.err <- NULL
            output$survout$cns$lower <- NULL
            output$survout$cns$upper <- NULL

            output$survout$pns$n <- NULL
            output$survout$pns$surv <- NULL
            output$survout$pns$cumhaz <- NULL
            output$survout$pns$std.err <- NULL
            output$survout$pns$lower <- NULL
            output$survout$pns$upper <- NULL

          }
      },
      error = function(cond) {
        message(paste("Error when regressing", depend_var))
        message("Original error message:")
        message(cond)
        message('Skipping for now...')
        return(NULL) }
    )

      })
    }

  ## severe
  if (depend_var=='severe'){
    tryCatch({
      # if time to severe is <= tcut, delta = 1
      df$delta=(df$time_to_severe<=tcut)
      # if time_to_neuro <= tcut and neuro_post = "Peripheral", delta = 1
      df$delta.pns=(df$time_to_neuro<=tcut)*(df$neuro_post=="Peripheral")
      # if time_to_neuro <= tcut and neuro_post = "Central", delta = 1
      df$delta.cns=(df$time_to_neuro<=tcut)*(df$neuro_post=="Central")
      # df$delta.none=(df$time_to_neuro<=tcut)*(df$neuro_post=="None") #this is 0 because these patients don't have a time_to_neuro

      # remove the neuro_post var
      ind_vars.new=ind_vars[-1]
      independ_vars.new <- paste(ind_vars.new, collapse = ' + ')
      print('create covariate model')
      covariate <- model.matrix(as.formula(paste("delta", '~',
                                              independ_vars.new)),data=df)[,-1] #[-1] removes intercept
      data=data.frame( cbind('delta'=df$delta,'delta.pns'=df$delta.pns,
                             'delta.pns.both'=df$delta*df$delta.pns,
                             'delta.cns'=df$delta.cns,
                             'delta.cns.both'=df$delta*df$delta.cns,covariate) )

      # fit/predict
      print('fit models')
      fit=glm(as.formula(paste('delta~',
                               paste(colnames(covariate),collapse='+'))),
              family='binomial',data=data)
      pred=predict(fit,newdata=data,type = 'response')
      fit=glm(as.formula(paste('delta.pns~',
                               paste(colnames(covariate),collapse='+'))),
              family='binomial',data=data)
      pred.pns=predict(fit,newdata=data,type = 'response')
      fit=glm(as.formula(paste('delta.pns.both~',
                               paste(colnames(covariate),collapse='+'))),
              family='binomial',data=data)
      pred.pns.both=predict(fit,newdata=data,type = 'response')
      fit=glm(as.formula(paste('delta.cns~',
                               paste(colnames(covariate),collapse='+'))),
              family='binomial',data=data)
      pred.cns=predict(fit,newdata=data,type = 'response')
      fit=glm(as.formula(paste('delta.cns.both~',
                               paste(colnames(covariate),collapse='+'))),
              family='binomial',data=data)
      pred.cns.both=predict(fit,newdata=data,type = 'response')

      #pointwise mutal information
      print('calculate pmi')
      pmi.pns=mean(pred.pns.both/(pred*pred.pns))
      print(pmi.pns)
      pmi.pns.naive=mean(data$delta.pns.both)/(mean(data$delta)*mean(data$delta.pns))
      print(pmi.pns.naive)
      pmi.cns=mean(pred.cns.both/(pred*pred.cns))
      print(pmi.cns)
      pmi.cns.naive=mean(data$delta.cns.both)/(mean(data$delta)*mean(data$delta.cns))
      print(pmi.cns.naive)
      pmi=data.frame(pmi.pns,pmi.pns.naive,pmi.cns,pmi.cns.naive)
      print(pmi)

      # bootstrap confidence intervals
      time.start=Sys.time()
      print('starting bootstrap')
      rep=500
      print('create vv')
      vv=gen.bootstrap.weights(2022,n=nrow(data),num.perturb=rep)
      #print(vv)
      print('create a bootstrapped object')
      boot=apply(X = vv, MARGIN = 2, FUN = resam, data, covariate)
      print('end time')
      time.end=Sys.time()
      time=time.end-time.start
      print(time)
      print('calculating standard errors')
      se=apply(log(boot),1,mad) # mad = mean absolute difference

      output <- list('log.pmi'=log(pmi),'se'=se)
      print('calculating confidence intervals')
      output$lower_bound = exp(output$log.pmi-1.96*output$se)
      output$upper_bound = exp(output$log.pmi+1.96*output$se)
    },
    error = function(cond) {
      message(paste("Error when regressing", depend_var))
      message("Original error message:")
      message(cond)
      message('Skipping for now...')
      return(NULL) }
    )
  }
  output
}


VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


gen.bootstrap.weights=function(data.num, n, num.perturb=500){
  set.seed(data.num)
  sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))
  #weights = apply(index,2,function(x) tabulate(x,nbins=n))
  #list(index=index,weights=weights)
}

resam=function(vv,data,covariate){
  data=data[vv,]

  fit1=glm(as.formula(paste('delta~',
                           paste(colnames(covariate),collapse='+'))),
          family='binomial',data=data)
  pred=predict(fit1,newdata=data,type = 'response')

  fit2=glm(as.formula(paste('delta.pns~',
                           paste(colnames(covariate),collapse='+'))),
          family='binomial',data=data)
  pred.pns=predict(fit2,newdata=data,type = 'response')

  fit3=glm(as.formula(paste('delta.pns.both~',
                           paste(colnames(covariate),collapse='+'))),
          family='binomial',data=data)
  pred.pns.both=predict(fit3,newdata=data,type = 'response')

  fit4=glm(as.formula(paste('delta.cns~',
                           paste(colnames(covariate),collapse='+'))),
          family='binomial',data=data)
  pred.cns=predict(fit4,newdata=data,type = 'response')

  fit5=glm(as.formula(paste('delta.cns.both~',
                           paste(colnames(covariate),collapse='+'))),
          family='binomial',data=data)
  pred.cns.both=predict(fit5,newdata=data,type = 'response')

  pmi.pns=mean(pred.pns.both/(pred*pred.pns))
  pmi.pns.naive=mean(data$delta.pns.both)/(mean(data$delta)*mean(data$delta.pns))
  pmi.cns=mean(pred.cns.both/(pred*pred.cns))
  pmi.cns.naive=mean(data$delta.cns.both)/(mean(data$delta)*mean(data$delta.cns))

  out=c(pmi.pns,pmi.pns.naive,pmi.cns,pmi.cns.naive)

}

run_hosps <- function(neuro_patients,
                      neuro_pt_post,
                      non_neuro_patients,
                      both_pts,
                      both,
                      mask_thres,
                      blur_abs,
                      include_race,
                      currSiteId,
                      readmissions,
                      demo_processed,
                      obs_raw,
                      neuro_icds,
                      index_scores_elix,
                      pca_covariates) {

  ### CNS vs PNS
  message("Start CNS vs PNS Analysis")
  binary = FALSE
  print(binary == FALSE)

  neuro_types <- c("None", "Peripheral", "Central")

  demo_df <- demo_processed %>%
    filter(!patient_num %in% both_pts$patient_num) %>%
    left_join(distinct(select(neuro_patients, patient_num, neuro_type, time_to_neuro)),
      by = "patient_num"
    ) %>%
    replace_na(list(neuro_type = "None")) %>%
    mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types))

  scores_unique <- index_scores_elix %>%
    right_join0(demo_df, by = "patient_num") %>%
    left_join(pca_covariates, by = "patient_num")

  obfus_tables <- get_tables(
    neuro_types,
    demo_df,
    scores_unique,
    comorb_names_elix,
    blur_abs,
    mask_thres
  ) %>%
    lapply(function(x) mutate(x, site = currSiteId))

  ## -------------------------------------------------------------------------
  # Create individual tables for those who met outcome on admission and are
  # excluded from main survival analysis

  surv_exclude_pts <- function(df, time_to_outcome) {

    if (time_to_outcome == "time_to_severe") {
      demo_subset_df <- df %>%
        filter(time_to_severe == 0 | time_to_death == 0)
    } else if (time_to_outcome == "time_to_death") {
      demo_subset_df <- df %>%
        filter(time_to_death == 0)
    } else if (time_to_outcome == "time_to_first_discharge") {
      demo_subset_df <- df %>%
        mutate(death_before_outcome = case_when(time_to_death <= 0 ~ 1,
                                                TRUE ~ 0)) %>%
        filter(time_to_first_discharge == 0 & death_before_outcome == 0)
    } else if (time_to_outcome == "time_to_last_discharge") {
      demo_subset_df <- df %>%
        mutate(death_before_outcome = case_when(time_to_death <= 0 ~ 1,
                                                TRUE ~ 0)) %>%
        filter(time_to_last_discharge == 0 & death_before_outcome == 0)
    }

    scores_unique <- index_scores_elix %>%
      right_join0(demo_subset_df, by = "patient_num") %>%
      left_join(pca_covariates, by = "patient_num")


    obfus_tables <- tryCatch(
      {
        get_tables(
          neuro_types,
          demo_subset_df,
          scores_unique,
          comorb_names_elix,
          blur_abs,
          mask_thres
        ) %>%
        lapply(function(x) mutate(x, site = currSiteId))
      },
      error = function(cond) {
        message("Original error message:")
        message(cond)
        message("No data to subset. Skipping for now...")
        return(NULL) # return NA in case of error
      }
    )
    return(obfus_tables)
  }

  severe_adm <- surv_exclude_pts(demo_df, "time_to_severe")
  death_adm <- surv_exclude_pts(demo_df, "time_to_death")
  first_adm <- surv_exclude_pts(demo_df, "time_to_first_discharge")
  last_adm <- surv_exclude_pts(demo_df, "time_to_last_discharge")


  ## ---run-survival-models--------------------------------------------------------------------

  # Adults
  # only run the adult analysis if the site is not a pediatric only hospital
  if(currSiteId != c("BCH", "GOSH")) {
    tryCatch({
      surv_results30 <- run_coxregressions(df = scores_unique, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE)
      surv_results60 <- run_coxregressions(df = scores_unique, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE)
      surv_results90 <- run_coxregressions(df = scores_unique, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE)
    },
    error = function(cond) {
      message("Original error message:")
      message(cond)
      message("No data to subset. Skipping for now...")
      return(NULL) # return NA in case of error
    }
    )
  }
  # Pediatric
  surv_results_peds30 <- run_coxregressions(df = scores_unique, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE)
  surv_results_peds60 <- run_coxregressions(df = scores_unique, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE)
  surv_results_peds90 <- run_coxregressions(df = scores_unique, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE)

  ## ----save-results---------------------------------------------------------
  cpns_results <- c(obfus_tables, surv_results30, surv_results60, surv_results90, surv_results_peds30, surv_results_peds60, surv_results_peds90)

  results <- list(
    icd_tables = icd_tables,
    cpns_results = cpns_results,
    severe_adm = severe_adm,
    death_adm = death_adm,
    first_adm = first_adm,
    last_adm = last_adm,
    both = both
  )

  return(results)
}

get_elix_mat <- function(obs_raw, icd_version, t1 = -365, t2 = -15, map_type = "elixhauser") {
  ## -------------------------------------------------------------------------
  # for elixhauser
  comorb_names_elix <- get_quan_elix_names()

  # t1: earliest time point to consider comorbidities
  # t2: latest time point to consider comorbidities
  # example <- t1 = -365, and t2 = -1 will map all all codes up to
  # a year prior but before admission (admission = day 0)

  comorb_elix <- map_char_elix_codes(
    df = obs_raw,
    comorb_names = comorb_names_elix,
    icd_version = icd_version,
    t1 = t1,
    t2 = t2,
    map_type = map_type
  )

  index_scores_elix <- comorb_elix$index_scores %>%
    rename("elixhauser_score" = van_walraven_score)
  # van Walraven is a modification of Elixhauser comorbidity measure
  # doi.org/10.1097/MLR.0b013e31819432e5

  mapped_codes_table <- comorb_elix$mapped_codes_table

  comorb_list <- list(index_scores_elix = index_scores_elix,
                      mapped_codes_table = mapped_codes_table)

  return(comorb_list)
}

# return icd codes for patients during their first hospitalization as determined via 'first_out' or who are still 'in_hospital' when the data was pulled
temporal_neuro <- function(comp_readmissions, obs_raw, neuro_icds, readmissions, in_hospital) {

  obs_first_hosp <- comp_readmissions %>%
    left_join(., in_hospital, by = "patient_num") %>%
    filter(first_out | still_in_hospital == 1) %>%
    # days since admission the patient is out of hospital
    # here we just rename 'days_since_admission' to 'dsa'
    transmute(patient_num, dsa = days_since_admission) %>%
    right_join(obs_raw, by = "patient_num") %>%
    # keep the observations that occurred after or on the 'dsa'
    filter(days_since_admission <= dsa) %>%
    select(-dsa)


  list(
    obs_first_hosp = obs_first_hosp
  )
}
