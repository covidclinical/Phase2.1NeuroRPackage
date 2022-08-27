# get independent variables
get_ind_vars <- function(df, include_race, elix = "LPCA") {
  unique_cols <- apply(df, 2, function(x) length(unique(x)))

  comorb_names_elix <- get_quan_elix_names()

  if(elix == "LPCA") {
    ind_vars <- setdiff(
      c(
        "neuro_post", "sex", "age_group",
        "pre_admission_cns", "pre_admission_pns",
        paste0(".fittedPC", 1:10)
      ),
      names(unique_cols)[unique_cols == 1]
    )} else if(elix == "score") {
      ind_vars <- setdiff(
        c(
          "neuro_post", "sex", "age_group",
          "pre_admission_cns", "pre_admission_pns",
          "elixhauser_score"
        ),
        names(unique_cols)[unique_cols == 1]
      )
    } else if(elix == "ind") {
      ind_vars <- setdiff(
        c(
          "neuro_post", "sex", "age_group",
          "pre_admission_cns", "pre_admission_pns",
          comorb_names_elix$Abbreviation
        ),
        names(unique_cols)[unique_cols == 1]
      )
    }

  if (include_race) {
    ind_vars <- c(ind_vars, "race")
  }
  ind_vars
}

# 8.15.2022: Updated function from Xuan
# tcut = the censor time
run_coxregressions <- function(analysis, df, include_race = TRUE, tcut=60, blur_abs, mask_thres, is_pediatric = NULL, elix, survtable_interval) {

  print(analysis)

  # convert age group to character to remove factor levels
  df <- df %>%
    mutate(age_group = as.character(age_group))

  if(is_pediatric == FALSE){
    df <- df %>% filter(age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus"))
  } else if(is_pediatric == TRUE){
    df <- df %>% filter(age_group %in% c("00to02", "06to11", "12to17"))
  }

  ind_vars <- get_ind_vars(df, include_race, elix)

  deceased_reg_elix <-
    run_coxregression(df, 'deceased', ind_vars, tcut=tcut, blur_abs, mask_thres, survtable_interval)

  time_first_discharge_reg_elix <-
    run_coxregression(df, "time_to_first_discharge", ind_vars, tcut=tcut, blur_abs, mask_thres, survtable_interval)

  severe_reg_elix <-
    run_coxregression(df, 'severe', ind_vars, tcut=tcut, blur_abs, mask_thres, survtable_interval)

  list(
    severe_reg_elix = severe_reg_elix,
    deceased_reg_elix = deceased_reg_elix,
    time_first_discharge_reg_elix= time_first_discharge_reg_elix
  )
}


run_coxregression <-function(df, depend_var, ind_vars, tcut=60, blur_abs, mask_thres, survtable_interval) {
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
    # if time_to_first_discharge == time_to_first_death, make time_to_first_discharge == NA
    df <- df %>% mutate(time_to_first_discharge = ifelse(time_to_first_discharge == time_to_death,
                                                          time_to_first_discharge == NA,
                                                          time_to_first_discharge))
    # which.min will identify the column integer with the minimum time
    # (e.g. if time_to_first_discharge is the smallest, delta = 1)
    df$delta=apply(cbind(df$time_to_first_discharge,df$time_to_death,df$c),1,which.min)
    # remove those who meet outcome (discharge or death) on day of admission
    df <- df %>% filter(!(time == 0 & delta == 1))
    df <- df %>% filter(!(time == 0 & delta == 2))

  } else if (depend_var=='severe'){
    # if missing, temporarily set 1000 days
    df$time_to_severe[is.na(df$time_to_severe)]=1000
    df$time_to_neuro[is.na(df$time_to_neuro)]=1000
    # remove those who meet outcome on day of admission
    df <- df %>% filter(!time_to_severe == 0)
  }

  if (depend_var!='severe'){

    # fit initial survival model in order to create our table counts with the survminer package
    tryCatch(
      {
        fit = survival::survfit(Surv(time, delta)~neuro_post, data = df)

        # evaluate latest time point to prevent errors with using the 'break.time.by' function
        if(max(fit[['time']]) < survtable_interval) {
          max_time = max(fit[['time']])
          survtable_interval = floor(max_time / 10)*10
        } else {
          survtable_interval
        }

        survtable = ggrisktable(fit,
                                break.time.by = survtable_interval)

        # mask and blur for obfuscation
        message("obfuscating KM survtable")
        survtable_obfs <- blur_it(survtable[['data']] %>%
                                    select(strata, time, n.risk, cum.n.event, cum.n.censor, strata_size),
                                  vars = c("n.risk", "cum.n.event", "cum.n.censor", "strata_size"), blur_abs, mask_thres)

        #survtable_output = list(survtable = survtable_obfs)
      },
      error = function(cond) {
        message(paste("Error when regressing", depend_var))
        message("Original error message:")
        message(cond)
        message("Skipping for now...")
        return(NULL) # return NA in case of error
      }
    )


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

          # to print plot
          #t=survout[[i]]$time
          #plot(t,surv[[i]],col=i,type='l',ylim=c(0,1))
          #par(new=T)
        }
        # plot legend
        #legend('bottomleft',c('None','Peripheral','Central'),lwd=rep(1,length(newdata)),
        #       col = 1:length(newdata))

        # name output
        names(survout)=c('none','pns','cns')
        names(surv)=c('none','pns','cns')

        # save survival model results to a list
        output=list('fit'=fit,
                    #'survout' = survout, # this is not helpful due to all of the patient level probs that were averaged
                    # also survout time to event is not taking the strata into account so it's not helpful
                    'surv_avg'= surv,
                    'survtable' = survtable_obfs)

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
                             'delta.cns.both'=df$delta*df$delta.cns,covariate))

      # tabulate severe counts
      severe_table <- data.frame("not_severe" = nrow(data) - sum(data$delta),
                                 "severe" = sum(data$delta),
                                 "pns" = sum(data$delta.pns),
                                 "severe_pns" = sum(data$delta.pns.both),
                                 "cns" = sum(data$delta.cns),
                                 "severe_cns" = sum(data$delta.cns.both),
                                 "not_neuro" = nrow(data) - sum(data$delta.cns) - sum(data$delta.pns))

      message("blurring severe_table counts")
      severe_table_obfs <- blur_it(severe_table,
                              vars = c("not_severe", "severe", "pns", "severe_pns",
                                       "cns", "severe_cns", "not_neuro"), blur_abs, mask_thres)

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
      # naive is on the raw data; whereas the other metrics on on the predictions
      # note: NAs are returned when counts are 0s
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

      # save confidence intervals to output
      output <- list('log.pmi'=log(pmi),'se'=se) # standard error is in log - so we should take the log fist when calculating CIs
      output$pmi = pmi
      output$lower_bound = exp(output$log.pmi-1.96*output$se)
      output$upper_bound = exp(output$log.pmi+1.96*output$se)
      output$severe_table_obfs = severe_table_obfs
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

run_hosps <- function(both_pts,
                      both_counts,
                      mask_thres,
                      blur_abs,
                      include_race,
                      currSiteId,
                      readmissions,
                      demo_processed,
                      obs_raw,
                      neuro_icds) {


  ### CNS vs PNS
  message("Start CNS vs PNS Analysis")

  neuro_types <- c("None", "Peripheral", "Central")

  ##############################################################################

  ## process dataframes by pediatric or adult populations
  demo_df_adults <- demo_processed %>%
    filter(!patient_num %in% both_pts$patient_num,
           adult_ped == "adult") %>%
    left_join(distinct(select(neuro_patients, patient_num, neuro_type, time_to_neuro)),
              by = "patient_num"
    ) %>%
    replace_na(list(neuro_type = "None")) %>%
    mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types))

  demo_df_pediatrics <- demo_processed %>%
    filter(!patient_num %in% both_pts$patient_num,
           adult_ped == "pediatric") %>%
    left_join(distinct(select(neuro_patients, patient_num, neuro_type, time_to_neuro)),
              by = "patient_num"
    ) %>%
    replace_na(list(neuro_type = "None")) %>%
    mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types))

  scores_unique_adults <- comorb_adults$index_scores_elix %>%
    right_join0(demo_df_adults, by = "patient_num") %>%
    left_join(comorb_adults$pca_covariates, by = "patient_num")

  scores_unique_pediatrics <- comorb_pediatrics$index_scores_elix %>%
    right_join0(demo_df_pediatrics, by = "patient_num") %>%
    left_join(comorb_pediatrics$pca_covariates, by = "patient_num")

  # for elixhauser
  comorb_names_elix <- get_quan_elix_names()

  ## create demographic tables
  tableone_adults <- get_tables(
    neuro_types,
    demo_df_adults,
    scores_unique_adults,
    comorb_names_elix,
    blur_abs,
    mask_thres
  ) %>%
    lapply(function(x) mutate(x, site = currSiteId))

  tableone_pediatrics <- get_tables(
    neuro_types,
    demo_df_pediatrics,
    scores_unique_pediatrics,
    comorb_names_elix,
    blur_abs,
    mask_thres
  ) %>%
    lapply(function(x) mutate(x, site = currSiteId))


  ## corbidity tables
  # pivot comorbidity tables
  comorb_adults_pivot <- comorb_adults$index_scores_elix %>%
    select(-elixhauser_score) %>%
    pivot_longer(!patient_num, names_to = "Comorbidity") %>%
    filter(value == 1) %>%
    left_join(., comorb_adults$index_scores_elix %>%
                select(patient_num, elixhauser_score),
              by = 'patient_num') %>%
    select(-value)

  comorb_pediatrics_pivot <- comorb_pediatrics$index_scores_elix %>%
    select(-elixhauser_score) %>%
    pivot_longer(!patient_num, names_to = "Comorbidity") %>%
    filter(value == 1) %>%
    left_join(., comorb_pediatrics$index_scores_elix %>%
                select(patient_num, elixhauser_score),
              by = 'patient_num') %>%
    select(-value)

  # return demographic tables
  tableone_comorbidity_adults <- get_tables(
    neuro_types = c("None", "Peripheral", "Central"),
    demo_df = left_join(demo_df_adults, comorb_adults_pivot %>%
                          select(patient_num, Comorbidity), by = 'patient_num'),
    scores_unique = right_join0(comorb_adults_pivot, demo_df_adults, by = "patient_num"),
    comorb_names_elix = comorb_names_elix,
    blur_abs = blur_abs,
    mask_thres = mask_thres,
    group_var = "Comorbidity"
  )

  tableone_comorbidity_pediatrics <- get_tables(
    neuro_types = c("None", "Peripheral", "Central"),
    demo_df = left_join(demo_df_pediatrics, comorb_pediatrics_pivot %>%
                          select(patient_num, Comorbidity), by = 'patient_num'),
    scores_unique = right_join0(comorb_pediatrics_pivot, demo_df_pediatrics, by = "patient_num"),
    comorb_names_elix = comorb_names_elix,
    blur_abs = blur_abs,
    mask_thres = mask_thres,
    group_var = "Comorbidity"
  )

  ## -------------------------------------------------------------------------
  # Create individual tables for those who met outcome on admission and are
  # excluded from main survival analysis

  surv_exclude_pts <- function(df, time_to_outcome, index_scores_elix, pca_covariates) {

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

  # create demographic tables of excluded patients
  severe_adm_adults <- surv_exclude_pts(demo_df_adults, "time_to_severe", index_scores_elix = comorb_adults$index_scores_elix, pca_covariates = comorb_adults$pca_covariates)
  death_adm_adults <- surv_exclude_pts(demo_df_adults, "time_to_death", index_scores_elix = comorb_adults$index_scores_elix, pca_covariates = comorb_adults$pca_covariates)
  first_adm_adults <- surv_exclude_pts(demo_df_adults, "time_to_first_discharge", index_scores_elix = comorb_adults$index_scores_elix, pca_covariates = comorb_adults$pca_covariates)

  severe_adm_pediatrics <- surv_exclude_pts(demo_df_pediatrics, "time_to_severe", index_scores_elix = comorb_pediatrics$index_scores_elix, pca_covariates = comorb_pediatrics$pca_covariates)
  death_adm_pediatrics <- surv_exclude_pts(demo_df_pediatrics, "time_to_death", index_scores_elix = comorb_pediatrics$index_scores_elix, pca_covariates = comorb_pediatrics$pca_covariates)
  first_adm_pediatrics <- surv_exclude_pts(demo_df_pediatrics, "time_to_first_discharge", index_scores_elix = comorb_pediatrics$index_scores_elix, pca_covariates = comorb_pediatrics$pca_covariates)

  # save these patients to one list
  demo_excluded_on_admission <- list(
    severe_adm_adults = severe_adm_adults,
    death_adm_adults = death_adm_adults,
    first_adm_adults = first_adm_adults,
    severe_adm_pediatrics = severe_adm_pediatrics,
    death_adm_pediatrics = death_adm_pediatrics,
    first_adm_pediatrics = first_adm_pediatrics
    )

  # create demographic tables of non-excluded patients?
  # (actually this is kind of hard since patients could be in multiple analyses)



  ## ---run-survival-models--------------------------------------------------------------------

  # Adults
  # only run the adult analysis if the site is not a pediatric only hospital
  if(currSiteId != c("BCH", "GOSH")) {
    tryCatch({
      # LPCA
      surv_results_adults_lpca_30 <- run_coxregressions(analysis = "Adults_30_days_lpca", df = scores_unique_adults, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "LPCA", survtable_interval = 10)
      surv_results_adults_lpca_60 <- run_coxregressions(analysis = "Adults_60_days_lpca", df = scores_unique_adults, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "LPCA", survtable_interval = 20)
      surv_results_adults_lpca_90 <- run_coxregressions(analysis = "Adults_90_days_lpca", df = scores_unique_adults, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "LPCA", survtable_interval = 30)
      # elix score
      surv_results_adults_score_30 <- run_coxregressions(analysis = "Adults_30_days_score", df = scores_unique_adults, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "score", survtable_interval = 10)
      surv_results_adults_score_60 <- run_coxregressions(analysis = "Adults_60_days_score", df = scores_unique_adults, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "score", survtable_interval = 20)
      surv_results_adults_score_90 <- run_coxregressions(analysis = "Adults_90_days_score", df = scores_unique_adults, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "score", survtable_interval = 30)
      # independent predictors
      surv_results_adults_ind_30 <- run_coxregressions(analysis = "Adults_30_days_ind", df = scores_unique_adults, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "ind", survtable_interval = 10)
      surv_results_adults_ind_60 <- run_coxregressions(analysis = "Adults_60_days_ind", df = scores_unique_adults, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "ind", survtable_interval = 20)
      surv_results_adults_ind_90 <- run_coxregressions(analysis = "Adults_90_days_ind", df = scores_unique_adults, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = FALSE, elix = "ind", survtable_interval = 30)
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
  tryCatch({
    # LPCA
    surv_results_pediatrics_lpca_30 <- run_coxregressions(analysis = "Pediatrics_30_days_lpca", df = scores_unique_pediatrics, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "LPCA", survtable_interval = 10)
    surv_results_pediatrics_lpca_60 <- run_coxregressions(analysis = "Pediatrics_60_days_lpca", df = scores_unique_pediatrics, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "LPCA", survtable_interval = 20)
    surv_results_pediatrics_lpca_90 <- run_coxregressions(analysis = "Pediatrics_90_days_lpca", df = scores_unique_pediatrics, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "LPCA", survtable_interval = 30)
    # elix score
    surv_results_pediatrics_score_30 <- run_coxregressions(analysis = "Pediatrics_30_days_score", df = scores_unique_pediatrics, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "score", survtable_interval = 10)
    surv_results_pediatrics_score_60 <- run_coxregressions(analysis = "Pediatrics_60_days_score", df = scores_unique_pediatrics, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "score", survtable_interval = 20)
    surv_results_pediatrics_score_90 <- run_coxregressions(analysis = "Pediatrics_90_days_score", df = scores_unique_pediatrics, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "score", survtable_interval = 30)
    # independent predictors
    surv_results_pediatrics_ind_30 <- run_coxregressions(analysis = "Pediatrics_30_days_ind", df = scores_unique_pediatrics, include_race, tcut = 30, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "ind", survtable_interval = 10)
    surv_results_pediatrics_ind_60 <- run_coxregressions(analysis = "Pediatrics_60_days_ind", df = scores_unique_pediatrics, include_race, tcut = 60, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "ind", survtable_interval = 20)
    surv_results_pediatrics_ind_90 <- run_coxregressions(analysis = "Pediatrics_90_days_ind", df = scores_unique_pediatrics, include_race, tcut = 90, blur_abs = blur_abs, mask_thres = mask_thres, is_pediatric = TRUE, elix = "ind", survtable_interval = 30)
    },
  error = function(cond) {
    message("Original error message:")
    message(cond)
    message("No data to subset. Skipping for now...")
    return(NULL) # return NA in case of error
  }
  )

  ## ----save-results---------------------------------------------------------
  tryCatch({

  tableone_results <- list(tableone_adults = tableone_adults,
                           tableone_pediatrics = tableone_pediatrics
                        )

  tableone_comorbidity_results <- list(tableone_comorbidity_adults = tableone_comorbidity_adults,
                                       tableone_comorbidity_pediatrics = tableone_comorbidity_pediatrics
  )

  survival_results <- list(
    # adult results
    surv_results_adults_lpca_30  = surv_results_adults_lpca_30,
    surv_results_adults_lpca_60 = surv_results_adults_lpca_60,
    surv_results_adults_lpca_90 = surv_results_adults_lpca_90,
    surv_results_adults_score_30 = surv_results_adults_score_30,
    surv_results_adults_score_60 = surv_results_adults_score_60,
    surv_results_adults_score_90 = surv_results_adults_score_90,
    surv_results_adults_ind_30 = surv_results_adults_ind_30,
    surv_results_adults_ind_60 = surv_results_adults_ind_60,
    surv_results_adults_ind_90 = surv_results_adults_ind_90,
    # pediatric results
    surv_results_pediatrics_lpca_30 = surv_results_pediatrics_lpca_30,
    surv_results_pediatrics_lpca_60 = surv_results_pediatrics_lpca_60,
    surv_results_pediatrics_lpca_90 = surv_results_pediatrics_lpca_90,
    surv_results_pediatrics_score_30 = surv_results_pediatrics_score_30,
    surv_results_pediatrics_score_60 = surv_results_pediatrics_score_60,
    surv_results_pediatrics_score_90 = surv_results_pediatrics_score_90,
    surv_results_pediatrics_ind_30 = surv_results_pediatrics_ind_30,
    surv_results_pediatrics_ind_60 = surv_results_pediatrics_ind_60,
    surv_results_pediatrics_ind_90 = surv_results_pediatrics_ind_90
    )

  },
  error = function(cond) {
    message("Original error message:")
    message(cond)
    message("No data to subset. Skipping for now...")
    return(NULL) # return NA in case of error
  }
  )

  icd_tables = list(icd_tables_adults = comorb_adults$icd_tables,
                    icd_tables_pediatrics = comorb_pediatrics$icd_tables)

  # save all result objects
  results <- list(
    icd_tables = icd_tables,
    tableone_results = tableone_results,
    tableone_comorbidity_results = tableone_comorbidity_results,
    survival_results = survival_results,
    both_counts = both_counts,
    demo_excluded_on_admission = demo_excluded_on_admission
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

## process_comorb_data() computes LPCA by the specified patient population (pediatric or adult)
# df - dataframe - adult_obs or ped_obs depending on the population to analyze
# adult_pediatric - TRUE/FALSE value, where TRUE indicates to analyze the pediatric population
process_comorb_data <- function(df, icd_version, is_pediatric) {

  comorb_list <- get_elix_mat(df, icd_version)

  index_scores_elix <- comorb_list$index_scores_elix

  # ensure data is formatted correctly
  index_scores_elix$patient_num <- as.character(index_scores_elix$patient_num)

  index_scores_elix <- index_scores_elix %>%
    right_join0(select(demo_raw, patient_num, age_group), by = "patient_num")

  # filter by age
  if(is_pediatric == FALSE){
    index_scores_elix <- index_scores_elix %>%
      filter(age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus"))
  } else if (is_pediatric == TRUE) {
    index_scores_elix <- index_scores_elix %>%
      filter(age_group %in% c("00to02", "06to11", "12to17"))
  }

  # filter by age
  if(is_pediatric == FALSE){
    nstay_filtered <- nstay_df %>%
      filter(age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus"))
  } else if (is_pediatric == TRUE) {
    nstay_filtered <- nstay_df %>%
      filter(age_group %in% c("00to02", "06to11", "12to17"))
  }

  mapped_codes_table <- comorb_list$mapped_codes_table

  elix_mat <- cor(select(
    index_scores_elix,
    -c(patient_num, elixhauser_score, age_group)
  ))

  # individual comorbidity matrixes
  cns_pts <- neuro_patients %>%
    filter(pns_cns == "Central") %>%
    distinct(patient_num)

  pns_pts <- neuro_patients %>%
    filter(pns_cns == "Peripheral") %>%
    distinct(patient_num)

  index_scores_elix_cns <- index_scores_elix %>%
    filter(patient_num %in% cns_pts$patient_num)

  index_scores_elix_pns <- index_scores_elix %>%
    filter(patient_num %in% pns_pts$patient_num)

  elix_mat_cns <-
    cor(select(
      index_scores_elix_cns,
      -c(patient_num, elixhauser_score, age_group)
    ))

  elix_mat_pns <-
    cor(select(
      index_scores_elix_pns,
      -c(patient_num, elixhauser_score, age_group)
    ))

  elix_pca <- index_scores_elix %>%
    select(-elixhauser_score, -age_group) %>%
    tibble::column_to_rownames("patient_num") %>%
    as.matrix()

  ## logicstic pca
  # k = 10 principal components, m is solved for
  lpca_fit <- logisticPCA::logisticPCA(elix_pca, k = min(nrow(elix_pca), 10), m = 0)

  deviance_expl <- lpca_fit$prop_deviance_expl

  pca_covariates <- lpca_fit$PCs %>%
    data.frame() %>%
    `colnames<-`(paste0(".fittedPC", 1:10)) %>%
    tibble::rownames_to_column("patient_num")

  ind_covariates <- index_scores_elix %>%
    select(-elixhauser_score, -age_group)

  comorb_names_elix <- get_quan_elix_names()

  icd_tables <- get_tables(
    c("None", "Peripheral", "Central"),
    nstay_filtered,
    right_join0(index_scores_elix %>% select(-age_group), nstay_filtered, by = "patient_num"),
    comorb_names_elix,
    blur_abs,
    mask_thres,
    "concept_code"
  )

  ## obfuscate comorbidity table
  mapped_codes_table_obfus <- blur_it(mapped_codes_table, vars = 'n_patients', blur_abs, mask_thres)
  #mapped_codes_table_obfus <- mask_it(mapped_codes_table_obfus, var = 'n_patients', blur_abs, mask_thres)

  # remove categories with 0 patients
  mapped_codes_table_obfus <- mapped_codes_table_obfus %>%
    filter(!n_patients == 0)

  # remove age_group from index_scores_elix
  index_scores_elix <- index_scores_elix %>%
    select(-age_group)

  # save output to a list
  comorb_results <- list(icd_tables = icd_tables,
                         elix_mat = elix_mat,
                         elix_mat_cns = elix_mat_cns,
                         elix_mat_pns = elix_mat_cns,
                         deviance_expl = deviance_expl,
                         mapped_codes_table_obfus = mapped_codes_table_obfus,
                         # to be deleted after analysis
                         pca_covariates = pca_covariates,
                         index_scores_elix = index_scores_elix)

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
