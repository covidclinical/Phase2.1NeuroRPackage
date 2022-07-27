run_coxregressions <- function(df, include_race = TRUE, tcut=60) {
  ind_vars <- get_ind_vars(df, include_race)
  
  deceased_reg_elix <-
    run_coxregression(df, 'deceased', ind_vars, tcut=tcut)
  
  
  time_first_discharge_reg_elix <-
    run_coxregression(df, "time_to_first_discharge", ind_vars, tcut=tcut)
  
  severe_reg_elix <-
    run_coxregression(df, 'severe', ind_vars, tcut=tcut)
  
  list(
    severe_reg_elix = severe_reg_elix,
    deceased_reg_elix = deceased_reg_elix,
    time_first_discharge_reg_elix= time_first_discharge_reg_elix
  )
}


run_coxregression <-function(df, depend_var, ind_vars, tcut=60) {
    if (length(unique(df[, depend_var, drop = T])) <= 1)
      return(NULL)
    independ_vars <- paste(ind_vars, collapse = ' + ')
    
    if (depend_var=="deceased"){
    # df$delta=df$deceased
    df$c=pmin(df$days_since_admission,tcut)    
    df$time_to_death[is.na(df$time_to_death)]=1000
    df$time=apply(cbind(df$time_to_death,df$c),1,min)
    df$delta=as.numeric(df$time_to_death<=df$c)
    }else if (depend_var=='time_to_first_discharge'){
    df$c=pmin(df$days_since_admission,tcut) 
    df$time_to_first_discharge[is.na(df$time_to_first_discharge)]=1000
    df$time_to_death[is.na(df$time_to_death)]=1000
    df$time=apply(cbind(df$time_to_first_discharge,df$time_to_death,df$c),1,min)
    df$delta=apply(cbind(df$time_to_first_discharge,df$time_to_death,df$c),1,which.min)
    }else if (depend_var=='severe'){
      df$time_to_severe[is.na(df$time_to_severe)]=1000
      df$time_to_neuro[is.na(df$time_to_neuro)]=1000
    }
    
    if (depend_var!='severe'){
    tryCatch(
      {
        # fit=coxph(as.formula(paste("Surv(time,delta==1)", '~', independ_vars)), data = df)
        covariate=model.matrix(as.formula(paste("Surv(time,delta==1)", '~', 
                                                  independ_vars)),data=df)[,-1]
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
          # surv[[i]]=apply(survout[[i]]$surv, 1, mean)
          # t=survout[[i]]$time
          # plot(t,surv[[i]],col=i,type='l',ylim=c(0,1))
          # par(new=T)
          }
        # legend('bottomleft',c('None','Peripheral','Central'),lwd=rep(1,length(newdata)),
        #        col = 1:length(newdata))
        names(survout)=c('none','pns','cns')
        
        output=list('fit'=fit,'survout'=survout)
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
      df$delta=(df$time_to_severe<=tcut)
      df$delta.pns=(df$time_to_neuro<=tcut)*(df$neuro_post=="Peripheral")
      df$delta.cns=(df$time_to_neuro<=tcut)*(df$neuro_post=="Central")
      # df$delta.none=(df$time_to_neuro<=tcut)*(df$neuro_post=="None") #this is 0
      
      ind_vars.new=ind_vars[-1]
      independ_vars.new <- paste(ind_vars.new, collapse = ' + ')
      covariate=model.matrix(as.formula(paste("delta", '~', 
                                              independ_vars.new)),data=df)[,-1]
      
      data=data.frame( cbind('delta'=df$delta,'delta.pns'=df$delta.pns, 
                             'delta.pns.both'=df$delta*df$delta.pns,
                             'delta.cns'=df$delta.cns,
                             'delta.cns.both'=df$delta*df$delta.cns,covariate) )
      
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
      
      pmi.pns=mean(pred.pns.both/(pred*pred.pns))
      pmi.pns.naive=mean(data$delta.pns.both)/(mean(data$delta)*mean(data$delta.pns))
      pmi.cns=mean(pred.cns.both/(pred*pred.cns))
      pmi.cns.naive=mean(data$delta.cns.both)/(mean(data$delta)*mean(data$delta.cns))
      pmi=data.frame(pmi.pns,pmi.pns.naive,pmi.cns,pmi.cns.naive)
      
      # time.start=Sys.time()
      rep=500
      vv=gen.bootstrap.weights(2022,n=nrow(data),num.perturb=rep)
      boot=apply(vv,2,resam,data)
      # time.end=Sys.time()
      # time=time.end-time.start
      # print(time)
      se=apply(log(boot),1,mad)
      
      output=list('log.pmi'=log(pmi),'se'=se)
      # exp(output$log.pmi+1.96*output$se)
      # exp(output$log.pmi-1.96*output$se)
      },
      error = function(cond) {
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


resam=function(vv,data){
  data=data[vv,]
  
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
  
  pmi.pns=mean(pred.pns.both/(pred*pred.pns))
  pmi.pns.naive=mean(data$delta.pns.both)/(mean(data$delta)*mean(data$delta.pns))
  pmi.cns=mean(pred.cns.both/(pred*pred.cns))
  pmi.cns.naive=mean(data$delta.cns.both)/(mean(data$delta)*mean(data$delta.cns))
  
  out=c(pmi.pns,pmi.pns.naive,pmi.cns,pmi.cns.naive)
}
