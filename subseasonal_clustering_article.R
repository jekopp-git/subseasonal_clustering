library(evd) #for clusters function
library(data.table) # for efficient data manipulation
library(tictoc) #for execution time
library(hydroTSM) #for function time2season
library(ggplot2) #for plotting
library(bit64) #to handle HYBAS_ID field which is a large integer
library(zoo) #for specific function to replace NA values
library(foreach);library(iterators);library(parallel);library(doParallel) #for parallelisation

detectCores() #give the number of available cores for parallelisation
 
###Compute quantiles for all elements of thresholds
###then compute clusters for all run.lengths and thresholds combinations using clusters function from evd package
###Input: aggregated precipitation data PR_BASIN_XX (files are generated using HydroBASINS_agg_prec.ipynb)
###Output: "quantiles_XX.csv" and "clusters_XX.csv" (XX = 1 of the 9 HydroBASINS regions)

#Parameters
thresholds <- as.array(c(0.95,0.98,0.99,0.995)) #Threshold values for empirical quantiles
rl <- c(1,2,3,4,5) # Run length parameter in days
run.lengths <- array(rl)
rownames(run.lengths) <- as.character(rl)

#Function to compute quantiles for a vector of thresholds
quants <- function(z,th) {
  q <- as.array(quantile(z,th))
  names(q) <- paste0(rep("t",length(th)),
                     gsub(".","",as.character(th*100),fixed = TRUE))
  return(q)
}

#load precipitation files
files_PR <- list.files(path = ".",pattern="PR_",full.names = TRUE)

#compute quantiles
for (i in 1:length(files_PR)){
  name_PR <- files_PR[i]
  PR_data_t <- fread(name_PR)
  subname <- substring(files_PR[i],10,11)
  PR_data_t$time <- as.Date(PR_data_t$time, "%Y-%m-%d")
  setorder(PR_data_t,HYBAS_ID,time)
  qnts<-tapply(PR_data_t$TP,PR_data_t$HYBAS_ID, function(z){quants(z,thresholds)})
  qnts.unl <- unlist(qnts)
  tq <- strsplit(names(qnts.unl),"[.]")
  vq <- unname(qnts.unl)
  df_qnts <- data.frame(matrix(unlist(tq),nrow=length(tq),byrow = TRUE),stringsAsFactors=FALSE)
  df_qnts$X1 <- as.integer64(df_qnts$X1)
  df_qnts <- cbind(df_qnts,vq)
  df_qnts <- setNames(df_qnts,c("HYBAS_ID","Threshold","TP"))
  df_qnts <- data.table(df_qnts)
  #save file
  savename <- paste0("quantiles_",subname,".csv")
  fwrite(df_qnts, file = savename, row.names = FALSE)
  rm(qnts)
  rm(qnts.unl)
  rm(tq)
  rm(vq)
  rm(df_qnts)
  gc()
  
  #compute clusters
  cl.th.rl <- tapply(PR_data_t$TP,PR_data_t$HYBAS_ID,
                     function(z){apply(run.lengths,1,
                                       function(y){apply(quants(z,thresholds),1,
                                                         function(x){clusters(z,x,y)})})})
  
  #Output of clusters function is a list, arrange clusters data in data frame before saving it
  cl.unl <- unlist(cl.th.rl)
  t <- strsplit(names(cl.unl),"[.]")
  v <- unname(cl.unl)
  temp_df <- data.frame(matrix(unlist(t),nrow=length(t),byrow = TRUE),stringsAsFactors=FALSE)
  temp_df$X5 <- PR_data_t[as.numeric(temp_df$X5),2]
  temp_df$X1 <- as.integer64(temp_df$X1)
  temp_df$X2 <- as.integer(temp_df$X2)
  temp_df <- cbind(temp_df,v)
  temp_df <- setNames(temp_df,c("HYBAS_ID","Run_length","Threshold","Cluster_id","time","TP"))
  rm(cl.th.rl)
  rm(cl.unl)
  rm(t)
  rm(v)
  gc()
  
  #save file
  temp_df <- data.table
  savename <- paste0("clusters",subname,".csv")
  fwrite(temp_df, file = savename, row.names = FALSE)
}

###Create files with precipitation accumulations and number of extreme of events for all time windows
###Each run length r is treated separately to avoid large files
###Input: aggregated precipitation data PR_BASIN_XX (files are generated using HydroBASINS_agg_prec.ipynb)
###and clusters data (clusters_XX)
###Output: tw_rY_XX_midlat (XX = 1 of the 9 HydroBASINS regions, Y = run length in days)

files_PR <- list.files(path = ".",pattern="PR_",full.names = TRUE)
files_CL <- list.files(path = ".",pattern="clusters_",full.names = TRUE)

#Apply only to subselection of catchments (see HydroBASINS_agg_prec.ipynb)
midlat_TP10 <- fread("Catchments_midlat_TP10.csv")

#Define time window lengths
windows <- c(14,21,28)

#Set some columns names
cols <- c("t95","t98","t99","t995")
cols2 <- c("t95n","t98n","t99n","t995n")

#Compute number of extreme events and precipitation accumulations in time window w
#Repeat for each run length r
r <- 1
for (i in 1:length(files_PR)){
  name_PR <- files_PR[i]
  name_CL <- files_CL[i]
  subname <- substring(files_CL[i],12,13)
  
  #Load data
  PR_data_t <- fread(name_PR)
  PR_data_t$time <- as.Date(PR_data_t$time, "%Y-%m-%d")
  setorder(PR_data_t,HYBAS_ID,time)
  
  df <- fread(name_CL)
  df$time <- as.Date(df$time, "%Y-%m-%d")
  
  #Select only midlatitude catchments
  df_midlat <- merge(df,midlat_TP10[,c("HYBAS_ID")],by=c("HYBAS_ID"))
  PR_data_midlat <- merge(PR_data_t,midlat_TP10[,c("HYBAS_ID")],by=c("HYBAS_ID"))
  rm(df)
  rm(PR_data_t)
  gc()
  
  #Compute cumulated precipitations over window lengths (windows)
  for(i in 1:length(windows)) {
    PR_data_midlat[,paste0("TPw",windows[i]):=
                     Reduce(`+`, shift(TP, -(windows[i]-1):0, fill=0)),by=.(HYBAS_ID)]
  }
  
  #Select run length of n days and only first peak of cluster
  peak_time <- dcast(df_midlat[Run_length==r,.SD[1],by=.(HYBAS_ID,Threshold,Cluster_id)],
                     HYBAS_ID + time ~ Threshold, value.var = c("time"))
  
  #Merge the identified peaks and the cumulated precipitations
  tw <- merge(PR_data_midlat,peak_time,by=c("HYBAS_ID","time"),all.x=TRUE)
  
  #Create new columns identifying peaks dates by 1 and NAs with 0
  tw[,paste0(cols,"n"):= lapply(.SD,function(x){as.integer(ifelse(is.na(x),"0","1"))}),.SDcols=cols]
  
  #Remove the initial columns (no longer needed)
  tw[,(cols):=NULL]
  
  #Compute the number of peaks over window lengths for each threshold (cumulated peaks)
  for(j in 1:length(windows)) {
    tw[,paste0("PK",cols2,windows[j]):=
         lapply(.SD,function(x){Reduce(`+`, shift(x, 0:(windows[j]-1),type="lead", fill=0))}),
       .SDcols=cols2,by=.(HYBAS_ID)]
  }
  
  #find season
  tw[,Season:=time2season(time)]
  
  #save file
  savename <- paste0("tw_r",r,"_",subname,"_midlat.csv")
  fwrite(tw, file = savename, row.names = FALSE)
}


###Algorithm to select non overlapping clustering episodes and compute cluster metrics
###Input: files with precipitation accumulations and number of extreme of events for all time windows
###tw_rY_XX_midlat (XX = 1 of the 9 HydroBASINS regions, Y = run length in days)
###Output: files with all clustering episodes and metrics by catchment
###One file for each parameters combination (thresholds, run length, time window)
###("eps_rY_tZ_wW_XX", Z = threshold, W = time window)
###Repeat for each run length r

files <- list.files(path = ".",pattern="tw_r2",full.names = TRUE)
r <- substring(files[1],7,7)

N = 50 #number of clustering episodes to extract
r <- substring(files[1],7,7)
windows <- c(14,21,28)
thresholds <- c(98,99)
Inc_list = c(15,20,30,50) #number of episodes to compute metrics based on incenter method

#Incenter function (cf. Sitarz, 2013)
incenter <- function(inc_N) {
  inc_1 <- 1
  inc_2 <- sqrt(2)+1
  inc_3 <- ((sqrt(2)+1)*(sqrt(3)+2)-(sqrt(3)+1))
  v <- c(length=inc_N)
  for (i in inc_N:1)
    if (i==inc_N) {
      v[i] <- inc_1
    } else if (i==inc_N-1) {
      v[i] <- inc_2
    } else if (i==inc_N-2) {
      v[i] <- inc_3
    } else {
      v[i] <- 3*v[i+1]-3*v[i+2]+v[i+3]
    }
  s <- v/v[1]
  return(s)
}

for (w in windows) {
  k <- 2*w*N
  for (th in thresholds) {
    for (file in files){
      name <- file
      region <- substring(file,9,10)
      tw <- fread(name)
      tw$time <- as.Date(tw$time, "%Y-%m-%d")
      v <- c("HYBAS_ID","time",paste0("TPw",w),paste0("PKt",th,"n",w))
      
      #Select the first N non-overlapping episodes with the largest precipitation accumulations
      tw_sel = tw[, v, with=FALSE]
      setorderv(tw_sel,c("HYBAS_ID",paste0("TPw",w)),c(1,-1))
      tw_sel <- tw_sel[, .SD[1:k],by=.(HYBAS_ID)]
      
      for (i in 1:N) {
        tw_sel[tw_sel[, .I[1], by=.(HYBAS_ID)]$V1, N1:=i]
        tw_sel <- tw_sel[,`:=`(timesup = time[1]+w-1 , timeinf = time[1]-w+1),by=.(HYBAS_ID)]
        if (i==1) {
          epsTP <- tw_sel[N1==i,]
        } else {
          epsTP <- rbind(epsTP,tw_sel[N1==i,])
        }
        tw_sel <- tw_sel[(time > timesup) | (time < timeinf),.SD,by=.(HYBAS_ID)]
      }
      
      epsTP[,`:=`(timeinf=NULL,timesup=NULL)]
      
      #Select the first N non-overlapping episodes with the largest number of events
      tw_sel = tw[, v, with=FALSE]
      setorderv(tw_sel,c("HYBAS_ID",paste0("PKt",th,"n",w),paste0("TPw",w)),c(1,-1,-1))
      tw_sel <- tw_sel[, .SD[1:k],by=.(HYBAS_ID)]
      
      for (i in 1:N){
        tw_sel[tw_sel[, .I[1], by=.(HYBAS_ID)]$V1, N2:=i]
        tw_sel <- tw_sel[,`:=`(timesup = time[1]+w-1 , timeinf = time[1]-w+1),by=.(HYBAS_ID)]
        if (i==1) {
          epsPK <- tw_sel[N2==i,]
        } else {
          epsPK <- rbind(epsPK,tw_sel[N2==i,])
        }
        tw_sel <- tw_sel[(time > timesup) | (time < timeinf),.SD,by=.(HYBAS_ID)]
      }
      
      epsPK[,`:=`(timeinf=NULL,timesup=NULL)]
      eps <- rbind(epsTP,epsPK,fill=TRUE)
      gc()
      
      ##S_f' calculations
      setorderv(eps,c("HYBAS_ID","N1"),c(1,1),na.last=TRUE)
      for (m in Inc_list){
        inc <- incenter(m)
        for (i in 1:m) {
          eps[N1==i, paste0("PK_INC1_",m) := get(paste0("PKt",th,"n",w))*inc[i]]
        }
        eps[N1<=m, paste0("Sfp_INC",m) := sum(get(paste0("PK_INC1_",m))) ,by=c("HYBAS_ID")]
        setnafill(eps,"locf",cols=c(paste0("Sfp_INC",m)))
      }
      
      ##S_f and S_r calculations
      setorderv(eps,c("HYBAS_ID","N2"),c(1,1),na.last=TRUE)
      for (m in Inc_list){
        inc <- incenter(m)
        for (i in 1:m) {
          eps[N2==i, paste0("PK_INC2_",m) := get(paste0("PKt",th,"n",w))*inc[i]]
        }
        eps[N2<=m, paste0("Sf_INC",m) := sum(get(paste0("PK_INC2_",m))) ,by=c("HYBAS_ID")]
        setnafill(eps,"locf",cols=c(paste0("Sf_INC",m)))
        eps[,paste0("Sr_INC",m) := get(paste0("Sfp_INC",m))/get(paste0("Sf_INC",m)),by=c("HYBAS_ID")]
      }
      
      #Maximum number of extreme events by catchment
      eps[,PK_Max := max(get(paste0("PKt",th,"n",w))),by=c("HYBAS_ID")]
      
      #remove temp variables
      eps[,`:=`(PK_INC1_15=NULL,PK_INC2_15=NULL,PK_INC1_20=NULL,PK_INC2_20=NULL,
                PK_INC1_30=NULL,PK_INC2_30=NULL,PK_INC1_50=NULL,PK_INC2_50=NULL)]
      
      #save file
      savename <- paste0("eps_r",r,"t",th,"w",w,"_",region,".csv")
      fwrite(eps, file = savename, row.names = FALSE)
      
    }
  }
}

###Compute index of dispersion for each catchment and make a significance test (DI > 1)
###at various confidence intervals based on resampling
###Input: files with precipitation accumulations and number of extreme of events for all time windows
###tw_rY_XX_midlat (XX = 1 of the 9 HydroBASINS regions, Y = run length in days)
###Output: files with all clustering episodes and metrics by catchment
###One file for each parameters combination (thresholds, run length, time window)
###("eps_rY_tZ_wW_XX", Z = threshold, W = time window)
###Repeat for each run length r

files <- list.files(path = ".",pattern="tw_r1",full.names = TRUE)
ri <- substring(files[1],7,7)

windows <- c(14,21,28)
thresholds <- c(98,99)
N <- 14699 #number of days in precipitation time series
rep <- 10000 #number of samples for resampling
nb_cores <- floor(detectCores()/3 - 1) #get number of cores for parallel computations of samples
registerDoParallel(cores=nb_cores) 

#compute dispersion index DI for each threshold, run length and time windows
for (w in windows) {
  for (file in files){
    name <- file
    region <- substring(file,9,10)
    temp <- fread(name)
    temp$time <- as.Date(temp$time, "%Y-%m-%d")
    th_fields <- sapply(thresholds,function(x){paste0("PKt",x,"n",w)})
    thn <- sapply(thresholds,function(x){paste0("t",x,"n")})
    v <- c("HYBAS_ID","time",th_fields)
    
    test <- temp[temp[, .I[seq(1, .N, w-1)], by=.(HYBAS_ID)]$V1,v,with=FALSE]
    test[, paste0("mean_",th_fields) := lapply(.SD, mean),.SDcols = th_fields,by=.(HYBAS_ID)]
    test[, paste0("var_",th_fields) := lapply(.SD, var),.SDcols = th_fields,by=.(HYBAS_ID)]
    for (th in th_fields) {
      test[, paste0("DI_",substr(th,4,5)) := get(paste0("var_",th))/get(paste0("mean_",th)),by=.(HYBAS_ID)]
    }
    
    savename <- paste0("DI_r",ri,"t",thresholds[1],"_",thresholds[2],"w",w,"_",region,".csv")
    fwrite(test, file = savename, row.names = FALSE)
  }
}
gc()


#function to generate N samples from Poisson distribution with same parameter lambda
make_DI <- function(h,lambda,N,rep,thresholds,windows,ri) {
  dt <- data.table(HYBAS_ID=integer64(), t=numeric(), w=numeric(), r=numeric(), DI=numeric())
  for (th in thresholds) {
    v <- paste0("m_",th)
    for (wi in windows) {
      process <- rpois(N*rep,as.numeric(lambda[HYBAS_ID==h,v,with=FALSE]))
      mat_proc <- matrix(process,nrow = N,ncol = rep)
      bin_size <- rep(1:ceiling(N/wi),each = wi,len = N)
      sum_bin <- rowsum(mat_proc,bin_size)
      means <- colMeans(sum_bin)
      vr <- apply(sum_bin,2,var)
      DI_vect <- vr/means
      temp <- data.table(HYBAS_ID=integer64(), t=numeric(), w=numeric(), r=numeric(), DI=numeric())[1:rep]
      temp[,`:=`(HYBAS_ID = h, t = th, w = wi, r = as.numeric(ri), DI = DI_vect)]
      dt <- rbindlist(list(dt,temp))
    }
  }
  return (dt)
}

#Create DI distribution based on N samples for all catchments
#Those files are large and are not kept
for (file in files){
  tic()
  dta <- data.table(HYBAS_ID=integer64(), t=numeric(), w=numeric(), r=numeric(), DI=numeric())
  name <- file
  region <- substring(file,9,10)
  tmp <- fread(name)
  for (th in thresholds) {
    tmp[,paste0("m_",th) := mean(get(paste0("t",th,"n"))), by=.(HYBAS_ID)]
  }
  lambda <- unique(tmp[,.(HYBAS_ID,m_98,m_99)])
  rm(tmp)
  gc()
  dta <- foreach(h=lambda[,HYBAS_ID],.packages = c('data.table','tictoc','bit64'),.combine = rbind) %dopar% {
    dt <- make_DI(h,lambda,N,rep,thresholds,windows,ri)
    return(dt)
  }
  savename <- paste0("dist_DI_r",ri,"_",region,".csv")
  fwrite(dta, file = savename, row.names = FALSE)
  gc()
  toc()
}

#Calculate chosen quantiles (p) of DI distributions and aggregate the result in one file (dist_disp_index_qt.csv)
#This file takes less memory than all distribution files (dist_DI_)
files <- list.files(path = ".",pattern="dist_DI",full.names = TRUE)
p <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5, 99.9, 100)/100
names(p) <- sapply(p,function(x) {paste0("qt",x)})
dta <- data.table(HYBAS_ID=integer64(), t=numeric(), w=numeric(), r=numeric())
dta[, names(p) := as.list(p)]
for (file in files){
  name <- file
  tmp <- fread(name)
  dt <- tmp[,names(p):= as.list(quantile(DI,p)),by=.(HYBAS_ID,t,w,r)]
  dt <- unique(tmp[,.SD,.SDcols = !c("DI")])
  dta <- rbindlist(list(dta,dt))
}
fwrite(dta, file = "dist_disp_index_qt.csv", row.names = FALSE)

#Merge actual dispersion index (DI_r...) with quantiles of distribution from resampling (dist_disp_index_qt)
#in a single dataframe, then test significance

dta <- fread("dist_disp_index_qt.csv")
files <- list.files(path = ".",pattern="DI_",full.names = TRUE)
thresholds <- c(98,99)
dt <- data.table(HYBAS_ID=integer64(), t=integer(), w=integer(), r=integer(), DI=numeric())

#merge and rearrange all data
for (file in files){
  name <- file
  tmp <- fread(name)
  ri <- as.integer(substring(file,7,7))
  wi <- as.integer(substring(file,15,16))
  tmp <- unique(tmp[,.SD,.SDcols = c("HYBAS_ID","DI_98","DI_99")])
  tmp <- tmp[,`:=`(w = wi,r = ri,t = thresholds[1])]
  tmp2 <- copy(tmp)
  tmp2[,t := thresholds[2]]
  tmp2[,DI_98 := NULL]
  setnames(tmp2,c("DI_99"),c("DI"))
  tmp[,DI_99 := NULL]
  setnames(tmp,c("DI_98"),c("DI"))
  tmp <- rbindlist(list(tmp,tmp2))
  dt <- rbindlist(list(dt,tmp),use.names = TRUE)
}

dta <- merge(dta,dt,on=.(HYBAS_ID,t,w,r,t),all.x = TRUE)

#make significance test (test if actual DI is greater and lower of chosen quantiles)
cols = names(dta)[5:21]
dta[, paste0(cols,"_g") := lapply(.SD,function(x){DI>x}), .SDcols = cols]
dta[, paste0(cols,"_l") := lapply(.SD,function(x){DI<x}), .SDcols = cols]
fwrite(dta, file = "disp_ind_all.csv", row.names = FALSE)

