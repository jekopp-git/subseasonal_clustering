#Resampling script for significance test of the clustering metric (S_cl)

library(evd) #for clusters function
library(data.table) # for efficient data manipulation
library(bit64) #to handle HYBAS_ID field which is a large integer
library(tictoc)
library(foreach);library(iterators);library(parallel);library(doParallel) #for parallelisation
library(doRNG) #to use the same seed for reproducibility

# make one precipitation quantile file
files_qt <- list.files(path = ".",pattern="quantiles_",full.names = TRUE)
qt_data <- data.table(HYBAS_ID=integer64(), Threshold=character(), TP=numeric())
for (file in files_qt) {
   tmp <- fread(file)
   qt_data <- rbindlist(list(qt_data,tmp),use.names = TRUE)
}
fwrite(qt_data, file="qt_all.csv", row.names = FALSE)

#make one precipitation file by catchment
 
files_PR <- list.files(path = ".",pattern="PR_",full.names = TRUE)
nb_cores <- detectCores() -1 #get number of cores for parallel computations of samples
registerDoParallel(cores=nb_cores)

for (file in files_PR){
  temp <- fread(file, colClasses = c("integer64","character","numeric"))
  catchments <- temp[order(HYBAS_ID), head(.SD, 1), by = HYBAS_ID][[1]]
  foreach (i=1:length(catchments), .packages = c('data.table','bit64')) %dopar% {
    tmp <- temp[temp$HYBAS_ID==catchments[i],]
    fwrite(tmp, file=paste0("pr_",as.character(catchments[i]),".RData"), row.names = FALSE)
  }
}

detectCores()
#parameters
th <- 99
rl <- 2
w <- 21
N <- 50 #number of clustering episodes to extract
k <- 2*w*N
nb_cores <- floor(detectCores()/3 - 1) #get number of cores for parallel computations of samples
registerDoParallel(cores=nb_cores) 
trials <- 1000
set.seed(2021) #to use the same seed for reproducibility

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

#Function for calculating the metrics

calc_metrics <- function(prec_data,basin,th,rl,w,qt,k,N) {
  #Compute high-frequency clusters
  cl_temp <- clusters(data=prec_data$TP, qt, r=rl)
  t <- strsplit(names(unlist(cl_temp)),"[.]")
  cl <- data.table(matrix(unlist(t),nrow=length(t),byrow = TRUE),stringsAsFactors=FALSE)
  rm(cl_temp,t)
  cl$V2 <- prec_data[as.numeric(cl$V2),1]
  cl <- setNames(cl,c("cluster_id","time"))
  
  #Select first day/extreme of each high-frequency cluster (declustering)
  cl_init <- cl[,.SD[1],by=.(cluster_id)]
  rm(cl)
  
  #Identify independant extremes with 1
  prec_data <- merge(prec_data,cl_init,by=c("time"),all.x=TRUE)
  rm(cl_init)
  prec_data[, "ext" := as.integer(ifelse(is.na(cluster_id),"0","1"))]
  
  #Count number of independant extremes over time window w
  prec_data[,"n":= Reduce(`+`, shift(ext, 0:(w-1),type="lead", fill=0))]
  
  #Compute accumulation over window w
  prec_data[,"acc":= Reduce(`+`, shift(TP, 0:(w-1), type="lead", fill=0))]
  
  ###Algorithm to select non overlapping clustering episodes and compute cluster metrics
  prec_data[,`:=`(cluster_id=NULL,TP=NULL,ext=NULL)]
  
  #Select the first N non-overlapping episodes with the largest precipitation accumulations
  setorderv(prec_data,"acc",-1)
  acc_sel <- prec_data[, .SD[1:k]]
  cl_acc <- data.table(time=as.Date(character()), n=integer(), acc=numeric(), acc_ep=integer(),
                       timesup=as.Date(character()), timeinf=as.Date(character()))
  for (i in 1:N) {
    acc_sel[1, acc_ep:=i]
    acc_sel <- acc_sel[,`:=`(timesup = time[1]+w-1 , timeinf = time[1]-w+1)]
    temp <- acc_sel[acc_ep==i,]
    cl_acc <- rbindlist(list(cl_acc,temp))
    acc_sel <- acc_sel[(time > timesup) | (time < timeinf),.SD]
  }
  
  #Select the first N non-overlapping episodes with the largest number of events
  setorderv(prec_data,c("n","acc"),c(-1,-1))
  n_sel <- prec_data[, .SD[1:k]]
  cl_n <- data.table(time=as.Date(character()), n=integer(), acc=numeric(), n_ep=integer(),
                     timesup=as.Date(character()), timeinf=as.Date(character()))
  for (i in 1:N){
    n_sel[1, n_ep:=i]
    n_sel <- n_sel[,`:=`(timesup = time[1]+w-1 , timeinf = time[1]-w+1)]
    temp <- n_sel[n_ep==i,]
    cl_n <- rbindlist(list(cl_n,temp))
    n_sel <- n_sel[(time > timesup) | (time < timeinf),.SD]
  }
  ##Metrics calculations
  inc <- incenter(N)
  n <- cl_acc[["n"]]
  Sacc <- sum(n*inc)
  n <- cl_n[["n"]]
  Sn <- sum(n*inc)
  Sr <- Sacc/Sn
  result <- data.table(HYBAS_ID = basin, S_acc = Sacc, S_n = Sn, S_r = Sr)
  
  return(result)
}

#Make 1000 (ntrials) samples of precipitation time series for each catchment
#Then calculate metrics for each sample

qt_data <- fread("qt_all.csv")
catchments_list <- fread("Catchments_midlat_TP10.csv")
qt_data <- qt_data[HYBAS_ID %chin% catchments_list$HYBAS_ID]
dt_perm <- data.table(HYBAS_ID=integer64(), S_acc=numeric(), S_n=numeric(), S_r=numeric())
dt_obs <- data.table(HYBAS_ID=integer64(), S_acc=numeric(), S_n=numeric(), S_r=numeric())

files_pr <- list.files(pattern="pr_",full.names = FALSE)
c_list <- paste0("pr_",as.character(catchments_list[["HYBAS_ID"]]),".RData")
files_pr <- files_pr[files_pr %in% c_list]

for (i in 1:length(files_pr)) {
  h <- as.integer64(substr(files_pr[i],4,13))
  qt <- qt_data[qt_data$HYBAS_ID==h & qt_data$Threshold==paste0("t",th),][[3]]
  prec_data <- fread(files_pr[i],select = c("time","TP"))
  prec_data$time <- as.Date(prec_data$time, "%Y-%m-%d")
  obs <- calc_metrics(prec_data,h,th,rl,w,qt,k,N)
  dt_obs <- rbindlist(list(dt_obs,obs))
  perm <- foreach(icount(trials), .packages = c('evd','data.table','tictoc','bit64'), .combine=rbind) %dorng% {
    prec <- sample(prec_data$TP)
    prec_data[,TP:=prec]
    res <- calc_metrics(prec_data,h,th,rl,w,qt,k,N)
  }
  dt_perm <- rbindlist(list(dt_perm,perm))
}

fwrite(dt_obs, file="obs_metrics_all.RData", row.names = FALSE)
fwrite(dt_perm, file="perm_metrics_all.RData", row.names = FALSE)

stopImplicitCluster()

#Calculate empirical p-value (pval) using empirical cumulative distribution function of S_cl (S_n)

perm <- fread("perm_metrics_all.RData")
obs <- fread("obs_metrics_all.RData")
cols <- c("S_acc","S_n","S_r")
perm[, paste0(cols, "_m") := lapply(.SD, mean), by = .(HYBAS_ID), .SDcols = cols]
perm[, paste0(cols, "_sd") := lapply(.SD, sd), by = .(HYBAS_ID), .SDcols = cols]

catchments <- obs[,HYBAS_ID]
for (i in 1:length(catchments)) {
  data <- perm[HYBAS_ID == catchments[i],S_n]
  val <- obs[HYBAS_ID == catchments[i],S_n]
  obs[HYBAS_ID == catchments[i], pval:=1-ecdf(data)(val)]
}

fwrite(obs, file="obs_metrics_all.RData", row.names = FALSE)