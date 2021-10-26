# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("cancerdata")

library(dplyr)
library(cancerdata)
library(doParallel)
library(data.table)
library(MASS)
library(e1071)
library(randomForest)
library(gbm)
library(xgboost)
library(glmnet)

rm(list=ls())
{
  data(YOUNG1); check_valid <- T            # comment out when you want to validate accuracy on "VEER1"
  # data(VEER); raw_dat <- VEER; check=T        # use 25,000 genes
  data(VEER1); raw_dat <- VEER1; check=F    # use 5,000 genes
  veer_dat <- raw_dat %>% as_tibble()
  veer_dat$y <- ifelse(veer_dat$class=="NODM",0,1)
  if(check) dat <- veer_dat[,c(1:24481,24484)] else dat <- veer_dat[,c(1:4948,4951)] 
  
  ii <- 0
  b_na <- apply(dat,2,function(x){ii <<- ii+1; if(any(is.na(x))) return(ii)}) %>% unlist %>% as.vector
  dat <- dat[,-b_na] # remove na columns

  set.seed(1234)
  test_idx <- sample(1:78,size=78*0.3)
  x <- dat[,1:(ncol(dat)-1)]; y <- dat$y
  
  if(check_valid){
    train_x <- x
    test_x <- (YOUNG1 %>% as.tibble)[,YOUNG1 %>% as.tibble %>% colnames %in% (x %>% colnames)]
    train_y <- y
    test_y <- ifelse((YOUNG1 %>% as.tibble)$class=="NODM",0,1)
  } 
  else{
    train_x <- x[-test_idx,]
    test_x <- x[test_idx,]
    train_y <- y[-test_idx]
    test_y <- y[test_idx]
  }
  
  
  logit_v <- function(x){return(log(x/(1-x)))}
  expit_v <- function(x){return(exp(x)/(exp(x)+1))}
  
  library(doParallel)
  mycl <- makeCluster(6)
  registerDoParallel(mycl)
  
  pval <- foreach(i=seq(ncol(x)), .combine=c) %dopar% {
    fit<-cor(train_x[,i],train_y)
    return(abs(fit))
  }
  
  p_pred <- rep(0,14)
  
  ## select top 5,10,15,...,70 genes
  for(k in seq(5,70,by=5)){
    
    # sel.index<-order(pval)[seq(k)]
    sel.index<-order(pval,decreasing = T)[seq(k)]
    fdat<-as.data.frame(cbind(train_y,train_x[,sel.index]))
    
    if(check_valid) n <- 78 else n <- (78-length(test_idx))
    print(k)
    pr.est <- foreach(j =1:n, .combine=c) %dopar% {
      
      fit.final<-glm(train_y~.,family="binomial",data=fdat[-j,])
      xnew<-fdat[j,]
      return(predict(fit.final,xnew))
    }
    

    p_pred[k/5] <- table( (expit_v(pr.est)>0.5), train_y) %>% diag %>% sum / n
  }
  
  dt <- data.frame(x=seq(5,70,by=5),y=p_pred)
  stopCluster(mycl)
}

dt

d_threshold <- function(pr.est,y,type="logit"){
  check <- T
  threshold <- 0.7
  while(check){
    switch(type,
      logit=tb <- table( (expit_v(pr.est)>threshold), y),
      boost=tb <- table( (pr.est>threshold), y),
      da=tb <- table( (pr.est$posterior[,2]>threshold), y)
    )
    if(tb[1,2]/(tb[,2] %>% sum)<=0.1 | threshold<=0.1) check <- F
    else {check <- T; threshold <- threshold-0.05}
    
  }
  
  return(threshold)
}

# debug(d_threshold)
# undebug(d_threshold)

{
  set.seed(1234)
  ll_res <- vector("list",7)
  max_dt <- dt[which(dt$y==max(dt$y)),1] %>% min
  sel.index<-order(pval,decreasing = T)[seq(max_dt)]
  fdat<-as.data.frame(cbind(train_y,train_x[,sel.index]))
  
  fit.final<-glm(train_y~.,family=binomial,data=fdat)
  fdat2<-test_x[,sel.index]
  pr.est2<-predict(fit.final,fdat2)
  thr <- d_threshold(pr.est2,test_y)
  ll_res[[1]] <- table( (expit_v(pr.est2)>thr), test_y)
  
  # randomforest
  fit.final <- randomForest(factor(train_y)~.,data=fdat)
  pr.est2<-predict(fit.final,fdat2,type="prob")
  thr <- d_threshold(pr.est2[,2],test_y,type="boost")
  ll_res[[2]] <- table( pr.est2[,2]>thr, test_y)
  
  # gradient boosting
  fit.final <- gbm(train_y~.,data=fdat, distribution = "bernoulli", cv.folds = 10)
  pr.est2<-predict(fit.final,fdat2)
  thr <- d_threshold(pr.est2,test_y,type="logit")
  ll_res[[3]] <- table( pr.est2 %>% expit_v > thr, test_y)
  
  # xgboost
  params <- list(objective   = "binary:logistic",
                 eval_metric = "error",
                 max_depth   = 5,
                 eta         = 0.01)
  dtrain <- xgb.DMatrix(data = as.matrix(fdat[,-1]),label=train_y)
  dtest <- xgb.DMatrix(data = as.matrix(fdat2),label=test_y)
  # fit.final <- xgboost(data = dtrain, label = train_y, max.depth = 5,
  #         eta = 0.01, nthread = 2, nrounds = 500,objective = "binary:logistic")
  fit.cv <- xgb.cv(params,data=dtrain,nrounds = 500, nfold = 10, obje)
  fit.final <- xgb.train(params = params,data=dtrain,nrounds = 500, print_every_n = 10,watchlist=list(train=dtrain))
  pr.est2<-predict(object = fit.final,newdata = dtest)
  thr <- d_threshold(pr.est2,test_y,type="boost")
  ll_res[[4]] <- table( pr.est2>thr, test_y)
  
  fit.final<-lda(train_y~.,data=fdat)
  fdat2<-test_x[,sel.index]
  pr.est2<-predict(fit.final,fdat2)
  thr <- d_threshold(pr.est2,test_y,type="da")
  ll_res[[5]] <- table( pr.est2$posterior[,2]>thr, test_y)
  
  fit.final<-qda(train_y~.,data=fdat)
  fdat2<-test_x[,sel.index]
  pr.est2<-predict(fit.final,fdat2)
  thr <- d_threshold(pr.est2,test_y,type="da")
  ll_res[[6]] <- table( pr.est2$posterior[,2]>thr, test_y)
  
  fit.final<-naiveBayes(train_y~.,data=fdat,laplace = 1)
  fdat2<-test_x[,sel.index]
  pr.est2<-predict(fit.final,fdat2,type="raw")
  thr <- d_threshold(pr.est2[,2],test_y,type="boost")
  ll_res[[7]] <- table( pr.est2[,2]>thr, test_y)
  
  names(ll_res) <- c("logistic","rf","gbm","xgb","lda","qda","nb")
  print(ll_res)
  lapply(ll_res,function(x){print(diag(x) %>% sum / sum(x))}) %>% invisible
}


mycl <- makeCluster(6)
registerDoParallel(mycl)

ll_reg <- foreach(alp=seq(0,1,by=0.05), .packages = c("dplyr","glmnet")) %dopar% {
  cvelastic.fit <- cv.glmnet(fdat %>% dplyr::select(-train_y) %>% as.matrix, fdat$train_y,alpha=alp,parallel = T)
  best.lam <- cvelastic.fit$lambda.min
  # plot(cvelastic.fit)
  pr.est2<-predict(cvelastic.fit,s=best.lam,newx=fdat2 %>% as.matrix)
  thr <- d_threshold(pr.est2,test_y,type="boost")
  return(table(pr.est2>thr,test_y))
}
mt_reg <- matrix(c(seq(0,1,by=0.05),ll_reg %>% sapply(function(x){return(diag(x) %>% sum / sum(x))}) %>% invisible),
                 ncol=2)
max_alp <- mt_reg[which(mt_reg[,2]==max(mt_reg[,2])) %>% max,1]
ll_reg[[which(mt_reg[,2]==max(mt_reg[,2])) %>% max]] # lasso
stopCluster(mycl)

