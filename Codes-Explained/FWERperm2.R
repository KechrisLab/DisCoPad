#' Estimation of permutation-based metabolome-wide significance level (MWSL) and corresponding effective number of tests (ENT)
#'
#' Parametric approximation methods based on the (log-)multivariate Normal distribution are employed to estimate a stable permutation-based significance level to be used for multiple tests adjustement in the case of correlated metabolomic variates.
#' @param \code{outcome} a vector or a data frame of \code{n} samples. Possible outcome types are: continuous (both symmetric and skewed), discrete binary (0/1), discrete countable, or survival time-to-event
#' @param \code{features} a data frame of \code{n} samples (rows) and \code{M} features (columns)
#' @param \code{confounders} an optional data frame of \code{n} samples (rows) and \code{P} fixed effects confounders (columns), default=\code{NULL}
#' @param \code{n.permutation} an optional integer number, default=10000
#' @param \code{method} an optional string with possible values \code{identity} (no transformation), \code{mN} (multivariate Normal), or \code{mlogN} (multivariate log-Normal), default=\code{mN}
#' @param \code{alpha} an optional probability value, default=0.05
#' @param \code{verbose} an optional logic value to suppress some output status messages, default=TRUE
#'
#' @return \code{matPvals} the matrix of p-values
#' @return \code{q} the vector of minimum p-values
#' @return \code{res} the vector of result estimates which includes:
#' \code{MWSL} = metabolome-wide significance level (MWSL),
#' \code{MWSL_CI.up} = upper value confidence interval MWSL,
#' \code{MWSL_CI.low} = lower value confidence interval MWSL,
#' \code{ENT} = effective number of tests (ENT),
#' \code{ENT_CI.up} = upper value confidence interval ENT,
#' \code{ENT_CI.low} = lower value confidence interval ENT,
#' \code{R.percent} = ENT/M
#' @return \code{t1err.percent} the type I error rate
#'
#' @export
#'
FWERperm2 <- function(outcome,features,confounders=NULL,n.permutation=NULL,method=NULL,alpha=NULL,verbose=TRUE,seed=11032021){
  Y <- outcome
  n <- nrow(features)
  M <- ncol(features)
  K <- ifelse(is.null(n.permutation)==T,10000,n.permutation)
  method <- ifelse(is.null(method)==T,'mN',method)
  alpha <- ifelse(is.null(alpha)==T,0.05,alpha)
  
  ### real X
  if (method=='identity'){Xmat <- features}
  
  ### simulated X ~ mN(0,Sigma*)
  if (method=='mN'){Xmat <- mvrnorm(n=n,mu=rep(0,M),Sigma=cov.shrink(features,verbose=verbose))}
  
  ### simulated X ~ mlogN(0,Sigma*)
  if (method=='mlogN'){Xmat <- mvrnorm(n=n,mu=rep(0,M),Sigma=cov.shrink(apply(features+abs(range(features)[1])+1, 2, function(x) {log(x)}),verbose=verbose))}
  
  if (verbose==T){
    print(paste("n.permutation =", K))
    print(paste("n.samples =", n))
    print(paste("n.features =", M))
    print(paste("method =", method))}
  
  # Y
  if (length(Y)<n) {
    time <- Y[,1]
    status <- Y[,2]
    model_type <- 'Cox'
    if (verbose==T){print(paste("model_type =", model_type))}
  } else {
    outcome_type <- ifelse(length(Y) != length(which(Y==as.integer(Y))) | length(table(Y))>180,'continuous','discrete')
    if (verbose==T){print(paste("outcome_type =", outcome_type))}
    
    if (outcome_type=='continuous') {
      model_type <- 'OLS'
      if (verbose==T){print(paste("model_type =", model_type))}
    }
    if (outcome_type=='discrete') {
      outcome_type.discrete <- ifelse(length(table(Y))>2,'count','binary')
      if (verbose==T){print(paste("outcome_type.discrete =", outcome_type.discrete))}
      model_type <- 'GLM'
      if (verbose==T){print(paste("model_type =", model_type))}
      
      if (outcome_type.discrete=='binary') {family_type <- "binomial"}
      if (outcome_type.discrete=='count') {family_type <- ifelse(var(Y)/mean(Y)>1.5,"negative.binomial","poisson")}
      
      if (verbose==T){print(paste("family_type =", family_type))}
    }
  }
  
  
  if (is.null(confounders)==F){
    
    YZmat <- cbind(Y,confounders)
    
    cl <- makeCluster(parallel::detectCores())
    registerDoParallel(cl)
    
    set.seed(seed)
    YZmat_permuted <-
      foreach(k=1:K) %dorng% {
        rowPerm = sample(nrow(YZmat), replace=F)
        YZmat[rowPerm,]
      }
    
    tic = proc.time()
    Pmatrix <-
      foreach(k=1:K, .combine='cbind') %:%
      foreach(m=1:M, .combine='c') %dopar% {
        new_data <- cbind(YZmat_permuted[[k]],Xmat[,m])
        names(new_data)[names(new_data)=="Xmat[, m]"] <- "X"
        
        if (model_type == 'OLS') {
          fmla <- as.formula(paste("Y ~ X +", paste(colnames(confounders),collapse="+")))
          reg.out <- lm(formula=fmla,data=as.data.frame(new_data))
          pval <- coef(summary(reg.out))["X","Pr(>|t|)"]}
        
        if (model_type == 'GLM') {
          fmla <- as.formula(paste("Y ~ X +", paste(colnames(confounders),collapse="+")))
          if (family_type %in% c('binomial')) {
            reg.out <- glm(fmla,data=as.data.frame(new_data),family="binomial")
            pval <- coef(summary(reg.out))["X","Pr(>|z|)"]}
          if (family_type %in% c('poisson')) {
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=stats::poisson(link="log"))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|z|)"]))}
          if (family_type =='negative.binomial'){
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=MASS::negative.binomial(theta=mean(Y)))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|t|)"]))}}
        
        if (model_type == 'Cox') {
          if (ncol(Y)>2){
            fmla <- as.formula(paste("survival::Surv(time1,time2,status) ~ X +", paste(colnames(confounders),collapse="+")))
          } else {
            fmla <- as.formula(paste("survival::Surv(time,status) ~ X +", paste(colnames(confounders),collapse="+")))}
          reg.out <- survival::coxph(formula=fmla,data=as.data.frame(new_data))
          pval <- summary(reg.out)$coefficients["X","Pr(>|z|)"]}
        pval
      }
    toc = proc.time() - tic
    
    parallel::stopCluster(cl)
  
 }
  
  if (is.null(confounders)==T){
    
    cl <- makeCluster(parallel::detectCores()-1)
    registerDoParallel(cl)
    
    set.seed(seed)
    Y_permuted <-
      foreach(k=1:K) %dorng% {
        sample(Y,replace=F)
      }
    
    Pmatrix <-
      foreach(k=1:K, .combine='cbind') %:%
      foreach(m=1:M, .combine='c') %dopar% {
        new_data <- cbind(Y_permuted[[k]],Xmat[,m])
        
        if (model_type == 'OLS') {
          colnames(new_data) <- c('Y','X')
          fmla <- as.formula(paste("Y ~ X"))
          reg.out <- lm(fmla, data=as.data.frame(new_data))
          pval <- coef(summary(reg.out))["X","Pr(>|t|)"]}
        
        if (model_type == 'GLM') {
          colnames(new_data) <- c('Y','X')
          fmla <- as.formula(paste("Y ~ X"))
          if (family_type %in% c('binomial')) {
            reg.out <- glm(fmla,data=as.data.frame(new_data),family="binomial")
            pval <- coef(summary(reg.out))["X","Pr(>|z|)"]}
          if (family_type %in% c('poisson')) {
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=stats::poisson(link="log"))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|z|)"]))}
          if (family_type =='negative.binomial'){
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=MASS::negative.binomial(theta=mean(Y)))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|t|)"]))}}
        
        if (model_type == 'Cox') {
          if (ncol(Y)>2){
            colnames(new_data) <- c('time1','time2','status','X')
            fmla <- as.formula(paste("survival::Surv('time1','time2','status') ~ X"))
          } else {
            colnames(new_data) <- c('time','status','X')
            fmla <- as.formula(paste("survival::Surv(time, status) ~ X"))
          }
          reg.out <- survival::coxph(formula=fmla,data=as.data.frame(new_data))
          pval <- summary(reg.out)$coefficients["X","Pr(>|z|)"]}
        pval
      }
    
    parallel::stopCluster(cl)
    
  }
  
  q = apply(Pmatrix, 2, min)
  t1err.percent = sum(Pmatrix <= alpha)*100/(M*K)
  
  q <- sort(q)
  ceil = ceiling(alpha*K)
  MWSL <- q[ceil]
  MWSL_CI.up <- q[ceil + sqrt(ceil*(1-alpha))]
  MWSL_CI.low <- q[ceil - sqrt(ceil*(1-alpha))]
  
  ENT <- alpha/MWSL
  ENT_CI.up = alpha/MWSL_CI.low
  ENT_CI.low = alpha/MWSL_CI.up
  
  R.percent = (ENT/M)*100
  
  res.ENT <- cbind(MWSL, MWSL_CI.up, MWSL_CI.low, ENT, ENT_CI.up, ENT_CI.low, R.percent)
  
  #t1err.percent <- (mean(t1err)/K)*100
  
  res <- list()
  class(res) = "FWERperm"
  res$matPvals <- Pmatrix
  res$q <- q
  res$res <- res.ENT
  res$t1err.percent <- t1err.percent
  return(res)
}
