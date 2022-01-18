## "Homemade" GEE function (jebsgee) for continuous outcomes solved using IRLS
## with an exchangeable (ex) or independence (in) working correlation matrix
## Example showing how generalized estimating equations work...
## AUTHOR: Francis L. Huang, PhD
## University of Missouri, huangf@missouri.edu

## To accompany: Huang, F. (2021). Analyzing cross-sectionally clustered
## data using generalized estimating equations. Journal of Educational
## and Behavioral Statistics.

## 2020.11.07: allows for unsorted data, can match corstr
## 2020.11.09: tidying up
## 2021.04.06: added stars to the output

## Two helper functions included: 
## 1) exchR: to compute the exchangeable working correlation matrix
## 2) geticc: to compute the ICC based on residuals
## main function is jebsgee
## Example usage is included at the end of the syntax

## only for creating an exchangeable R matrix
exchR <- function(icc, maxsize, y, NG, gpsz, ...){
  wr1 <- matrix(icc, nrow = maxsize, ncol = maxsize)
  diag(wr1) <- 1
  ## converting to a covariance matrix
  wcv <- wr1 * var(y) #save it, used when getting the RSE
  wcl <- list() #empty list
  for (i in 1:NG){ #making several covariance matrices
    GS <- gpsz[i]  #depending on how many units per cluster
    tmp <- wcv[1:GS, 1:GS]
    wcl[[i]] <- tmp
  }
  vm2 <- solve(Matrix::bdiag(wcl)) #create block diagonal (this is V^-1)
  return(list(vm2 = vm2, wcv = wcv)) #return the inverse of the variance matrix
}

geticc <- function(data, cluster, r){
  
  scalep <- mean(r^2) #dispersion parameter
  rpg <- split(r, data[,cluster]) #individual residuals
  nm <- names(table(data[,cluster])) #names of clusters
  coll <- numeric() #empty container
  
  ### p. 63 Hardin and Hilbe
  ### get ICC based on residuals per cluster
  ### right now, need the data to be sorted by cluster
  
  multresid <- function(x){
    r2 <- rpg[[x]] #extract resid per group (rpg)
    egeg <- r2 %*% t(r2) #e %*% t(e)
    coll[x] <- sum(egeg[lower.tri(egeg)]) #only lower diag
  }
  
  tst <- sum(sapply(nm, multresid)) #how many per group 
  ns <- sapply(rpg, length) #add up how many were added
  den <- sum((ns * (ns - 1 )) / 2) #how many products were added
  icc.model <- (tst / den) * (1 / scalep) #the icc
  return(icc.model)
  
}

# HOMEGROWN GEE FUNCTION using iterative reweighted least squares (IRLS)
jebsgee <- function(fml, data, cluster, corstr = 'independence'){
  ## extract data
  tmp <- cbind(data, cluster = data[,cluster]) #dataframe with cluster
  tmp <- tmp[order(tmp$cluster), ] #sorting by cluster
  fml <- formula(fml)
  df <- model.frame(fml, tmp)
  X <- model.matrix(fml, df)
  y <- model.response(df)

  if(sum(is.na(df)) > 0) (stop("You have missing data."))
  
  gpsz <- table(data[, cluster]) #how many per group; group size
  NG <- length(gpsz) #how many groups
  maxsize <- max(gpsz) #what's the biggest group size
  
  CS <- c('independence', 'exchangeable')
  cs <- pmatch(corstr, CS, -1) #allow to match corstr by keywords
  if (cs == -1) (stop("Currently can only use an independence or exchangeable correlation structure"))
  
  corstr <- CS[cs] #put in the whole word
  
  # STEP #1
  firstrun <- glm(formula(fml), data = df) #just a regular glm
  
  # STEP #2
  r <- resid(firstrun, 'pearson') #get residuals
  betas <- coef(firstrun) #get initial coefficients
  
  if (corstr == 'exchangeable') {

  ### setup iterations need for exchangeable structure
    dev <- 0
    delta.dev <- 1
    tol <- 1e-5 #can make this bigger or smaller
    maxiter <- 50 #number of iterations, can make this bigger
    i = 1 #starting at iteration 1
    
    cat("Iteration: ") 
    while(abs(delta.dev > tol & i < maxiter)){ #when change in deviance is small, stop
      cat(i, "::")
      
      # after iteration, this is STEP #5
      r <- y - X %*% betas #residuals / use Pearson if non-identity link
      # STEP #3
      icc <- geticc(data = tmp, 'cluster', r) #compute new iccs
      results <- exchR(icc, maxsize = maxsize, y = y, NG, gpsz) #compute new weight matrix
      vm2 <- results$vm2
      # STEP #4
      betas <- solve(t(X) %*% vm2 %*% X) %*% t(X) %*% vm2 %*% y #update betas
      dev0 <- dev #get prior dev
      dev <- sum((y - X %*% betas)^2) #new deviance
      delta.dev <- dev - dev0 #change in deviance
      i = i + 1 #add one to the iteration
    }
    
    cat("\nFinal alpha:", icc, "\n")
  }
  
  ### STEP #7: computing Liang and Zeger SEs
  re <- as.numeric(y - X %*% betas) #get residuals
  k <- ncol(X) #how many predictors (including intercept)
  cdata <- data.frame(cluster = tmp$cluster, r = re) #data with cluster and residuals
  gs <- names(table(cdata$cluster)) #names of the clusters

  u <- matrix(NA, nrow = NG, ncol = k) #empty matrix
  gpsv <- tmp$cluster
  
  if (corstr == 'independence') wcv <- vm2 <- diag(nrow(X)) #if independence
  if (corstr == 'exchangeable') wcv <- results$wcv
  
  for(i in 1:NG){
      tmp <- nrow((df[gpsv == gs[i], ]))
      u[i,] <- t(cdata$r[cdata$cluster == gs[i]]) %*% 
        solve(wcv[1:tmp, 1:tmp]) %*% X[gpsv == gs[i], 1:k]
  }
  
  mt <- crossprod(u) #t(u) %*% u :: meat
  br <- solve(t(X) %*% vm2 %*% X) #bread matrix
  clvc <- br %*% mt %*% br #LZ robust vcov matrix
  
  ### putting it all together
  
  se <- as.numeric(sqrt(diag(clvc))) #standard error
  b <- as.numeric(betas) #betas
  wald <- (b / se)^2 #wald
  pv <- pchisq(wald, 1, lower.tail = F) #p value
  stars <- cut(pv, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
               labels = c("***", "**", "*  ", ".  ", " "), 
               include.lowest = TRUE)
  
  res <- data.frame(estimates = b, se, wald, pv = round(pv, 4), s = stars)
  row.names(res) <- colnames(X) #getting the names of the coefficients
  
  cat("Working correlation structure:", corstr, "\n")
  print(res) #output results
}

### TESTING with some toy examples
### Compare with results using geeglm: they are the same

library(geepack) #for geeglm function

## Examples 1 and 2
## using mtcars, cyl is the clustering variable
### COMPARE RESULTS

mtcars2 <- mtcars[order(mtcars$cyl), ] #needs to be sorted for geeglm
jebsgee(mpg ~ wt + am + qsec + hp + vs, 
        data = mtcars, 
        cluster = 'cyl', 
        corstr = 'ex')
summary(geeglm(mpg ~ wt + am + qsec + hp + vs, 
        id = cyl, 
        corstr = 'ex', 
        data = mtcars2)
)

## Using an independence working correlation matrix structure
jebsgee(mpg ~ wt + am + qsec + hp + vs, 
        data = mtcars, 
        cluster = 'cyl', 
        corstr = 'in')
summary(geeglm(mpg ~ wt + am + qsec + hp + vs, 
        id = cyl, 
        corstr = 'ind', 
        data = mtcars2)) #with an independence working R matrix

### Example 3: using High School and Beyond
library(mlmRev)
data(Hsb82)
jebsgee(mAch ~ sector + meanses + cses + sector * cses, 
        data = Hsb82, 
        cluster = 'school',
        corstr = 'ex')
summary(geeglm(mAch ~ sector + meanses + cses + sector * cses, 
        id = school, 
        corstr = 'ex', 
        data = Hsb82))

### END:: 2020.11.11
