
#' Linear regression models using an S3 class and the calculations are done using QR Decomposition
#'
#' @param formula formula
#' @param data a dataframe
#'
#' @return Returns an object of the class linreg
#' @export
#'
#' @examples  linreg_mod <-linreg(Petal.Length~Sepal.Width+Sepal.Length, data=iris)
#'  summary(linreg_mod)
#'
linreg <- function(formula, data){
  Z <- all.vars(formula)
  for(i in 1:length(Z)){
    if(!(Z[i] %in% names(data))) stop("wrong column name!")
  }

  X <- model.matrix(formula,data)

  y <- Z[1]

  cl <- match.call()

  lnrg <- list()
  class(lnrg) <- "linreg"

  new_qr_x <- qr(X)
  lnrg$regcoe <- solve(qr.R(new_qr_x)) %*% t(qr.Q(new_qr_x)) %*% data[,y]
  lnrg$fitval <- qr.fitted(new_qr_x,data[,y])
  lnrg$residu <- qr.resid(new_qr_x,data[,y])
  lnrg$degfre <- nrow(data)-ncol(X)
  lnrg$resvar <- (t(lnrg$residu) %*% lnrg$residu) /lnrg$degfre
  sca_resvar <- lnrg$resvar[1,1]

  regcoemat1 <- chol2inv(qr.R(new_qr_x))
  regcoemat2 <- vector(length = nrow(regcoemat1))
  for (i in 1:nrow(regcoemat1)) {
    regcoemat2[i] <- regcoemat1[i,i]
  }

  lnrg$varregcoe <- regcoemat2 * sca_resvar

  tvalue <- vector(length = nrow(lnrg$regcoe))
  for (i in 1:nrow(lnrg$regcoe)) {
    tvalue[i] <- lnrg$regcoe[i,1]/sqrt(lnrg$varregcoe[i])
  }
  lnrg$tval <- tvalue
  lnrg$pval <- pt(abs(lnrg$tval),lnrg$degfre, lower.tail = FALSE)
  lnrg$call <- cl

  return(lnrg)
}


resid<- function(x){
  return(x$residu)
}

pred<- function(x){
  return(x$fitval)
}
coef.linreg <- function(x){
  named_coe <- vector(length = length(x$regcoe))
  vector_name <- vector(length = length(x$regcoe))
  for (i in 1:nrow(x$regcoe)) {
    named_coe[i] <- x$regcoe[i,1]
    vector_name[i] <- names(x$regcoe[i,1])
  }
  names(named_coe) <- vector_name
  return(named_coe)
}

summary.linreg <- function(x){
  if(length(coef(x))){
    cat("Coefficients:\n")

    # the output has 5 columns, the first four are numeric, the fifth is character
    #  change summ from matrix to data.frame
    #  the fifth column has no name, so its name is placed by ""
    stderror <- sqrt(x$varregcoe)
    summ <- data.frame(x$regcoe,stderror,x$tval,x$pval, c("***"),stringsAsFactors = F)
    colnames(summ) <- c("Estimate", "Std. Error", "t value", "P(>|t|)","")
    print(summ)
    cat("\n\nResidual standard error:", sqrt(x$resvar),"on",x$degfre, "degrees of freedom" )
  }
  else cat("No coefficients\n")
}


plot.linreg <- function(x){
  data1 <- cbind(x$fitval,x$residu)
  data1 <- as.data.frame(data1)
  names(data1) <- c("fitval","residu")
  p <- ggplot(data1, aes(x=fitval, y=residu)) + geom_point(shape=1)
  p <- p + labs(x="Fitted values", y="Residuals") + ggtitle("Residuals vs Fitted")
  p <- p + geom_smooth(method = lm)




  sca_resvar <- x$resvar[1,1]
  st <- sqrt(sca_resvar)
  abs_residu <- abs(x$residu)/st
  stadarres <- sqrt(abs_residu)
  data2 <- cbind(x$fitval,stadarres)
  data2 <- as.data.frame(data2)
  names(data2) <- c("fitval","stadarres")
  p2 <- ggplot(data2, aes(x=fitval, y=stadarres)) + geom_point(shape=1)
  p2 <- p2+ labs(x="Fitted values", y="|Standardized residuals|") + ggtitle("Scale-Location")

  return(p)



}

print.linreg <- function(x, digits=max(3,getOption("digits")-3)){
  cat("\nCall:\n", deparse(x$call),"\n\n", sep = "")
  if(length(coef(x))){
    cat("Coefficients:\n")
    print.default(format(coef(x),digits = digits), print.gap = 2,quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)

}
