#' ---
#' title: "A robust model for bayesian material value estimation"
#' author: "Shimamura _et al_."
#' 
#' output: 
#'   html_document:
#'      toc: true # table of content true
#'      depth: 3  # upto three depths of headings (specified by #, ## and ###)
#'      number_sections: true  ## if you want number sections at each table header
#'      theme: united  # many options for theme, this one is my favorite.
#'      highlight: tango  # specifies the syntax highlighting style
#'      toc_float: true
#'      self_contained: false
#'      keep_md: true
#'      
#' ---

library(ggplot2)
source("./tools.r")

#' *******************************************************************************************
#' # Experiment

#' ## Extract statistical info. 
#' Load CD and MD  results and 
# /* tate == CD == "20140128132334_crv.csv"
# tate == CD == "20140203110917_crv.csv" ??? */
# expt.CD.filename <- choose.files(caption = 'Select CD(==Tate) "*_crv.csv" file')
(expt.CD.filename <- "./inputdata/heimen-20140128132334_crv.csv")
# /* yoko == MD == "20140204094221_crv.csv" */
# expt.MD.filename <- choose.files(caption = 'Select MD(==Yoko) "*_crv.csv" file')
(expt.MD.filename <- "./inputdata/heimen-20140204094221_crv.csv")

expt.CD.all <- parse_crv(expt.CD.filename)
expt.MD.all <- parse_crv(expt.MD.filename)
expt.CD.all$No <- as.factor(expt.CD.all$No)
expt.MD.all$No <- as.factor(expt.MD.all$No)
pairs(expt.CD.all[,c("Hennikei", "Kajuu", "Idouryou")], pch = ".", main="CD")
pairs(expt.MD.all[,c("Hennikei", "Kajuu", "Idouryou")], pch = ".", main="MD")

#' Trimmig
num.expt = 10 # No. 1 to (num.expt) experiments result are used
expt.CD = crop(part(expt.CD.all, num.expt), "Idouryou", 1, 2)
expt.MD = crop(part(expt.MD.all, num.expt), "Idouryou", 1, 2)
expt.CD$No <- factor(expt.CD$No)
expt.MD$No <- factor(expt.MD$No)
pairs(expt.CD[,c("Hennikei", "Kajuu", "Idouryou")], pch = ".", main="CD (trimmed)")
pairs(expt.MD[,c("Hennikei", "Kajuu", "Idouryou")], pch = ".", main="MD (trimmed)")

coefs.CD <- t(sapply(1:num.expt, function(no) lm(Kajuu ~ Hennikei, data=expt.CD[expt.CD$No == no, ])$coefficients))
coefs.MD <- t(sapply(1:num.expt, function(no) lm(Kajuu ~ Hennikei, data=expt.MD[expt.MD$No == no, ])$coefficients))
coefs <- rbind(cbind(data.frame(const='CD'), data.frame(coefs.CD)), cbind(data.frame(const='MD'), data.frame(coefs.MD)))
colnames(coefs) <- c('Const', 'Intercept', 'Slope')

ggplot(coefs, aes(Intercept, Slope, colour=Const)) + geom_point()

#' *****************************************************************************************
#' # Simulation

# filename.tate <- choose.files(caption = 'Select a "zentansaResult_tate*.txt" file')
(filename.tate <- paste(c("./inputdata", "zentansaResult_tate24.txt"), collapse = '/'))
# filename.yoko <- choose.files(caption = 'Select a "zentansaResult_yoko*.txt" file')
(filename.yoko <- paste(c("./inputdata", "zentansaResult_yoko24.txt"), collapse = '/'))
result.ANSYS.tate <- read.csv(filename.tate, sep=" ", header=FALSE, col.names=c("E1", "E2", "U"))
result.ANSYS.yoko <- read.csv(filename.yoko, sep=" ", header=FALSE, col.names=c("E1", "E2", "U"))
result.ANSYS.tate$U <- -result.ANSYS.tate$U
result.ANSYS.yoko$U <- -result.ANSYS.yoko$U

pairs(result.ANSYS.tate)
pairs(result.ANSYS.yoko)

(N.tate <- nrow(result.ANSYS.tate))
(N.yoko <- nrow(result.ANSYS.yoko))

force = 20
# /* Y.tate <- c(1.50312,1.5008,1.49858,1.49641,1.49836,1.50323,1.50776,1.50772,1.51209,1.51558) */
(Y.tate <- force/coefs[coefs=='CD', 'Slope'])
# /* Y.yoko <- c(1.63114,1.63824,1.61835,1.61404,1.61802,1.61603,1.61552,1.61263,1.6158,1.61752) */
(Y.yoko <- force/coefs[coefs=='MD', 'Slope'])

(L.tate <- length(Y.tate))
(L.yoko <- length(Y.yoko))
EEU.tate <- as.matrix(result.ANSYS.tate)
EEU.yoko <- as.matrix(result.ANSYS.yoko)

fit.lm.tate <- lm(U ~ I(log(E1)^2) + I(log(E2)^2) + log(E1)*log(E2) + log(E1) + log(E2) + 1, data=result.ANSYS.tate)
fit.lm.yoko <- lm(U ~ I(log(E1)^2) + I(log(E2)^2) + log(E1)*log(E2) + log(E1) + log(E2) + 1, data=result.ANSYS.yoko)
(C.tate <- coefficients(fit.lm.tate))
(C.yoko <- coefficients(fit.lm.yoko))

# Visualize
response_surfaces = function(e1, e2 ,coefs)
  coefs[[1]]*1 + coefs[[2]]*I(log(e1)^2) + coefs[[3]]*I(log(e2)^2) + coefs[[4]]*log(e1) + coefs[[5]]*log(e2) + coefs[[6]]*log(e1)*log(e2)
library('functional')
response_surface.tate = Curry(response_surfaces, coefs=fit.lm.tate$coefficients)
response_surface.yoko = Curry(response_surfaces, coefs=fit.lm.yoko$coefficients)
#From: http://stackoverflow.com/questions/18147595/plot-3d-plane-true-regression-surface
my_surface <- function(f, n=10, ...) {
  ranges <- rgl:::.getRanges()
  x <- seq(ranges$xlim[1], ranges$xlim[2], length=n)
  y <- seq(ranges$ylim[1], ranges$ylim[2], length=n)
  z <- outer(x,y,f)
  surface3d(x, y, z, ...)
}
## library(rgl)
## knitr::knit_hooks$set(rgl = hook_rgl)
#+ rgl=TRUE
## plot3d(EEU.tate[, 'E1'], EEU.tate[, 'E2'], EEU.tate[, 'U'], type="p", col="black", xlab="e1", ylab="e2", zlab="u", site=5, lwd=15)
## my_surface(response_surface.tate, alpha=.2 )
## par3d(zoom = 0.5)
## #+ rgl=TRUE
## plot3d(EEU.yoko[, 'E1'], EEU.yoko[, 'E2'], EEU.yoko[, 'U'], type="p", col="black", xlab="e1", ylab="e2", zlab="u", site=5, lwd=15)
## # plot3d(yoko$e1, yoko$e2, yoko$u, type="p", col="blue", xlab="e1", ylab="e2", zlab="u", site=5, lwd=15)
## my_surface(response_surface.yoko, alpha=.2 )
## par3d(zoom = 0.5)
#-

e1.bounds <- c(max(min(EEU.tate[, 1]), min(EEU.yoko[, 1])), min(max(EEU.tate[, 1]), max(EEU.yoko[, 1])))
e2.bounds <- c(max(min(EEU.tate[, 2]), min(EEU.yoko[, 2])), min(max(EEU.tate[, 2]), max(EEU.yoko[, 2])))
ut.bounds <- range(EEU.tate[, 3])
uy.bounds <- range(EEU.yoko[, 3])
(Bounds <- rbind(e1.bounds, e2.bounds, ut.bounds, uy.bounds))

data.combi <- list(Nt=N.tate, Ny=N.yoko, EEUt=EEU.tate, EEUy=EEU.yoko, Lt=L.tate, Ly=L.yoko, Yt=Y.tate, Yy=Y.yoko, Ct=C.tate, Cy=C.yoko, Bounds=Bounds)

#' **********************************************************************************
#' # The Simplest Rubust Model
#+stan, cache=T

library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

stanmodel <- stan_model(file='./model/robust_combi01.stan')
pars.show = c("sigma", "lambda", "e1", "e2", "ut", "uy", "y_pred_t", "y_pred_y", "lp__")
pars.output = c(pars.show, c("log_lik")) # log lik. is used to compute ICs later.
(iter <- 1e3) # Length of each marcov chain. Burn-in is the half of the length by default.
fit.combi <- sampling(stanmodel, data=data.combi, pars=pars.output, seed=1234, iter=iter, chains=4)

#-

library(ggmcmc)
ggs_pairs(ggs(fit.combi), lower = list(continuous = "density"), family=paste(pars.show, collapse="|"))


#' ## visually diagnose MCMC computation

#' Taraces of each chains
invisible(lapply(pars.show, function(para) print(ggs_traceplot(ggs(fit.combi, inc_warmup=TRUE, stan_include_auxiliar=TRUE), family=para))))

#' Comparison of the whole chain with the latest 10 % part 
invisible(lapply(pars.show, function(para) print(ggs_density(ggs(fit.combi, inc_warmup=TRUE, stan_include_auxiliar=TRUE), family=para))))

#' Autocorrelation. Closer to 0, better.
invisible(lapply(pars.show, function(para) print(ggs_autocorrelation(ggs(fit.combi, inc_warmup=TRUE, stan_include_auxiliar=TRUE), family=para))))


# shell.exec2("./output/fit-traceplot.pdf")
# ggmcmc(ggs(fit.combi), file='./output/fit-ggmcmc-combi.pdf', family=paste(pars.show, collapse="|"))
# shell.exec2("./output/fit-ggmcmc.pdf")

ms.combi <- rstan::extract(fit.combi)

#' ## Posterior predictive check
#' Check the fittness of the predictive density with data.  
#' Note these graphs can NOT used to detect overfitting, only underfitting could indicate. 
ggplot(data.frame(ms.combi$y_pred_t), aes(x=ms.combi.y_pred_t)) + geom_line(stat="density") + expand_limits(y=0) + geom_histogram(data = data.frame(Y.tate), aes(Y.tate), fill="cornsilk", colour="grey60", size=.2) + ggtitle("CD")  + xlab("y") + ylab("Obs / Pred")
ggplot(data.frame(ms.combi$y_pred_y), aes(x=ms.combi.y_pred_y)) + geom_line(stat="density") + expand_limits(y=0) + geom_histogram(data = data.frame(Y.tate), aes(Y.yoko), fill="cornsilk", colour="grey60", size=.2) + ggtitle("MD") + xlab("y") + ylab("Obs / Pred")

#-

#' **************************************************************************  
#' ## Information Criteria 
#' In contrast to the Posterior vs observation check, various information criteria such as WAIC and PSIS,
#' indicate overfitting.  
#' See [[Matsuura's homepage]](http://statmodeling.hatenablog.com/entry/comparison-of-LOOCV-and-WAIC) for explanations.  
#'
#' PSIS [Vehtari _et al_. 16] is 
loo::loo(ms.combi$log_lik)$looic/(2*(L.tate+L.yoko))
#' And WAIC [Watanabe10] is
loo::waic(ms.combi$log_lik)$waic/(2*(L.tate+L.yoko))

quantile(ms.combi$e1, probs = c(0.05, 0.95))
quantile(ms.combi$e2, probs = c(0.05, 0.95))

#' *****************************************************************************
#' # References
#' * Papers
#'     + [Watanabe 10]Watanabe, Sumio. (2010). "Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory." Journal of Machine Learning Research 11, 3571–3594.
#'     + [Vehtari _et al_. 16]Vehtari, Aki, Andrew Gelman, and Jonah Gabry. "Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC." Statistics and Computing (2016): 1-20.
#'

