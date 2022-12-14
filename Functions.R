lmfake <- function(x, y){
  n = nrow(y)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  rr = y - yhat
  rsq = (1 - sum(rr^2)/((n-1)*var(y)))
  print(rr)
  betahat
}

predict <- function(x,y, newx) {
  betah = lmfake(x,y)
  betah[1] + betah[-1] %*% newx
}

confint <- function(x, y, alpha = 0.05, type = c("conf", "pref")) {
  new.Xvec = c(1, predict(x,y,newx))
  n = nrow(y)
  k = ncol(x) - 1
  betah = lmfake(x,y)
  sigma2 = t(y-x%*%betah)%*%(y-x%*%betah) / (n-k-1)
  cv = qt(alpha/2, n-k-1, lower.tail=FALSE)
  if (type == "conf") {
    se.conf = sqrt(sigma2*(t(new.Xvec) %*% solve(t(x)%*%x) %*% new.Xvec) )
    c(pred-cv*se.conf, pred+cv*se.conf)
  } else {
    se.pred = sqrt(sigma2*(1+t(new.Xvec) %*% solve(t(x)%*%x) %*% new.Xvec) )
    c(pred-cv*se.pred, pred+cv*se.pred)
  }
}

anovareg <- function(x, y) {
  yyhat = x %*% betahat
  yybar = mean(y)
  SSreg = sum((yyhat - yybar)^2)
  SSresid = sum(rr^2)
  SStot = sum((yy - yybar)^2)
  f_stat = SSreg/SSresid
  print(c(SSreg, SSresid, SStot, f_stat))
}

variancereg <- function(x, y, hyp = c(0, -1, 1)) {
  n = length(y)
  k = ncol(x)-1
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  rr = y - yhat
  sig2hat = sum(rr^2) / (n-k-1)
  var.beta = solve(t(x) %*% x) * sig2hat
  t_stat = beta/var.beta
  p_value = 1-pt(t_stat, n-k)
  t(hyp) %*% var.beta %*% hyp
  print(c(t_stat,p_value))
}

welch <- function(A, B, twot = TRUE) {
  ssA = var(A); ssB=var(B); R = ssA/ssB*(length(B)/length(A));
  df = (R/(1+R))^2/(length(A)-1)+(1/(length(B)-1))/((1+R)^2); df = 1/df;
  TT = (mean(A)-mean(B))/(ssA/length(A)+ssB/length(B))^.5;
  if (twot == TRUE) {
    pValue = 2*(1-pt(abs(TT), df))
  } else {
    pValue = (1-pt(abs(TT), df))
  }
  print(c(df, TT, pValue))
}

equal_t <- function(A,B, twot = TRUE) {
  ssPool = ((length(B)-1)*var(B) + (length(A)-1)*var(A))/(length(A)+length(B)-2)
  tValue = (mean(A) - mean(B))/((1/length(A)+1/length(B))*ssPool)^.5
  if (twot == TRUE) {
    pValue = 2*(1 - pt(abs(tValue), (length(A)+length(B) - 2)))
  } else {
    pValue = (1 - pt(abs(tValue), (length(A)+length(B) - 2)))
  }
  print(c(ssPool, tValue, pValue))
}

simulation <- function(A,B,reps = 10000, two.tail = TRUE) {
  set.seed(1)
  zz = c(A,B)
  d.obs = mean(A) - mean(B); dd = rep(0, reps);
  for(i in 1:reps) { 
    ind = sample(length(zz), length(A))
    dd[i] = mean(zz[ind]) - mean(zz[-ind])
  }
  count = (sum(dd>d.obs) + 0.5*sum(d.obs==dd))
  p.temp = count/reps
  if (two.tail == TRUE) {
    pvalue = 2*min(p.temp, 1-p.temp)
  } else {
    pvalue = p.temp
  }
  print("d.obs and pvalue are")
  print(c(d.obs, pvalue))
}
simulation(aa,bb,50000,FALSE)

y = c(yy1, yy2, yy3, yy4, yy5)
trtt =
  as.factor(rep(1:5, c(
    length(yy1),
    length(yy2),
    length(yy3),
    length(yy4),
    length(yy5)
  )))

dd = data.frame(x = y, trt = trtt)

N = length(y)

ANOVAA <- function(trt = trtt, x = y, data = dd) {
  y.bar = mean(y)
  yi.bar = tapply(x, INDEX=factor(trt), FUN=mean)
  N = nrow(data)
  ni = summary(dd$trt)
  k = length(ni)
  xxbar = mean(x)
  SS.tot = sum(x**2) - N*(xxbar)**2
  SS.trt = sum(ni*(yi.bar)**2) - N*(y.bar)**2
  SS.err = SS.tot - SS.trt
  MS.trt = SS.trt / (k-1)
  MS.err = signmahat = SS.err / (N-k)
  ft = MS.trt/MS.err
  p.value = pf(ft, (k-1), N-k, lower.tail=F)
  print(c(k-1, SS.trt, MS.trt, ft))
  print(c(N-k, SS.err, MS.err))
  print(c(N-1, SS.tot))
  print(p.value)
  signmahat
}

Bonferoni <- function(trt, MS.err = ANOVAA(), datap = 64, alpha = 0.05) {
  yi.bar = tapply(y, INDEX=factor(trt), FUN=mean)
  groups = levels(trt)
  k = length(groups)
  pairs = combn(k, 2)
  name1 = groups[pairs[1,]]
  name2 = groups[pairs[2,]]
  pair.names = paste(name1, name2, sep="vs")
  m = choose(k, 2)
  diff.means = apply(pairs, MARGIN=2, function(x) -diff(yi.bar[x]))
  se = sqrt(MS.err*(1/(datap/k)+1/(datap/k))); se = rep(se, m)
  cv.bonferroni = qt(alpha/(m*2), N-k, lower.tail=FALSE)
  LB.Bonf = diff.means - se*cv.bonferroni
  UB.Bonf = diff.means + se*cv.bonferroni
  pairwise.res = data.frame(diff.means, LB.Bonf, UB.Bonf)
  row.names(pairwise.res) = pair.names
  print(round(pairwise.res, 4))
}

Tukey <- function(trt = trtt, MS.err = ANOVAA(), alpha = 0.05) {
  yi.bar = tapply(y, INDEX=factor(trt), FUN=mean)
  groups = levels(trt)
  k = length(groups)
  pairs = combn(k, 2)
  name1 = groups[pairs[1,]]
  name2 = groups[pairs[2,]]
  pair.names = paste(name1, name2, sep="vs")
  m = choose(k, 2)
  diff.means = apply(pairs, MARGIN=2, function(x) -diff(yi.bar[x]))
  se = sqrt(MS.err*c(1/ni[1] + 1/ni[2],
                     1/ni[1] + 1/ni[3],
                     1/ni[1] + 1/ni[4],
                     1/ni[1] + 1/ni[5],
                     1/ni[2] + 1/ni[3],
                     1/ni[2] + 1/ni[4],
                     1/ni[2] + 1/ni[5],
                     1/ni[3] + 1/ni[4],
                     1/ni[3] + 1/ni[5],
                     1/ni[4] + 1/ni[5]))
  cv.tukey = qtukey(1 - alpha,k,N-k) / sqrt(2)
  LB.Tukey = diff.means - se*cv.tukey
  UB.Tukey = diff.means + se*cv.tukey
  pairwise.res <- data.frame(diff.means, LB.Tukey, UB.Tukey)
  row.names(pairwise.res) = pair.names
  print(round(pairwise.res, 3))
}

powerTest <- function(n, k) {
  tao <- yi.bar - y.bar
  bar.tau = sum(ni*tao)/n
  kk = k-1
  delta = sum(ni*(tao-bar.tau)^2) / (10*MS.err)
  qq = qf(0.95, kk, n-k); 
  print(c("critical value=", round(qq,4)))
  power = pf(qq, kk, n-k, delta, lower.tail=F); 
  print(c("power=",round(power, 4)))
}

powerSim <- function(alpha, k = 4, nn = 1:100, delta, power = 0.8) {
  qq = qf(1-alpha, k-1, nn*k-k)   # Critical value of central F
  pwrs = sapply(seq(length(nn)),  # Calculate power for each N=n*k in vector
                function(i) pf(qq[i],k-1,nn[i]*k-k,delta[i],lower.tail=FALSE))
  min_n_ind = which(pwrs >= power)[1]  # Find index of minimum n needed
  per_grp = nn[min_n_ind]  # Minimum n needed
  total = nn[min_n_ind] * k
}


eta <- function(n = c(5,5,5,5), alpha) {
  k = length(n)
  e = sum(yy)/sum(n)
  signmadsq = (MS.trt-MS.err)/n[1]
  vare = signmadsq/n[1] + MS.err/sum(n)
  conf = c(e-qt(1-alpha/2,k)*sqrt(vare), e+qt(1-alpha/2, k)*sqrt(vare))
  print(c(e, conf))
}

