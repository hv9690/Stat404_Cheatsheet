lmfake <- function(x, y){
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  print(betahat)
  yhat = x%*%betahat
  rr = y - yhat
  print(rr)
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

variancereg <- function(x, y, theta, hyp = c(0, 0, -2, 1)) {
  n = length(y)
  sig2hat = sum(rr^2) / (n - length(hyp))
  var.beta = solve(t(x) %*% x) * sig2hat
  vv = as.matrix(hyp)
  var.theta = t(vv) %*% var.beta %*% vv
  t_stat = theta/var.theta
  p_value = 1-pt(t_stat, n-length(hyp))
  print(c(var.theta,t_stat,p_value))
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
  set.seed(2022)
  zz = c(A,B)
  d.obs = mean(A) - mean(B); count = 0;
  for(i in 1:reps) { 
    ind = sample(1:(length(A)+length(B)), length(A))
    TT = mean(zz[ind]) - mean(zz[-ind])
    if (TT < d.obs) {
      count = count + 1
    } else if (TT == d.obs) {
      count = count + 0.5
    }
  }
  p.temp = count/reps
  if (two.tail == TRUE) {
    pvalue = 2*min(p.temp, 1-p.temp)
  } else {
    pvalue = min(p.temp, 1-p.temp)
  }
  print("d.obs and pvalue are")
  print(c(d.obs, pvalue))
}

ANOVAA <- function(trt, x, data) {
  N = nrow(data)
  k = length(levels(trt))
  xxbar = mean(x)
  SS.tot = sum(x**2) - N*(xxbar)**2
  yi.bar = tapply(x, INDEX=factor(trt), FUN=mean)
  ni = table(y);
  SS.trt = sum(ni*(yi.bar)**2) - N*(grandMean)**2
  SS.err = SS.tot - SS.trt
  MS.trt = SS.trt / (k-1)
  MS.err = SS.err / (N-k)
  ft = MS.trt/MS.err
  p.value = pf(ft, (k-1), N-k, lower.tail=F)
  print(c(k-1, SS.trt, MS.trt, ft))
  print(c(N-k, SS.err, MS.err))
  print(c(N-1, SS.tot))
}

Bonferoni <- function(trt, yi.bar, datap = 64, alpha = 0.05) {
  groups = levels(trt)
  k = length(groups)
  pairs = combn(k, 2)
  name1 = groups[pairs[1,]]
  name2 = groups[pairs[2,]]
  pair.names = paste(name1, name2, sep="vs")
  m = choose(k, 2)
  diff.means = apply(pairs, MARGIN=2, function(x) +diff(yi.bar[x]))
  se = sqrt(MS.err*(1/(datap/k)+1/(datap/k))); se = rep(se, m)
  cv.bonferroni = qt(alpha/(m*2), N-k, lower.tail=FALSE)
  LB.Bonf = diff.means - se*cv.bonferroni
  UB.Bonf = diff.means + se*cv.bonferroni
  pairwise.res = data.frame(diff.means, LB.Bonf, UB.Bonf)
  row.names(pairwise.res) = pair.names
  print(round(pairwise.res, 3))
}

Tukey <- function(trt, yi.bar, datap = 64, alpha = 0.05) {
  groups = levels(trt)
  k = length(groups)
  pairs = combn(k, 2)
  name1 = groups[pairs[1,]]
  name2 = groups[pairs[2,]]
  pair.names = paste(name1, name2, sep="vs")
  m = choose(k, 2)
  diff.means = apply(pairs, MARGIN=2, function(x) +diff(yi.bar[x]))
  se = sqrt(MS.err*(1/(datap/k)+1/(datap/k))); se = rep(se, m)
  cv.tukey = qtukey(alpha,k,datap-k) / sqrt(2)
  LB.Tukey = diff.means - se*cv.tukey
  UB.Tukey = diff.means + se*cv.tukey
  pairwise.res <- data.frame(diff.means, LB.Tukey, UB.Tukey)
  row.names(pairwise.res) = pair.names
  print(round(pairwise.res, 3))
}

powerTest <- function(n, k, delta) {
  qq = round(qf(0.95, k-1, n*k-k), 4); 
  print(c("critical value=", qq))
  power = pf(qq, k-1, n*k-k, delta, lower.tail=F); 
  print(c("power=",round(power, 4)))
}

delt <- function(n = c(5,5,5,5), tao, sigma) {
  delta = sum(n*(tao-mean(tao)**2)/sigma**2)
  signmadsq = (MSS.trt-MSS.e)/n[1]
  print(c(delta, sigmadsq))
}

eta <- function(n = c(5,5,5,5), alpha) {
  e = sum(yy)/sum(n)
  vare = signmadsq/n[1] + MSS.e/sum(n)
  conf = c(e-qt(alpha, n[1]-1)*sqrt(vare), e+qt(alpha, n[1]-1)*sqrt(vare))
  print(c(e, conf))
}

