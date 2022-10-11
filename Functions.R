welch <- function(A, B, twot = TRUE) {
  ssA = var(A); ssB=var(B); R = ssA/ssB*(length(B)/length(A));
  df = (R/(1+R))^2/(length(A)-1)+1/(R+1)^2/(length(B)-1); df = 1/df;
  TT = (mean(A)-mean(B))/(ssA/length(A)+ssB/length(B))^.5;
  if (twot == TRUE) {
    pValue = 2*(1-pt(abs(TT), df))
  } else {
    pValue = (1-pt(abs(TT), df))
  }
  print(c(df, TT, pValue))
}

equal_t <- function(A,B, twot = TRUE) {
  ssPool = (length(B)*var(B) + length(A)*var(A))/(length(A)+length(B))
  tValue = (mean(A) - mean(B))/((1/length(A)+1/length(B))*ssPool)^.5
  if (twot == TRUE) {
    pValue = 2*(1 - pt(abs(tValue), (length(A)+length(B) - 2)))
  } else {
    pValue = (1 - pt(abs(tValue), (length(A)+length(B) - 2)))
  }
  print(c(ssPool, tValue, pValue))
}

simulation <- function(A,B,reps = 10000) {
  set.seed(2022)
  zz = c(A,B)
  d.obs = abs(mean(A) - mean(B)); dd = rep(0, reps)
  for(i in 1:reps) { 
    ind = sample(1:length(zz), length(A)); dd[i] = mean(zz[ind]-zz[-ind])}
  pvalue = mean(abs(d.obs) < abs(dd)) + 
    0.5*mean(abs(d.obs) == abs(dd))
  print("d.obs and pvalue are")
  print(c(d.obs, pvalue))
}

ANOVAA <- function(y, a, b, c, d, groups = 4) {
  aabar = mean(a); bbbar = mean(b); 
  ccbar = mean(c); ddbar = mean(d); 
  yybar = mean(y);
  SS.trt = length(a)*((aabar - yybar)^2+(bbbar - yybar)^2 + (ccbar - yybar)^2 + (ddbar - yybar)^2)
  MSS.trt = SS.trt/(groups-1)
  SS.e = sum((a - aabar)^2)+sum((b-bbbar)^2)+sum((c-ccbar)^2) + sum((d-ddbar)^2)
  MSS.e = SS.e/(length(y)-groups)
  SS.tot = sum( (y - mean(y))^2)
  ft = MSS.trt/MSS.e
  p.value = pf(ft, (groups-1), length(y)-groups, lower.tail=F)
  print(ft, p.value)
}

Bonferoni <- function(a,b,c,d, groups = 4, alpha = 0.05, pairs = 6) {
  yy = rbind(a, b, c, d)
  yy.bar = rowMeans(yy)
  AB = yy.bar[1] - yy.bar[2]
  AC = yy.bar[1] - yy.bar[3]
  AD = yy.bar[1] - yy.bar[4]
  BC = yy.bar[2] - yy.bar[3]
  BD = yy.bar[2] - yy.bar[4]
  CD = yy.bar[3] - yy.bar[4]
  mu.diff = c(AB, AC, AD, BC, BD, CD)
  sigmahat = MSS.e^.5
  denominator = sigmahat * (1/length(a) + 1/length(a))^.5
  tt = mu.diff/denominator
  print(tt)
  print(qt(1-alpha/pairs/2, groups*(length(a)-1)))
  
  low.limit = mu.diff - qt(1-alpha/pairs/2, groups*(length(a)-1))*denominator
  upper.limit = mu.diff + qt(1-alpha/pairs/2, groups*(length(a)-1))*denominator
  print(round(cbind(low.limit, upper.limit), 3))
}

Tukey <- function(a,b,c,d, groups = 4, alpha = 0.05, pairs = 6) {
  yy = rbind(a, b, c, d)
  yy.bar = rowMeans(yy)
  AB = yy.bar[1] - yy.bar[2]
  AC = yy.bar[1] - yy.bar[3]
  AD = yy.bar[1] - yy.bar[4]
  BC = yy.bar[2] - yy.bar[3]
  BD = yy.bar[2] - yy.bar[4]
  CD = yy.bar[3] - yy.bar[4]
  mu.diff = c(AB, AC, AD, BC, BD, CD)
  sigmahat = MSS.e^.5
  denominator = sigmahat * (1/length(a) + 1/length(a))^.5
  tt = mu.diff/denominator
  print(tt)
  print(qtukey(1-alpha,groups,length(yy)-groups)/sqrt(2))
  
  low.limit2 = mu.diff - qtukey(1-alpha,groups,length(yy)-groups)/sqrt(2)*denominator
  upper.limit2 = mu.diff + qtukey(1-alpha,groups,length(yy)-groups)/sqrt(2)*denominator
  print(round(cbind(low.limit2, upper.limit2), 3))
}

