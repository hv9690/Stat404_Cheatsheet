welch <- function(A, B) {
  ss1 = var(A); ss2=var(B); R = ss1/ss2*(length(B)/length(A));
  df = (R/(1+R))^2/(length(A)-1)+1/(R+1)^2/(length(B)-1); df = 1/df;
  TT = (mean(A)-mean(B))/(ss1/length(A)+ss2/length(B))^.5;
  pValue = 2*(1-pt(abs(TT), df))
}