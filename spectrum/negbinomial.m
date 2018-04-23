function ret = negbinomial(k, mu, sigma2)
  if ((sigma2 - mu) < 1e-3)
    sigma2 = mu + 1e-3;
  endif
    
  p = (sigma2 - mu) / sigma2;
  r = mu *  mu / (sigma2 - mu);
  ret = exp(lgamma(k + r) - lgamma(k + 1) - lgamma(r) + r*log(1 - p) + k*log(p));
endfunction