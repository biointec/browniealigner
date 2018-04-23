function ret = negbinomialratio(k, mu1, sigma21, mu2, sigma22)
  
  if ((sigma21 - mu1) < 1e-3)
    sigma21 = mu1 + 1e-3;
  endif
  
  if ((sigma22 - mu2) < 1e-3)
    sigma22 = mu2 + 1e-3;
  endif
  
  p1 = (sigma21 - mu1) / sigma21;
  r1 = mu1 *  mu1 / (sigma21 - mu1);
  p2 = (sigma22 - mu2) / sigma22;
  r2 = mu2 *  mu2 / (sigma22 - mu2);
  ret = exp(lgamma(k + r1) - lgamma(r1) + r1*log(1 - p1) + k*log(p1) - lgamma(k + r2) + lgamma(r2) - r2*log(1 - p2) - k*log(p2));
endfunction