function ret = mixtureModel(spectrum, C, maxIt)

  numElements = size(spectrum,1);
  numComponents = size(C,1)
  
  weights = zeros(numElements, numComponents);
 
  # mixture model fitting
  for itCount = 1:maxIt

    # compute the weights 
    for j = 1:numComponents
      nom = ones(numElements,1);
      for k = 1:numComponents
        if (k != j)
          nom += C(k, 3)/C(j,3) * negbinomialratio(spectrum(:,1), C(k,1), C(k,2), C(j,1), C(j,2));
        endif          
      endfor
      weights(:,j) = 1./nom;
    endfor

    # compute the mean
    for j = 1:numComponents
      count = sum(weights(:,j) .* spectrum(:,2));
      total = sum(weights(:,j) .* spectrum(:,1) .* spectrum(:,2));

      C(j,1) = total / count;
      C(j,3) = count;
    endfor

    # compute the variance
    for j = 1:numComponents
      
      total = sum(weights(:,j) .* spectrum(:,2) .* (spectrum(:,1) - C(j,1)) .* (spectrum(:,1) - C(j,1)));
      count = sum(weights(:,j) .* spectrum(:,2));
     
      C(j,2) = total / count;

      if (C(j,2) - C(j,1) < 1e-3)
        C(j,2) = C(j,1) + 1e-3;
      endif
    endfor  

  endfor # of main loop
  
  ret = C;
endfunction