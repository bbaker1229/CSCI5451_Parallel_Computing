%%
function [nctr] = rand_ftr(V)
  [m,n] = size(V);
  tot = 0.0;
  nctr = zeros(1,n);
  for j=1:m
    t = rand();
    nctr = nctr + t*V(j,:);
    tot = tot+t;
  end
  nctr = nctr/tot;
  
