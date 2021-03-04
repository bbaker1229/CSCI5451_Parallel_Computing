  function [ctrs,idx] = myKmeans1(V,Nc,max_its)
% function [ctrs,idx] = myKmeans1(V,Nc,max_its)
%% K-means algorithm
%% IN:     
%% V       = V(m,n) = data matrix with m samples and n features
%% Nc      = number of clusters
%% max_its = max number of iterations allowed  
%% OUT:
%% ctrs = ctrs(Nc,n) centroids of  clisters
%% idx  = idx(m,1) = cluster labels of samples
  [m, n] = size(V);
  tol = 1.e-06;
%% m == number of samples, n = num pf features.
%%-------------------- initialize centers - random means
 for i=1:Nc
   ctrs(i,:) = rand_ftr(V);
 end
%%--------------------Main loop
 idx = zeros(m,1);
 for i=1:max_its
%%-------------------- getClosestCtrs
  counter = zeros(1,Nc);
  pt_sums = zeros(Nc,n);
  for i=1:m
    vt = V(i,:);
    for j=1:Nc
      dists(j) = sum((ctrs(j,:) -vt) .^2);
    end
%%-------------------- closest ctr to V(i,:) 
    [~, k] = min(dists);
    idx(i) = k;
    counter(k)=counter(k)+1;
    pt_sums(k,:) = pt_sums(k,:) + V(i,:);
  end
%%-------------------- second half: reset centroids
  not_conv = 0;          %% if Not_conv indicates not converged yet/.
%%-------------------- 
  for i=1:Nc
    if (counter(i)>0)
      new_c = pt_sums(i,:)/counter(i);
%%-------------------- this is a convergence test
      if (norm(new_c - ctrs(i,:)) > tol)
	ctrs(i,:) = new_c;
	not_conv = 1;
      end
    else
%%-------------------- cluster is empty - reset centriod      
      ctrs(i,:) = rand_ftr(V);
      not_conv = 1;
   end
  end
  if (not_conv == 0), break, end
end

