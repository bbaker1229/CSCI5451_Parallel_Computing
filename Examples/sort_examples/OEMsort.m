  function [ A ] = OEMsort(A)
%%---------- assumes n = power of 2
n = length(A);
%%-------------------- k = size of subsets
k = 1 ;
k2 = 2;
while (k2 <= n) 
    fprintf(1,' ** Outer loop * subsets size of %2d \n',k);
    num  = n/k2;
%%-------------------- show A at each outer loop
    A    
    for i=1:num
      j0 = (i-1)*k2+1;  j1 = j0+k;
fprintf(1,'interval [%2d : %2d ] [%2d : %2d] \n',j0,j0+k-1,j1,j1+k-1) 
      M = OEmerge(A(j0:j0+k-1),A(j1:j1+k-1));
      A(j0:j1+k-1) = M;
    end
    k = k2;
    k2 = 2*k;
end
