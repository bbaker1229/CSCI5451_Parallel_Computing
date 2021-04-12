 function [M] = OEmerge(A, B)
%%---------- Odd even merge 
n = length(A);
disp(' -->> OEmerge [A  B] = ')
disp([A, B])
if (n<=1) 
%%---------- sort if size <=2
  M = sort([A,B]) ;
  disp(M)
else
%%---------- E(A), O(B) 
  EA = A(1:2:n); 
  OB = B(2:2:n);
%%---------- O(A), E(B) 
  OA = A(2:2:n); 
  EB = B(1:2:n);
%%---------- recursive calls
  C = OEmerge(EA,OB);
  D = OEmerge(OA,EB);
  M = interleave(C, D) ;
  disp(' -->> Interleave C D ->  M=')
  disp(M)
  pause 
end
%%     M

