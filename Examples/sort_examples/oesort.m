   function [b] = oesort(a)
%% function [b] = oesort(a)  
%% odd-even sort. Sorts increasingly
  b = a;
  n = length(a);
%%-------------------- main loop (sequential) 
  for i=1:n
    start = 1+mod(i-1,2);
%%-------------------- parallel loop
    for j=start:2:n-1
%%-------------------- compare
      if (b(j)>b(j+1))
%%-------------------- exchange
	t = b(j);
	b(j) = b(j+1);
	b(j+1) = t;
      end
    end
    b
    pause 
    
  end
