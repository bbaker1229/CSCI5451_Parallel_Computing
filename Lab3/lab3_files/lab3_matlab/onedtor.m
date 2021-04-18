   function T = onedtor(w, e, n)
%% function T = onedtor(w, e, n)  
%% West, East, Center 
%% for vertical axis w means north (-1), e means south (+1)
%% off diagonal entries only! 
     for i=1:n
       pe = 0.25*(n-i+2)*e/n;
       pw = 0.25*(n+i-2)*w/n;
       ie = 1+mod(i,n);
       iw = 1+mod(i-2,n);
       T(i,ie) = pe;
       T(i,iw) = pw;
     end
     
