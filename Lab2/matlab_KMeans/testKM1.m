
clear; clc; close;

load("jumpDS");

[m,n,p] = size(X);

figure(1)
imagesc(X)

V = reshape(X,[m*n,p]);
V = double(V);

Nc = 3;
 
%%--------------------
[ctrs, idx] = myKmeans1(V,Nc,50);
disp('Kmeans out..')
%%-------------------- use indices as labels for plotting 
lab2 = fix(255*(idx-1)/Nc);
lab2 = reshape(lab2,m,n);

figure(2)
imagesc(lab2)
