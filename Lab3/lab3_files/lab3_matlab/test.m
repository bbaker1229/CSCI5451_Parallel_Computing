close
%%clear
ni = 20;
nj = 20;
%%-------------------- east-west interactions
Tj = onedtor(0.15, 0.3, nj);
%%-------------------- north-south  * acts on i index
Ti = onedtor(0.20, 0.20, ni);
%%  
T = kron(speye(ni,ni),Tj)+kron(Ti,speye(nj,nj));
%%-------------------- adjust diagonal
d = 1 - sum(T,2);
T = T + diag(d);
n = size(T,1);

norm(T*ones(n,1) - ones(n,1))
%%-------------------- actual matrix is transpose.
T = T';
%%-------------------- get perron vector [use eigs for large
%%                     matrices !]
[X, D] = eig(T);
[d1, idx] = sort(abs(diag(D)),'descend');
u=X(:,idx(1));
u = u/sum(u);
%%-------------------- RESHAPE 
uu = reshape(u,ni,nj);
%%-------------------- Need to transpose!
uu = uu';

csvwrite('data.txt', uu);

surf(uu)
axis tight
 
