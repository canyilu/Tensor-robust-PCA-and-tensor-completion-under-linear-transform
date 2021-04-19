% Tensor robust principal component analysis based on tensor nuclear norm 
% under linear transform
%
% version 1.0 - 04/03/2019
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
% References:
% Canyi Lu, Pan Zhou, Exact recovery of tensor robust principal component 
% analysis under linear transforms. arXiv preprint arXiv:1907.08288, 2019
%

addpath(genpath(cd))
clear
close all

n1 = 30;
n2 = n1;
n3 = n1;
r = 0.1*n1 % tubal rank
L1 = randn(n1,r,n3)/n1;
L2 = randn(r,n2,n3)/n2;

T = dct(eye(n3)); l=1;
% T = RandOrthMat(n3); l=1;

L = tprod_transform(L1,L2,T);

p = 0.1;
m = p*n1*n2*n3;
temp = rand(n1*n2*n3,1);
[B,I] = sort(temp);
I = I(1:m);
Omega = zeros(n1,n2,n3);
Omega(I) = 1;
E = sign(rand(n1,n2,n3)-0.5);

S = Omega.*E; % sparse part, or noises. S = P_Omega(E)

Xn = L+S;
lambda = 1/sqrt(l*max(n1,n2));

opts.tol = 1e-8;
opts.mu = 1e-4;
opts.rho = 1.1;
opts.DEBUG = 1;

tic
[Lhat,Shat] = trpca_tnn_transform(Xn,T,lambda,opts);

Lr = norm(L(:)-Lhat(:))/norm(L(:))
Sr = norm(S(:)-Shat(:))/norm(S(:))
rank = r
sparsity = m
sparsityhat = length(find(Shat~=0))

toc

