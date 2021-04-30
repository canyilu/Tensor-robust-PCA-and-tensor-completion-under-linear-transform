% Tensor robust principal component analysis based on tensor nuclear norm 
% under linear transform
%
% version 1.0 - 04/03/2019
% version 1.1 - 29/04/2021
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

n1 = 20;
n2 = n1;
n3 = 5;
r = 0.1*n1; % tubal rank
P = randn(n1,r,n3)/n1;
Q = randn(r,n2,n3)/n2;

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3; % not encourage to use, numerically unstable, may result to complex matrix
transform.L = dct(eye(n3)); transform.l = 1;
% % transform.L = RandOrthMat(n3); transform.l = 1;

L = tprod(P,Q,transform);
% L = real(L); % uncomment this line if transform dftmtx(n3) is used to avoid very small complex element

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
lambda = 1/sqrt(transform.l*max(n1,n2));

opts.tol = 1e-8;
opts.mu = 1e-4;
opts.rho = 1.1;
opts.DEBUG = 1;

[Lhat,Shat] = trpca_tnn(Xn,lambda,transform,opts);

Lr = norm(L(:)-Lhat(:))/norm(L(:))
Sr = norm(S(:)-Shat(:))/norm(S(:))
rankL = r
rankLhat = tubalrank(Lhat,transform)
sparsity = m
sparsityhat = length(find(Shat~=0))

