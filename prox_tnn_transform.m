function [X,tnn,trank] = prox_tnn_transform(Y,L,rho)

% The proximal operator of the tensor nuclear norm under linear transform
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% Y     -    n1*n2*n3 tensor
% L     -    n3*n3 matrix, L'*L=L*L'=l*I, l>0
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X
%
% version 1.0 - 01/02/2019
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
% References:
% Canyi Lu, Xi Peng, Yunchao Wei, Low-Rank Tensor Completion With a New Tensor 
% Nuclear Norm Induced by Invertible Linear Transforms. IEEE International 
% Conference on Computer Vision and Pattern Recognition (CVPR), 2019
%

[n1,n2,n3] = size(Y);
max12 = max(n1,n2);
X = zeros(n1,n2,n3);
Y = lineartransform(Y,L);

tnn = 0;
trank = 0;
for i = 1 : n3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    S = max(S-rho,0);
    tol = max12*eps(max(S));
    r = sum(S > tol);
    S = S(1:r);
    X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
    tnn = tnn+sum(S);
    trank = max(trank,r);
end
if strcmp(L,'fft')
    l = n3;
else
    l = norm(L(:,1))^2;
end
tnn = tnn/l;
X = inverselineartransform(X,L);
