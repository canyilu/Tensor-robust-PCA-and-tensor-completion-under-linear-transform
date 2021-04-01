function trank = tubalrank_transform(X,L,tol)

% The tensor tubal rank of a 3 way tensor under linear transform
%
% X     -    n1*n2*n3 tensor
% L     -    n3*n3 matrix, L'*L=L*L'=l*I, l>0
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

if strcmp(L,'fwt')
    X = nmodetransform(X,'fwt',3);
elseif strcmp(L,'fft')
    X = fft(X,[],3);
elseif ismatrix(L)
    X = tmprod(X,L,3);    
end

[n1,n2,n3] = size(X);
S = zeros(n3,min(n1,n2));

for i = 1 : n3
    s = svd(X(:,:,i),'econ');
    S(i,:) = s';
end

if strcmp(L,'fwt')
    J = log2(n3);
    S = ifwt(S(:,:),'db8',J,n3);
elseif strcmp(L,'fft')
    S = ifft(S,[],3);
elseif ismatrix(L)
    l = norm(L(:,1))^2;
    S = L'*S/l;
end

if nargin==2
   tol = max(n1,n2)*n3*eps(max(abs(S(:))));
end
trank = sum(sum(abs(S)) > tol);

