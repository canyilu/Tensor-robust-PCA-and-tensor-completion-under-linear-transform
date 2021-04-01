% Low-Rank Tensor Completion (LRTC) under linear transform

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

clear

n = 50;
n1 = n;
n2 = n;
n3 = 10;
r = 25 % tubal rank
% L = dftmtx(n3);
% L = 'fft'
L = dct(eye(n3));
% L = RandOrthMat(n3);

X = tprod_transform(randn(n1,r,n3)/n1,randn(r,n2,n3)/n2,L);
trankX = tubalrank(X,L)

dr = (n1+n2-r)*r*n3;
rho = 3
m = rho*dr;
p = m/(n1*n2*n3)

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

opts.DEBUG = 1;
Xhat = lrtc_tnn_transform(M,L,omega,opts);

trank = tubalrank(Xhat,L);
RSE = norm(X(:)-Xhat(:))/norm(X(:));

fprintf('\nsampling rate: %f\n', p);
fprintf('tubal rank of the underlying tensor: %d\n', tubalrank(X,L));
fprintf('tubal rank of the recovered tensor: %d\n', trank);
fprintf('relative recovery error: %.4e\n', RSE);


