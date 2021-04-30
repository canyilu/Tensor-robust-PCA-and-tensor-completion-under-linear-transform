% Low-Rank Tensor Completion (LRTC) based on TNN under linear transform
%
%
% version 1.0 - 01/02/2019
% version 1.1 - 29/04/2021
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
r = 5 % tubal rank

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% L = dftmtx(n3); transform.l = n3; transform.L = L;
% L = dct(eye(n3)); transform.l = 1; transform.L = L;
L = RandOrthMat(n3); transform.l = 1; transform.L = L;

X = tprod(randn(n1,r,n3)/n1,randn(r,n2,n3)/n2,transform);
trankX = tubalrank(X,transform)

dr = (n1+n2-r)*r*n3;
rho = 3
m = rho*dr;
p = m/(n1*n2*n3)

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

opts.DEBUG = 1;
Xhat = lrtc_tnn(M,omega,transform,opts);

trank = tubalrank(Xhat,transform);
RSE = norm(X(:)-Xhat(:))/norm(X(:));

fprintf('\nsampling rate: %f\n', p);
fprintf('tubal rank of the underlying tensor: %d\n', tubalrank(X,transform));
fprintf('tubal rank of the recovered tensor: %d\n', trank);
fprintf('relative recovery error: %.4e\n', RSE);


