function A = inverselineartransform(A,L)

% inverse linear transform along 3rd dim
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
    A = nmodetransform(A,'ifwt',3);
elseif strcmp(L,'fft')
    A = ifft(A,[],3);
elseif ismatrix(L)
    l = norm(L(:,1))^2;
    A = tmprod(A,L'/l,3);
end