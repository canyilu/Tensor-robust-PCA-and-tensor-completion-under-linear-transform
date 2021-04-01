function A = lineartransform(A,L)

% linear transform along 3rd dim
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
    A = nmodetransform(A,'fwt',3);
elseif strcmp(L,'fft')
    A = fft(A,[],3);
elseif ismatrix(L)
    n3 = size(A,3);
    [l1,l2] = size(L);
    if l1 ~= l2 || l1 ~= n3
        error('Inner tensor dimensions must agree.');
    end
    A = tmprod(A,L,3);
end