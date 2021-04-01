function C = tprod_transform(A,B,L)

% Tensor-tensor product of two 3 way tensors under linear transform
% C = A *_L B
% A - n1*n2*n3 tensor
% B - n2*l*n3  tensor
% C - n1*l*n3  tensor
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

[n1,n2,n3] = size(A);
[m1,m2,m3] = size(B);

if n2 ~= m1 || n3 ~= m3
    error('Inner tensor dimensions must agree.');
end

A = lineartransform(A,L);
B = lineartransform(B,L);
C = zeros(n1,m2,n3);

if strcmp(L,'fft')
    % first frontal slice
    C(:,:,1) = A(:,:,1)*B(:,:,1);
    % i=2,...,halfn3
    halfn3 = round(n3/2);
    for i = 2 : halfn3
        C(:,:,i) = A(:,:,i)*B(:,:,i);
        C(:,:,n3+2-i) = conj(C(:,:,i));
    end    
    % if n3 is even
    if mod(n3,2) == 0
        i = halfn3+1;
        C(:,:,i) = A(:,:,i)*B(:,:,i);
    end
else    
    for i = 1 : n3
        C(:,:,i) = A(:,:,i)*B(:,:,i);
    end
end

C = inverselineartransform(C,L);
