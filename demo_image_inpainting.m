% applying tnsor completion for image inpainting
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
    
X = double(imread('testimage.jpg')); 
X = X/255;
maxP = max(abs(X(:)));
dimX = size(X);

% sampling rate
p = 0.5

omega = find(rand(prod(dimX),1)<p);
M = zeros(dimX);
M(omega) = X(omega);

M2 = Frontal2Lateral(M); % each lateral slice is a channel of the image
omega2 = zeros(dimX);
Iones = ones(dimX);
omega2(omega) = Iones(omega);
omega2 = Frontal2Lateral(omega2);
omega2 = find(omega2==1);
n3 = size(M2,3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% L = dftmtx(n3); trransform.l = n3; transform.L = L;
% L = dct(eye(n3)); transform.l = 1; transform.L = L;
% L = RandOrthMat(n3); transform.l = 1; transform.L = L;

opts.DEBUG = 1;
Xhat2 = lrtc_tnn(M2,omega2,transform,opts);
Xhat2 = max(Xhat2,0);
Xhat2 = min(Xhat2,maxP);
Xhat2 = Lateral2Frontal(Xhat2); % each lateral slice is a channel of the image
psnr2 = PSNR(X,Xhat2,maxP)

