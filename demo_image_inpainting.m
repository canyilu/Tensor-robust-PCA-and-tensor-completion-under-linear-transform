% applying tnsor completion for image inpainting

clear
    
X = double(imread('testimage.jpg')); 
X = X/255;
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);

% sampling rate
p = 0.5

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

M2 = Frontal2Lateral(M); % each lateral slice is a channel of the image
omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
omega2(omega) = Iones(omega);
omega2 = Frontal2Lateral(omega2);
omega2 = find(omega2==1);

L = dct(eye(size(M2,3)));

Xhat2 = lrtc_tnn_transform(M2,L,omega2);
Xhat2 = max(Xhat2,0);
Xhat2 = min(Xhat2,maxP);
Xhat2 = Lateral2Frontal(Xhat2); % each lateral slice is a channel of the image
psnr2 = PSNR(X,Xhat2,maxP)

