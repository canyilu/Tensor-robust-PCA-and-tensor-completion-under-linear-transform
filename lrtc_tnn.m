function [X,obj,err,iter] = lrtc_tnn(M,omega,transform,opts)

% Low-Rank Tensor Completion (LRTC) based on Tensor Nuclear Norm
% (TNN) under linear transform
%
% min_X ||X||_*, s.t. P_Omega(X) = P_Omega(M)
%
% ---------------------------------------------
% Input:
%       M       -    d1*d2*d3 tensor
%       L       -    d3*d3 matrix, L'*L=L*L'=l*I, l>0
%       omega   -    index of the observed entries
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
%   transform   - a structure which defines the linear transform
%       transform.L: the linear transform of two types:
%                  - type I: function handle, i.e., @fft, @dct
%                  - type II: invertible matrix of size n3*n3
%
%       transform.inverseL: the inverse linear transform of transform.L
%                         - type I: function handle, i.e., @ifft, @idct
%                         - type II: inverse matrix of transform.L
%
%       transform.l: a constant which indicates whether the following property holds for the linear transform or not:
%                    L'*L=L*L'=l*I, for some l>0.                           
%                  - transform.l > 0: indicates that the above property holds. Then we set transform.l = l.
%                  - transform.l < 0: indicates that the above property does not hold. Then we can set transform.l = c, for any constant c < 0.
%
%
%
% Output:
%       X       -    d1*d2*d3 tensor
%       err     -    residual
%       obj     -    objective function value
%       iter    -    number of iterations
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

tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if nargin < 3
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

dim = size(M);
X = zeros(dim);
X(omega) = M(omega);
E = zeros(dim);
Y = E;

for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    % update X
    [X,tnnX] = prox_tnn(-E+M+Y/mu,1/mu,transform); 
    % update E
    E = M-X+Y/mu;
    E(omega) = 0;
 
    dY = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnX;
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
end
obj = tnnX;
err = norm(dY(:));

 