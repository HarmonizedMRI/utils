function [w, pinvA] = slgcal(ycal, z, c, K, varargin)
% Return slice GRAPPA weights for one coil and one target slice
%
% Inputs:
%   ycal     [nx ny mb nc]     individually-acquired slice data; 
%                              non-zero values define the calibration region
%   z        [1]               target slice
%   c        [1]               target coil
%   K        [2]               kernel size (e.g., [5 5] or [7 7])
%
% Options:
%   lam      [1]        Tikhonov regularization constant for calculating w. Default: 0
%   pinvA               pinv(A), where A is a matrix constructed from ycal inside calibration region, see below.
%                       Default: construct from ycal and return to caller
%   slweights  [mb+1]   Determines slice weights when forming A. Default: [1 zeros(1,mb)] (traditional slice GRAPPA)
%                       slweights(1) = weight for collapsed slices
%                       slweights(2:mb+1) = weight for individual slices (leakage block)
%                       Examples:
%                         slweights = [1 zeros(1,mb)]: traditional slice GRAPPA (no leakage block)
%                         slweights = [0 ones(1,mb)]: split slice GRAPPA (leakage block)
%
% Output:
%   w        [K(1) K(2) nc]      slice-GRAPPA weights (set of nc 2D convolution kernels)
%
% The model here is:
%    A*w(:) = ycal([mask],z,c)
% where
%    'A' contains ycal data inside calibration region,
%    size(A) = [N, prod(K)*nc],
%    N = sum(mask(:)) = number of target calibration points,
%    mask = calibration region (non-zero entries in ycal)
%
% To use w to generate k-space data (for slice z, coil c) from
% acquired (collapsed) SMS data 'ysms', do:
%   y = zeros(size(ysms(:,:,1)));
%   for c = 1:nc
%       y = y + conv2(ysms(:,:,c), w(:,:,c), 'same');
%   end

if strcmp(ycal, 'test')
    tic; sub_test(); toc;
    return
end

% matrix size 
[nx ny mb nc] = size(ycal);

arg.lam = [];    
arg.pinvA = [];
arg.slweights = [0 ones(1,mb)];
arg = toppe.utils.vararg_pair(arg, varargin);

assert(xor(isempty(arg.lam), isempty(arg.pinvA)), 'Either lam or pinvA (or neither) must be provided');

if isempty(arg.lam)
    lam = 0;
else
    assert(length(arg.lam)==1, 'lam must be a scalar value');
    lam = arg.lam;
end

% check inputs
assert(length(K) == 2, 'K must be a length-2 vector');
assert(all(mod(K,2)==1), 'kernel size must be odd');

% calibration region
mask = abs(ycal(:,:,z,c)) > 0;       % calibration region (target points)
N = sum(mask(:));                    % number of calibration points
assert(N < nx*ny, 'calibration region must be smaller than k-space matrix size');
assert(N >= prod(K)*nc, 'number of calibration points must be >= prod(K)*nc (nc = number of coils)');

% construct A
if isempty(arg.pinvA)
    A = sub_getA(ycal, mask, K, arg.slweights, lam);
    pinvA = pinv(A); 
else
    pinvA = arg.pinvA;
end

% calculate w
tmp = ycal(:,:,z,c);
b = arg.slweights(1)*tmp(mask);
for iz = 1:mb
    tmp = ycal(:,:,iz,c);
    b = [b; (iz==z)*tmp(mask)];
end
w = pinvA*[b; zeros(prod(K)*nc,1)];
w = reshape(w, K(1), K(2), nc);

return


%% construct calibration matrix A
function A = sub_getA(ycal, mask, K, slweights, lam)
% Inputs
%   ycal     [nx ny mb nc]
%
% Output
%   A        [N*(mb+1)  prod(K)*nc]

mb = size(ycal, 3);
nc = size(ycal, 4);

N = sum(mask(:));              % number of calibration points

A = zeros(N*(mb+1), prod(K), nc);

A(1:N,:,:) = sub_getsubA(squeeze(sum(ycal,3)), mask, K, slweights(1));
for z = 1:mb
    A((1:N)+N*z,:,:) = sub_getsubA(squeeze(ycal(:,:,z,:)), mask, K, slweights(z+1));
end

A = reshape(A, N*(mb+1), prod(K)*nc);

A = [A; lam*eye(prod(K)*nc)];  % add regularization

return

function A = sub_getsubA(ycal1, mask, K, slweight)
% Input
%   ycal1     [nx ny nc]  
%
% Output
%   A         [N prod(K) nc]

nc = size(ycal1, 3);

N = sum(mask(:));              % number of calibration points

A = zeros(N, prod(K), nc);

% loop over points in calibration region
[I, J] = find(mask==1);  % target/center point indices
for n = 1:N
    % collect points inside the kernel centered at (I(n), J(n))
    col = 1;
    for jk = J(n)+floor(K(2)/2) : -1 : J(n)-floor(K(2)/2)
        for ik = I(n)+floor(K(1)/2) : -1 : I(n)-floor(K(1)/2)
            A(n, col, :) = slweight*ycal1(ik, jk, :);
            col = col + 1;
        end
    end
end

return


%% Old development code

function sub_test

nx = 128; ny = 96; nc = 32;
mb = 6;

K  = [5 5];

% synthesize truth data
w = randn(prod(K), nc) + 1i*randn(prod(K), nc);
ycal = randn(nx,ny,mb,nc) + 1i*randn(nx,ny,mb,nc);

% define calibration region 
ncalx = 2*ceil(sqrt(prod(K)*nc)/2);
ncaly = ncalx + 4;  % make it non-square for testing
mask = false(nx, ny);
mask(nx/2-ncalx/2:nx/2+ncalx/2-1, ny/2-ncaly/2:ny/2+ncaly/2-1) = true;

% Calibration matrix
z = 4; c = 20;
A = sub_getA(ycal, mask, K, [1 zeros(1,mb)]);        % [sum(mask) prod(K)*nc]

keyboard

return



function sub_test2

nx = 128; ny = 96; nc = 32;
mb = 6;

K  = [5 5];

% synthesize truth data
w = randn(prod(K), nc) + 1i*randn(prod(K), nc);
ycal = randn(nx,ny,mb,nc) + 1i*randn(nx,ny,mb,nc);

% define calibration region 
ncalx = 2*ceil(sqrt(prod(K)*nc)/2);
ncaly = ncalx + 4;  % make it non-square for testing
mask = false(nx, ny);
mask(nx/2-ncalx/2:nx/2+ncalx/2-1, ny/2-ncaly/2:ny/2+ncaly/2-1) = true;

% Calibration matrix
A = sub_getA(ycal, z, c, mask, K);        % [sum(mask) prod(K)*nc]

% synthesize calibration data
ycal = zeros(nx, ny, mb, c);
ycal(mask) = A*w(:);

% Check that we get w back with A\ycal(mask)
tol = 1e-9;
w_hat = A\ycal(mask);
w_hat = reshape(w_hat, size(w));  % for plotting
assert(abs(norm(w(:)-w_hat(:))) < tol, 'Test failed: A\ycal(mask) ~= w_true');
subplot(121); im(abs(w-w_hat));  title('w_true - A\ycal(mask)'); colorbar;

% Check that we get w back by calling slgcal() with lam=0
lam = 0;
w_hat2 = slgcal(ysms, ycal, K, 'lam', lam);
w_hat2 = reshape(w_hat2, prod(K), nc); 
assert(abs(norm(w(:)-w_hat2(:))) < tol, 'Test failed: slgcal() ~= w_true');
subplot(122); im(abs(w-w_hat2)); title('w_true - slgcal()'); colorbar;

% Check that we get ycal back (inside calibration region) when convolving dsms with w
w = reshape(w, K(1), K(2), nc);
ycal_hat = 0;
for c = 1:nc
    ycal_hat = ycal_hat + conv2(ysms(:,:,c), w(:,:,c), 'same');
end
ycal_hat(~mask) = 0;
assert(abs(norm(ycal_hat(:)-ycal(:))) < tol, 'Test failed: ycal_true ~= conv2(dsms, w)');
figure; im(abs(ycal_hat - ycal)); colorbar;
title('ycal - sum_coils(conv2(ysms(:,:,c), w(:,:,c)))');

% Visually inspect effect of regularization on reconstructed data
for lam = [0 5e1 20e1]
    w_hat3 = slgcal(ysms, ycal, K, 'lam', lam);
    y_hat3 = 0;
    for c = 1:nc
        y_hat3 = y_hat3 + conv2(ysms(:,:,c), w_hat3(:,:,c), 'same');
    end
    figure; im(y_hat3); colorbar; title(sprintf('lam = %d', lam));
end

return
