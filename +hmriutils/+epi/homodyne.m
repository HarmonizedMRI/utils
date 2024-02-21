function [yhat, xhat] = homodyne(y, ny, varargin)
% function [yhat, xhat] = homodyne(y, ny, varargin)
%
% Homodyne partial Fourier recon for 2D input
%
% Inputs
%  y      [nx etl]   2D k-space
%  ny     [1]        y matrix size (full)
%
% Outputs
%  yhat   [nx ny]    2D k-space
%
% Example usage: see sub_test() function in this script. To run it:
%  >> hmriutils.epi.homodyne('test');

% development version in github/jfnielsen/misc/partialFourier

if strcmp(y, 'test')
    sub_test;
    return
end

kmax = ny/2;
[nx, etl] = size(y);
k0 = etl-kmax;   % sampled range is -k0:kmax-1
k_cent = kmax + 1;

% parse input options
arg.mask = ones(nx, ny);
arg.plot = false;
arg.trw_lp = 10;
arg.trw_hp_lo = 6;
arg.trw_hp_hi = 6;
arg = toppe.utils.vararg_pair(arg, varargin);

% low pass filter, range is -k0:k0-1
lp((-kmax:-k0-1) + k_cent) = 0;
lp((-k0:k0-1) + k_cent) = 1;
lp((k0:kmax-1) + k_cent) = 0;
if arg.plot
    figure; hold on;
    plot(lp, 'bx-'); 
    title(sprintf('low pass filter (ny=%d, etl=%d)', ny, etl));
end
x = linspace(0, pi/2, arg.trw_lp);
lp(-k0 + (0:arg.trw_lp-1) + k_cent) = sin(x).^2;
lp(k0-1 - arg.trw_lp + (1:arg.trw_lp) + k_cent) = 1 - sin(x).^2;
if arg.plot
    plot(lp, 'ro-'); 
    legend('before', 'after');
end

% high pass filter, range -k0:kmax-1
hp((-kmax:-k0-1) + k_cent) = 0;
hp((-k0:k0-1) + k_cent) = 1;
hp((k0:kmax-1) + k_cent) = 2;
if arg.plot
    figure; hold on;
    plot(hp, 'bx-'); 
    title(sprintf('high pass filter (ny=%d, etl=%d)', ny, etl));
end
x = linspace(0, pi/2, arg.trw_hp_lo);
hp(-k0 + (0:arg.trw_hp_lo-1) + k_cent) = sin(x).^2;
%hp(k0-1 + (-trw/2:trw/2-1) + k_cent) = 1 + sin(x).^2;
x = linspace(0, pi/2, arg.trw_hp_hi);
hp(k0-1 + (0:arg.trw_hp_hi-1) - arg.trw_hp_hi/2 + 1 + k_cent) = 1 + sin(x).^2;
if arg.plot
    plot(hp, 'ro-'); 
    plot(hp, 'o'); title('high pass filter');
    legend('before', 'after');
end

% pad data to full size
yfull = zeros(nx, ny);
yfull(:, ny-etl+1:ny) = y;
xfull = fftshift(ifftn(fftshift(yfull))); % for comparison
if arg.plot
    figure; subplot(131); im(abs(xfull.*arg.mask)); title('zero-padded'); 
    subplot(132); im(angle(xfull).*arg.mask); title('phase'); 
    subplot(133); im(abs(yfull).^0.3); title('k-space');
end

% get phase from low-pass filtered image
ylp = yfull .* (ones(nx,1)*lp(:).');
xlp = fftshift(ifftn(fftshift(ylp)));
%figure; im(abs(ylp).^0.3);
phi = angle(xlp);

% high-pass filtered image
yhp = yfull .* (ones(nx,1)*hp(:).');
xhp = fftshift(ifftn(fftshift(yhp)));

% homodyne recon
xhat = real(xhp .* conj(xlp)./abs(xlp));

if arg.plot
    figure; subplot(121); im(abs(yhp).^0.3); title('hp');
    subplot(122); im(abs(xhp).^0.3); title('hp');
    figure; subplot(121); im(xhat); title('xhat');
    subplot(122); im(angle(xhat));
end

yhat = fftshift(fftn(fftshift(xhat)));

return

function sub_test

    % true object
    n = 90;
    mag = phantom(n)';
    mag(mag==1) = 0.3;
    mag = mag/max(mag(:));
    mask = mag > 0.01*max(mag(:));
    [X, Y] = meshgrid((-n/2:n/2-1)/n/2, (-n/2:n/2-1)/n/2);
    phs = X.*Y.^2;
    phs = 64*pi*phs/norm(phs(:)) + 1.0*mag.^1.0;
    xtrue = mag.*exp(1i*phs);

    % partial Fourier sampling
    yfull = fftshift(fftn(fftshift(xtrue)));
    etl = 72;
    k0 = n/2 - (n-etl);
    y = yfull(:,(n-etl+1):n);
    %im(abs(y).^0.3);

    % homodyne recon
    yhat = hmriutils.epi.homodyne(y, n, 'mask', mask, 'plot', true);
    xhat = fftshift(ifftn(fftshift(yhat)));

    % compare with true and zero-padded
    yzp = yfull;
    yzp(:,1:n-etl) = 0;
    xzp = fftshift(ifftn(fftshift(yzp)));
    figure;
    subplot(211);
    im(abs(cat(1, xtrue, xhat, xzp))); colorbar;
    title('Left to right: true; homodyne; zero-filled');
    subplot(212);
    im(cat(1, xtrue-xtrue, abs(abs(xtrue)-abs(xhat)), abs(abs(xtrue)-abs(xzp))));
    colorbar;
    title('difference');

return
