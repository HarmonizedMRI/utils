function dc = rampsamp2cart(dr, kx, nx, fov)
%
% Interpolate ramp-sampled data (e.g., EPI) to Cartesian grid.
%
% Inputs:
%  dr      [nr, ...]   Ramp-sampled raw data along 1st ("x") dimension 
%  kx      [nr 1]      k-space sampling locations (cycles/cm)
%  nx      [1]         Number of pixels (along x) in reconstructed image
%  fov     [1]         FOV (cm) along x
%
% Output:
%  dc      [nx, ...]   Data interpolated onto Cartesian grid along x
%
% To test:
% >> rampsamp2cart('test');

if strcmp(dr, "test")
	sub_test;
	return
end

nd = ndims(dr);

nr = size(dr,1);  % number of samples

drSize = size(dr);

kx = kx(:);

if nr ~= size(kx,1)
    error('Length of k-space and data (along 1st dimenstion) must match');
end

% Interpolate to Cartesian grid
osf = 2;   % oversampling factor
kxc = [ (-nx/2+1/(2*osf)):(1/osf):(nx/2-1/(2*osf)) ] * 1/fov;  % cycles/cm
d1 = interp1(kx, dr, kxc, 'spline', 'extrap');

% crop fov
x1 = fftshift(ifft(fftshift(d1,1), [], 1),1);
x1 = x1( (nx-nx/2+1):(nx+nx/2), :);
x1 = reshape(x1, [nx drSize(2:end)]);
dc = fftshift(fft(fftshift(x1,1), [], 1),1);

return

% The following uses NUFFT and avoids apodization, but is ~20x slower

% reshape to 2D for looping
dr = reshape(dr, nr, []);
nm = size(dr,2);

% Get Gmri object
[~, A, dcf] = reconecho([], nx, [], [], kx, fov);

% Interpolate to Cartesian grid
dc = zeros(nx, nm);  % Cartesian data
for ii = 1:nm
    x = reconecho(dr(:,ii), nx, A, dcf);
    dc(:,ii) = fftshift(fft(fftshift(x,1), [], 1),1);
end

% Reshape
dc = reshape(dc, [nx drSize(2:end)]);

return


function sub_test

% test object
nx = 128;
p = phantom(nx);
x = p(:,end/2);
x = [x 2*x 3*x];
fov = 20;  % cm

% nonuniform k-space (trapezoidal gradient)
dt = 4e-6;             % gradient raster time (s)
gamma = 4257.6;        % Hz/G
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (half-area of each readout trapezoid)
gmax = 1/(fov*gamma*dt);    % Gauss/cm
gslew = 10;      % G/cm/ms
gx = toppe.utils.trapwave2(2*area, gmax, gslew, dt*1e3);
gx = gx(2:(end-1));
kx = gamma*dt*cumsum(gx);
kx = kx - max(kx)/2;
%kx = linspace(kx(1), kx(end), nx);

% Synthesize ramp-sampled data
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov*kx(:)],mask,'nufft',nufft_args);
yr = A*x;

% Interpolated onto Cartesian grid 
yc = rampsamp2cart(yr, kx, nx, fov);

% Do IFT and compare with true object
xhat = fftshift(ifft(fftshift(yc,1), [], 1),1);
hold off; plot(x); hold on; plot(abs(xhat),'o'); legend('x true', 'xhat');

return


% Check normalization
figure
yfft = fftshift(fft(fftshift(x,1), [], 1),1);
plot(abs(yr)); hold on; plot(abs(yfft),'ro')

figure
xyr = fftshift(ifft(fftshift(yr,1), [], 1),1);
xyr2 = A'*yr/nx;
xyfft = fftshift(ifft(fftshift(yfft,1), [], 1),1);
plot(abs(xyr)); hold on; plot(abs(xyfft),'ro')
plot(abs(xyr2), 'gx');

yc = rampsamp2cart(yr, kx, nx, fov);
xyc = fftshift(ifft(fftshift(yc,1), [], 1),1);
plot(abs(xyc), 'bx');

