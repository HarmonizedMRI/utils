function dc = rampsamp2cart(dr, kx, nx, fov)
%
% Convert ramp-sampled trapezoidal raw data to Cartesian data using NUFFT (from MIRT).
%
% Inputs:
%  dr      [nr, ...]   Ramp-sampled raw data along 1st ("x") dimension 
%  kx      [nr 1]      k-space sampling locations (cycles/cm)
%  nx      [1]         Number of pixels (along x) in reconstructed image
%  fov     [1]         FOV (cm) along x
%
% Output:
%  dc      [nx, ...]   Data interpolated onto Cartesian grid along x

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

return


function sub_test

% 1d test object
nx = 128;
p = phantom(nx);
x = p(:,end/2);
fov = 20;  % cm

% nonuniform sampling (trapezoidal gradient)
dt = 4e-6;             % gradient raster time (s)
gamma = 4257.6;        % Hz/G
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (area of each readout trapezoid)
gmax = 1/(fov*gamma*dt);    % Gauss/cm
gslew = 10;      % G/cm/ms
gx = toppe.utils.trapwave2(2*area, gmax, gslew, dt*1e3);
gx = gx(2:(end-1));
kx = gamma*dt*cumsum(gx);
kx = kx - max(kx)/2;

% Ramp-sampled data
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov*kx(:)],mask,'nufft',nufft_args);
yr = A*x(:);

% create multi-dimensional input
yr = [yr 2*yr];

yc = rampsamp2cart(yr, nx, fov, kx);
%function dc = rampsamp2cart(dr, nx, fov, kxo, kxe)

%hold off; plot(x); hold on; plot(abs(xhat),'o'); legend('x true', 'xhat');

return
