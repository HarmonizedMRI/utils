function th = getoephase(x)
% function ph = getoephase(x)
%
% Estimate odd/even EPI echo phase difference (linear and constant term)
% from one EPI echo train.
%
% Inputs
%  x    [nx etl nCoils]   EPI echo train of length 'etl' without phase-encoding.

verbose = true;

[nx etl nCoils] = size(x);

etl = size(x, 2);


% Estimate off-resonance phase from even echoes (center line)
% Only need one 'slice'
% Magnitude-squared coil weighting as in phase-contrast MRI
% Fit straight line and subtract from all echoes
th = zeros(1, etl/2);
for ic = 1:nCoils
    xe = x(end/2,2:2:etl,ic);  % even echoes
    tmp = unwrap(angle(xe));
    tmp = tmp - tmp(end/2);  % need common reference since we're combining coils
    th = th + abs(xe).^2 .* exp(1i*tmp);
end
if etl/2 < 10 
    maSpan = 1;
else
    maSpan = 5;   % moving-average span for smoothing
end
th = smooth(unwrap(angle(th)), maSpan);
X = [1:2:etl]';
B = [ones(length(th),1) X];
a = B\th(:);  % a(2) = off-resonance phase accrual per echo (radians)
XY = ones(nx,1) * [1:etl] - 1;  % [nx etl]
DPH = a(2)*repmat(XY, [1 1 nCoils]);  % linear fit to off-resonance phase accrual (radians)
xc = x.*exp(-1i*DPH);  % 'c' for 'corrected'

if verbose
    figure; subplot(121); im(angle(x(:,:,14)));
    subplot(122); im(angle(xc(:,:,14)));
    figure; ic = 14;
    plot(angle(x(end/2,:,ic)),'r'); hold on
    plot(angle(xc(end/2,:,ic)), 'g'); hold on
end

keyboard

return

pho = angle(x(end/2,1:2:end));
dph = angle(exp(1i*phe)./exp(1i*pho)); % / (2*pi*n*sys.raster);
offresPhase = mean(dph(2:end)); % exclude first echo
x(:,1:2:end) = x(:,1:2:end).*exp(1i*offresPhase);
%dCart(:,1:2:end,:) = dCart(:,1:2:end,:).*exp(1i*phOffset);

% get odd/even echo linear phase offset
phoe = angle(x(:,2:2:end)./x(:,1:2:end));
tmp = angle(x2./x1.*exp(1i*angle(x1./x3)/2));
th = th + abs(x1).^2.*exp(1i*tmp);
xsos = xsos + abs(x1).^2;

th = zeros(nx,1);
xsos = zeros(nx,1);  % sum-of-squares coil combined image (for mask)
	
for ic = 1:1:nCoils
    x1 = x(:,1:2:end);
    x2 = x(:,2:2:end);

	tmp = angle(x2./x1.*exp(1i*angle(x1./x3)/2));
	th = th + abs(x1).^2.*exp(1i*tmp);
	xsos = xsos + abs(x1).^2;
end

th = angle(th);
xsos = sqrt(xsos);

mask = xsos > 0.1*max(xsos(:));

% fit phase difference to 2d plane
H = [ones(sum(mask(:)),1) X(mask)];  % spatial basis matrix (2d linear)
ph(isl,:) = H\th(mask);  

thhat = embed(H*ph(isl,:)', mask);
figure; plot([th thhat th-thhat]); legend


