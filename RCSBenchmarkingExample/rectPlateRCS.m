function [sigmadB,sigma] = rectPlateRCS(a,b,f,phi,theta)
%Eq 2.62 Mahafza
lambda = physconst('lightspeed')/f;
k = 2*pi/lambda;
ph = phi*pi/180;
th = theta*pi/180;
arg1 = a*k*sin(th).*cos(ph);
arg2 = b*k*sin(th).*sin(ph);

sigma = ((4*pi*a^2*b^2)/lambda^2).*(sinc(arg1./pi) .* sinc(arg2./pi)).^2 .*((cos(th)).^2);

% sigma = ((4*pi*a^2*b^2)/lambda^2).* ((cos(th)).^2);
% sigma = ((4*pi*a^2*b^2)/lambda^2).*((sinc(arg1./pi)).^2).*(cos(th).^2);


% % eps = .000001;
% angle = k*a*sin(th);
% sigma = (4*pi*a^2*b^2/lambda^2).*(cos(th)).^2.*((sin(angle)./angle).^2);% + eps;
sigmadB = 10*log10(abs(sigma));