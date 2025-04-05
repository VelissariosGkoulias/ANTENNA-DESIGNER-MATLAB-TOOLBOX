function [sigmadB,sigma] = circPlateRCS(r,f,theta)

lambda = physconst('lightspeed')/f;
k = 2*pi/lambda;
theta = theta*pi/180;

if isequal(theta,0)
   sigma = (4 * pi^3 * r^4)/lambda^2;
else
    arg = 2*k*r*sin(theta);
    sigma = (pi*(k^2)*(r^4)).*(2*besselj(1,arg)./arg).^2 .* (cos(theta)).^2;
end
sigmadB = 10*log10(sigma);