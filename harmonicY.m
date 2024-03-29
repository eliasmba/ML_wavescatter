function Ylm = harmonicY(n,m,theta,phi)
% spherical harmonics, Condon-Shortley sign convention
% for 0 <= theta <= pi only
% phi and theta are in radians
% if theta and phi are both vectors, result is a matrix with
% theta changing down the columns, phi changing along the rows
%
% Ylm = Ylm(l,m,theta,phi)
%
% caution, no checking on validity of input variables
theta = theta(:);
phi = phi(:)';
Pn = legendreP(n,cos(theta));
% sign change required for odd positive 
if m >= 0
  Pn = (-1)^m*Pn(abs(m)+1,:);
else
  Pn = Pn(abs(m)+1,:);
end
Ylm = (1/sqrt(2*pi))*Pn'*exp(1i*m*phi);
% output of legendre(n,x) is the vector P(n,m,x) for 0<=m<=n.
% P(n,-m,x) = (-)^m P(n,m,x) by convention. 
% Sign changes based on Matlab behavior that legendre(n,x) contains
% the Condon-Shortley factor (-)^m, but legendre(n,x,'norm') does not.