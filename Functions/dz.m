function [dz_output] = dz(u,delta)

% deadzone nonlinearity

temp=abs(u)>delta;

dz_output = zeros(size(u));
dz_output = dz_output+temp.*(u-sign(u)*delta);
