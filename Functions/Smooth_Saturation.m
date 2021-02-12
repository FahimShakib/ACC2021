function [out] = Smooth_Saturation(u,par)

% Smooth saturation as a atan function

out = atan(par*u);