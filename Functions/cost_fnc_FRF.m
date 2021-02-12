% Author: Fahim Shakib
% Date:   Feb. 12, 2021
% Email:  m.f.shakib@tue.nl

function [J,sysr] = cost_fnc_FRF(G,S1,S2,L1,L2,sysr,mag0,wout0,v2)

% Inputs:
% G         To-be-optimized parameters
% S1,S2     Signal generator system matrix
% L1,L2     Signal generator output matrix
% sysr      Reduced-order model
% mag0      Magnitude of Sigma
% wout0     Frequencies corresponding to mag0
% v2        Order of internal reduced-order LTI model

% Outputs:
% J         Value of cost function
% sysr      Reduced-order model with updated matrices

% Define G and F
G1      = G(1:length(G)-v2);
G2      = G(length(G)-v2+1:end);
F1      = S1-G1*L1;
F2      = S2-G2*L2;

% Update sysr
sysr.A = [F1 zeros(length(G)-v2,v2);zeros(v2,length(G)-v2) F2];
sysr.B = [G1 0*G1; 0*G2 G2];

% Compute FRF of sysr at wout0
[magr,~,~] = bode(sysr,wout0);

% Compute J
J = norm(squeeze(mag0-magr));