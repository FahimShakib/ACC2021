%% Description
% This files contains the example in the ACC 2021 paper on model order
% reduction by moment matching for convergent Lur'e-type systems.

% Author: Fahim Shakib
% Date:   Feb. 12, 2021
% Email:  m.f.shakib@tue.nl

%% Initialization
clear all; clc
addpath('Functions')

% ACC paper example or new example?
load_data = 1;

%% Harmonic oscillator 1 dot{s1} = W1s1, u = L1w1 External
f01  = [0.01 0.1 0.38 1.32 4.1];                    % Base frequencies SG1
v1   = 2*length(f01);                               % State dimension SG1
tau01  = kron(ones(length(f01),1),[1; 0]);          % Initial condition SG1
% Construct S1
S1   = [];                                          
for k = 1:length(f01)
    S1 = blkdiag(S1,[0 1;-1 0]*2*pi*f01(k));
end
L1   = tau01';                                      % Output matrix SG1

% Check Assumption 2
if (rank(ctrb(S1,tau01)) + rank(obsv(S1,L1))) == 2*v1
    disp('SG1: Assumption 2 satisfied')
else
    disp('SG1: Assumption 2 not satisfied')
end

% State-space realization of SG1
SG1 = ss(S1,zeros(size(S1,1),1),L1,0);

%% Harmonic oscillator 2 dot{s2} = W2s2, r = L2w2 Internal
f02  = [0.01 0.1 0.38 1.27 4.1];                  % Base frequencies SG2
v2   = 2*length(f02);                             % State dimension SG2
tau02  = kron(ones(length(f02),1),[1; 0]);        % Initial condition SG2
% Construct S2
S2   = [];
for k = 1:length(f02)
    S2 = blkdiag(S2,[0 1;-1 0]*2*pi*f02(k));
end
L2   = tau02';                                    % Output matrix SG1

% Check Assumption 2
if (rank(ctrb(S2,tau02)) + rank(obsv(S2,L2))) == 2*v2
    disp('SG2: Assumption 2 satisfied')
else
    disp('SG1: Assumption 2 not satisfied')
end

% State-space realization of SG2
SG2 = ss(S2,zeros(size(S2,1),1),L2,0);

%% Input to MTF 
% MTF algorithm is used to simulate the Lur'e-type system below
% See [Pavlov A, Hunnekens BG, Wouw N, Nijmeijer H. Steady-state performance
% optimization for nonlinear control systems of Lur’e type. Automatica. 
% 2013 Jul 1;49(7):2087-97.] for the MTF algorithm

% MTF Settings
Tmtf      = 100;            % End time
ts        = 1e-3;           % Sampling time
tmtf      = 0:ts:Tmtf-ts;   % Time vector
max_iter  = 1000;           % # iterations in MTF
tol       = 1e-6;           % Convergence criterion of MTF

% Simulate SG1 to obtain input u
[u,~,~] = lsim(SG1,tmtf*0,tmtf,tau01);
u       = u';

%% True system is convergent with gamma = 1;
% First column B: External input
% Second column B: Internal input

if load_data
    load Data/Sys
else
    n       = 100;                      % State dimension
    sys0    = rss(n);                   % Generate random stable system
    sys0.D  = 0;
    sys0    = sys0/norm(sys0,inf)*0.9;  % Normalize its infinity gain to 0.9
    nrmB    = norm(sys0.B)/6;
    sys0.B  = [randn(n,1)*nrmB sys0.B]; % Take a random input matrix
    sys0.D  = [0 0];                    % Set D = 0 according

    delta   = 0.05;                     % Deadzone length
    NLfnc   = @(y)dz(y,delta);          % Define deadzone NL function
    % NLfnc   = @(y)Smooth_Saturation(y,delta);
end

% Assumption 2 holds automatically --> All poles in CLHP
if all(real(eig(sys0))<0)
    disp('Assumption 2 satisfied')
else
    disp('Assumption 2 not satisfied')
    return
end

%% Generate data NL System
% Perfrom simulation using MTG
y_NL    = MTF(sys0.A,sys0.B(:,2),sys0.B(:,1),sys0.C,max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);

%% ExoSystem (Linear Harmonic Oscillator + LTI Block)
% External
ex1  = ss([sys0.A sys0.B(:,1)*L1;zeros(v1,n) S1],zeros(n+v1,1),[sys0.C SG1.D*L1],0);
% Internal
ex2  = ss([sys0.A sys0.B(:,2)*L2;zeros(v2,n) S2],zeros(n+v2,1),[sys0.C SG2.D*L2],0);

%% Generate data LTI System 1: External
T1      = 100;              % End time
ts      = 1e-2;             % Sampling time
t1      = 0:ts:T1;          % Time vector
[y1,~,z1] = lsim(ex1,t1*0,t1,[zeros(n,1);tau01]);
x1      = z1(:,1:n);
w1      = z1(:,n+1:end);

%% Generate data LTI System 2: Internal
T2      = 100;              % End time
ts      = 1e-2;             % Sampling time
t2      = 0:ts:T2;           % Time vector
[y2,~,z2] = lsim(ex2,t2*0,t2,[zeros(n,1);tau02]);
x2      = z2(:,1:n);
w2      = z2(:,n+1:end);

%% Plot true response
figure
plot(t1,y1,t2,y2)
xlabel('Time [sec]')
ylabel('$y$')

%% Compute CPi (Simplified version of Theorem 6)
% External
tw1      = 20; % Time-window in seconds
CPi1     = y1(end+1-tw1/ts:end,:)'/w1(end+1-tw1/ts:end,:)';

% Internal
tw2      = 20; % Time-window in seconds
CPi2     = y2(end+1-tw2/ts:end,:)'/w2(end+1-tw2/ts:end,:)';

%% Initial reduced order models: Internal & External - LMI solution
gam  = 1;               % gamma^\star in paper

% Define sdp variables
Q1   = sdpvar(v1);
Q2   = sdpvar(v2);
X1   = sdpvar(v1,1);
X2   = sdpvar(v2,1);

% Define H
H1 = CPi1;
H2 = CPi2;

% Define constraints
C1 = [Q1*S1-X1*L1+S1'*Q1-L1'*X1'  gam*H1'*X2';
     gam*X2*H1                 S2'*Q2-L2'*X2'-gam*H2'*X2'+Q2*S2-X2*L2-gam*X2*H2] <= eps*eye(v1+v2);
 
C2 = [Q1*S1-X1*L1+S1'*Q1-L1'*X1'  gam*H1'*X2';
     gam*X2*H1                 S2'*Q2-L2'*X2'+gam*H2'*X2'+Q2*S2-X2*L2+gam*X2*H2] <= eps*eye(v1+v2);
C3 = Q2 >= eps*eye(v2);
C4 = Q1 >= eps*eye(v1);
constraints = C1+C2+C3+C4;

% Solve SDP the problem
sol             = solvesdp(constraints, []); 

% Retrieve solution
X2      = double(X2);
Q2      = double(Q2);
X1      = double(X1);
Q1      = double(Q1);
G1      = inv(Q1)*X1;
G2      = inv(Q2)*X2;
F1      = S1-G1*L1;
F2      = S2-G2*L2;

% Define reduced-order LTI models
sysr1   = ss(F1,G1,H1,0);
sysr2   = ss(F2,G2,H2,0);

%% Check for convergence
if norm(sysr2,inf) < 1
    display('Reduced-order model convergent')
else
    display('Reduced-order model not convergent')
end

%% Combined reduced order LTI model of order v1 + v2
sysrinit    = ss([F1 zeros(v1,v2);zeros(v2,v1) F2],[G1 0*G1; 0*G2 G2],...
                        [H1 H2],0);

%% Compute FRF of initial model
[mag0,~,wout0]          = bode(sys0);
[magrinit,~,woutrinit]  = bode(sysrinit);
Mag0                    = squeeze(20*log10(mag0));
Magrinit                = squeeze(20*log10(magrinit));

%% Optimize over G1 and G2
% Define cost function
J = @(G)cost_fnc_FRF(G,S1,S2,L1,L2,sysrinit,mag0,wout0,v2);

% Optimization settings
options                         = optimoptions('fminunc');
options.MaxFunctionEvaluations  = 10000;
% Perform optimization
Gopt                            = fminunc(J,[G1;G2],options);

% Construct reduced-order LTI model
[~,sysr] = J(Gopt);

%% Validate on training input - Simulate Lur'e-type model with MTF
psi    = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C,max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);

%% Validate on validation input
fv      = [0.02 0.25 2 4];                  % New Base frequencies
vv      = 2*length(fv);                     % State dimension
tau01v  = kron(ones(length(fv),1),[1; 0]);  % Initial condition

% Construct S
Sv   = [];
for k = 1:length(fv)
    Sv = blkdiag(Sv,[0 1;-1 0]*2*pi*fv(k));
end
Lv   = tau01v';                                      % Output matrix SG

% Define state-space representation of signal generator
SGv      = ss(Sv,zeros(size(Sv,1),1),Lv,0);
% Simulate signal generator
[uv,~,~] = lsim(SGv,tmtf*0,tmtf,tau01v);
uv       = uv';

% Perfrom simulation using MTF
y_NLv   = MTF(sys0.A,sys0.B(:,2),sys0.B(:,1),sys0.C,max_iter,tol,1,Tmtf,length(tmtf),uv,NLfnc);
psiv    = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C,max_iter,tol,1,Tmtf,length(tmtf),uv,NLfnc);

%% Combined plot training and validation test (Figure 2 in paper)
figure
subplot(211)
plot(tmtf,y_NL-mean(y_NL),tmtf,psi-mean(psi),'--')
hold all
plot(tmtf,y_NLv,tmtf,psiv-mean(psiv)+mean(y_NLv),'--')
l1 = legend('$\bar{y}_t$','$\bar{\psi}_t$','$\bar{y}_v$','$\bar{\psi}_v$',...
    'orientation','Horizontal');
xlabel('Time [sec]')
ylabel('Steady-State Output')
set(gca,'fontsize', 12)
set(findall(gcf,'type','line'),'linewidth',2)
l1.Position = [0.462 0.93 0.4438 0.0519];
ylim([-13 13])

subplot(212)
plot(tmtf,y_NL-mean(y_NL)-psi+mean(psi))
hold all
plot(tmtf,y_NLv-mean(y_NLv)-psiv+mean(psiv))
l2 = legend('$\bar{y}_t-\bar{\psi}_t$','$\bar{y}_v-\bar{\psi}_v$',...
    'orientation','Horizontal');
ylabel('Steady-State Error')
xlabel('Time [sec]')
set(gca,'fontsize', 12)
set(findall(gcf,'type','line'),'linewidth',2)
l2.Position = [0.57 0.455 0.3352 0.0519];
ylim([-0.3 0.3])

%% Plot bode (Figure 3 in paper)
% Compute FRF of reduced-order LTI model
[magr,~,woutr] = bode(sysr);
Magrf = squeeze(20*log10(magr));

% Display the FRF at base frequencies
disp(['Sys0(1) transfer: ' num2str(freqresp(sys0(1),2*pi*f01))])
disp(['Sysr(1) transfer: ' num2str(freqresp(sysr(1),2*pi*f01))])

disp(['Sys0(2) transfer: ' num2str(freqresp(sys0(2),2*pi*f02))])
disp(['Sysr(2) transfer: ' num2str(freqresp(sysr(2),2*pi*f02))])

% Compute FRF at base frequencies
magf1 = 20*log10(abs(squeeze(freqresp(sys0(1),2*pi*f01))));
magf2 = 20*log10(abs(squeeze(freqresp(sys0(2),2*pi*f02))));

figure;
subplot(121)
semilogx(wout0/2/pi,Mag0(1,:),'LineWidth',2);hold all
semilogx(woutrinit/2/pi,Magrinit(1,:),'LineWidth',2,'linestyle','--');
semilogx(woutr/2/pi,Magrf(1,:),'LineWidth',2,'linestyle',':');grid on;
semilogx(f01,magf1,'kx','MarkerSize', 12)
legend('$\Sigma$','$\Sigma_{r}^0$','$\Sigma_r$','$\Omega_0^1$',...
    'location','SW')
xlabel('Frequency [Hz]')
ylabel('Magnitude $y,\psi$ [dB]')
title('$u$')
set(gca,'fontsize', 12)
xlim([1e-2 1e1])
ylim([-50 12])
subplot(122)
semilogx(wout0/2/pi,Mag0(2,:),'LineWidth',2);hold all
semilogx(woutrinit/2/pi,Magrinit(2,:),'LineWidth',2,'linestyle','--');
semilogx(woutr/2/pi,Magrf(2,:),'LineWidth',2,'linestyle',':');grid on;
semilogx(f02,magf2,'kx','MarkerSize', 12)
legend('$\Sigma$','$\Sigma_{r}^0$','$\Sigma_r$','$\Omega_0^2$',...
    'location','SW')
title('$\varphi(y),\varphi(\psi)$')
xlabel('Frequency [Hz]')
xlim([1e-2 1e1])
ylim([-50 12])
set(gca,'fontsize', 12)