function y = MTF(A,B2,B1,C,max_iter,tol,~,T,n,w,NLfnc)

% State-space realization of the Lur'e-system:
% x_dot = A*x + B*u + H varphi(y)
% y     = C*x

G_yw = ss(A,B1,C,0);
G_yu = ss(A,B2,C,0);

% MTF parameters
y0            = zeros(n,1);   % initial guess for steady state response for e
yerror        = 1;            % set initial error > tol
k             = 0;            % set iteration index of while loop

% Memory allocation
y        = y0;                % output signal (error e in this case)
Y        = zeros(n,1);        % fft of error signal
Y_old    = zeros(n,1);        % used for determining convergence in MTF algorithm
f        = zeros(n,1);        % output of the nonlinearity (f = -phi(y))
U        = zeros(n,1);        % fft of output of the nonlinearity

% Corresponding frequency vector with frequency steps 1/T
f = (0 : 1 : n-1) * (1/T);

% Calculate input signal r(t) and Fourier coefficients R of r(t). Note that
% the spectrum repeats itself.
W  = fft(w).';

% Calculate Fourier coefficients of G_eu_bar and G_er_bar at frequencies
% [-f_N,f_N], where f_N is the Nyquist frequency (n/T)/2. Note that this
% spectrum repeats itself, see any signal-analysis book. (note that
% the negative frequencies are stored behind the positive ones,
% corresponding to the fft.m implementation.)
G_yu_vec = squeeze(freqresp(G_yu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));
G_yw_vec = squeeze(freqresp(G_yw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));

% Calculate initial response E0 of linear dynamics to input r with
% nonlinearity input u_bar set to zero: E = G_er_bar_vec*R + G_eu_bar_vec*U_bar =
% G_er_bar_vec*R + 0 = G_er_bar_vec*R. This is performed computationally efficient
% in frequency domain!
Y0 = G_yw_vec.*W;

while k < max_iter && yerror > tol
    if any(abs(Y)>1e10)
        display('Diverged')
        y = w'*0+nan;
        break;
    end
    u = -NLfnc(y);
    U = fft(u);                 % fft of output of nonlinearity
    Y = Y0 + G_yu_vec.*U;       % linear dynamics in frequency domain
    y = real(ifft(Y));          % use real value because there is a very small imaginary part present (1e-21)
    
    % error based on Fourier coefficients to check convergence
    yerror = norm(Y-Y_old)/norm(Y_old);
    nrm(k+1) = yerror;
    Y_old  = Y;
    
    % update the iteration index
    k     = k + 1;
end

if k == max_iter
    y = w'*0+nan;
    display('Max Iter reached')
end
