% H: M x N ?

% H_los: all elements equal to 1

% SNR: N=1 for all antenna -> SNR = P (remember to convert from dB to units)
% to find linear from dB: x = 10^(dB/10)

% dB = 10*log_10(x)

% generate H using randN

% ES 1

% Compute and plot the CDF of capacity for the 
% case of SISO AWGN( K = infinite ?) and SISO Rayleigh fading 
% channel (K=0) at 10 dB SNR.
% SISO -> N = 1, M = 1

% SISO Rayleigh
N = 1;
M = 1;
% H must be complex
H = 1/2*(randn(1, 1) +1i*randn(1,1))

% Normalize
H_norm = H / sqrt( sum( H.^2 ) )


Identity_N = eye(N);
Identity_M = eye(M);

P = 10;
Q = (P/M)*Identity_M;

trace_Q = trace(Q)

C = log2(det(Identity_N + H_norm*Q*transpose(H_norm)))



