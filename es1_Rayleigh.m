% H: N X M?

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

total_iterations = 200;

Capacities = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H = (1/sqrt(2))* (randn(N, M) + 1i*randn(N, M));
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    Q = (P/M)*Identity_M;
    
    trace_Q = trace(Q);
    
    Capacities(j) = log2(det(Identity_N + H*Q*ctranspose(H)));
end

cdfplot(Capacities)
xlabel("Capacity [bit/s]")
ylabel("Cumulative of capacity")
title('K = 0, SISO: M = 1, N = 1')





