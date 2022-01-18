% H: N X M?

% H_los: all elements equal to 1

% SNR: N=1 for all antenna -> SNR = P (remember to convert from dB to units)
% to find linear from dB: x = 10^(dB/10)

% dB = 10*log_10(x)

% generate H using randN

% ES 2
% Compute and plot separately the CDF of capacity for the case of SIMO fading
% (K = 0) channel at 10 dB SNR, with N = 3 receive antennas.

% SIMO fading
N = 3;
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
title('K = 0, SIMO: M = 1, N = 3')

