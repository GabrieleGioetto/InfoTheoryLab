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

% SISO AWGN
N = 1;
M = 1;

total_iterations = 200;

Capacities = zeros(1, total_iterations);

for i = 1:total_iterations
    %H = ones(M,N) + norm(M,N);
    % H = randn(M,N);

    H = ones(N, M);
    
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    Q = (P/M)*Identity_M;
    
    trace_Q = trace(Q);
    
    Capacities(1, i) = log2(det(Identity_N + H*Q*ctranspose(H)));
end

cdfplot(Capacities)
xlabel("Capacity")
ylabel("Cumulative of capacity")
title('K = infinite, SISO: M = 1, N = 1')


% values = unique([-20 Capacities 20]);             %// obtain values where CDF should be computed
% F = cumsum(histc(Capacities, values))/numel(Capacities);   %// compute CDF
% stairs(values, F)                        %// do the plot
% axis([min(values) max(values) -.2 1.2])  %// adjust axes
% 
% fontSize = 15;
% title('cdf of Capacities (AWGN)', 'FontSize', fontSize);
% ylabel('cdf', 'FontSize', fontSize);
% xlabel('Capacity', 'FontSize', fontSize);
% grid on;
% 
