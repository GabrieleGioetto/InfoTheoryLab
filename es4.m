% ES 4

% Compute and plot separately the CDF of capacity for the case of MISO fading
% (K = 0) channel at 10 dB SNR, with M = 3 transmit antennas and N = 3 receive
% antennas.
% First: with no CSIT.
% Second: With CSIT (in the same plot).

% MIMO fading
N = 3;
M = 3;

total_iterations = 200;
Capacities_no_feedback = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H = (1/sqrt(2)) * (randn(N, M) + 1i*randn(N, M));
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    Q = (P/M)*Identity_M;
    
    trace_Q = trace(Q);
    
    Capacities_no_feedback(j) = log2(det(Identity_N + H*Q*ctranspose(H)));
end


Capacities_feedback = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H = (1/sqrt(2)) * (randn(N, M) + 1i*randn(N, M));
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    u = zeros(1, min(M,N));
    
    % TODO: get the values of the iteration before the last
    trace_Q = 0;
    Q = zeros(min(N,M), min(N,M));

    while trace_Q < P

        prev_Q = Q;

        [U,S,V] = svd(H);

        % x = 0 if x < 0, x = x if x >= 0
        S = S .^ 2;

        % get A
        a = (u - ones(1, min(M,N)) / diag(S)');
        a = max(a,0);
        A = diag(a);
        Q = V * A * ctranspose(V);
        
        trace_Q = trace(Q);

        u = u + 0.001;
    end
    
    % We get last value of Q that does not exceed power requirement
    Q = prev_Q;
    
    Capacities_feedback(j) = log2(det(Identity_N + H*Q*ctranspose(H)));
end

cdfplot(Capacities_feedback)
hold on;
cdfplot(Capacities_no_feedback)

