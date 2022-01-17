% ES 5

% Compute and plot separately the CDF of capacity for the case of MIMO Ricean
% (K = 20) channel at 10 dB SNR, with M = 3 transmit antennas and N = 3
% receive antennas.
% First: with no CSIT.
% Second: With CSIT (in the same plot).

% MIMO fading
N = 3;
M = 10;

K = 20

total_iterations = 200;
Capacities_no_feedback = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H_r = (1/sqrt(2)) * (randn(N, M) + 1i*randn(N, M));
    H_los = ones(N, M);

    H = sqrt(K / (K + 1)) * H_los + sqrt(1 / (K + 1))*H_r;


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
    H_r = (1/sqrt(2)) * (randn(N, M) + 1i*randn(N, M));
    H_los = ones(N, M);

    H = sqrt(K / (K + 1)) * H_los + sqrt(1 / (K + 1))*H_r;
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    u = zeros(1, min(M,N));
    
    trace_Q = 0;
    Q = zeros(min(N,M), min(N,M));

    [U,S,V] = svd(H);
    S = S .^ 2;

    while trace_Q < P
        prev_Q = Q;


        % get A
        a = (u - ones(1, min(M,N)) ./ diag(S)');
        % a = 0 if a < 0, a = a if a >= 0
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
legend("CSIT", "No CSIT")
xlabel("Capacity")
ylabel("Cumulative of capacity")
title('K = 20, MIMO: M = N = 3')


