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
    u = 0;
    
    trace_Q = 0;
    Q = zeros(min(N,M), min(N,M));

    [U,S,V] = svd(H);
    S = S .^ 2;

    zero_matrix = zeros(M-N, M);

    while trace_Q < P

        prev_Q = Q;

        % get A
        A = u - (1 ./ S);

        % a = 0 if a < 0, a = a if a >= 0
        A = max(A, 0);
        A(A==Inf)=0;


        % add the zeros in case M > N
        if M > N
            A_updated = [A; zero_matrix];
        elseif M < N
            A_updated = A(1:M,1:M);
        else
            A_updated = A;
        end

        Q = V * A_updated * ctranspose(V);
        
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
xlabel("Capacity [bit/s]")
ylabel("Cumulative of capacity")
title('K = 0, MIMO: M = N = 3')

