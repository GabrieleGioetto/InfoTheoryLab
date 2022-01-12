% EX 1

% Compute and plot the CDF of capacity for the 
% case of SISO AWGN( K = infinite ?) and SISO Rayleigh fading 
% channel (K=0) at 10 dB SNR.

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

% SISO Rayleigh
N = 1;
M = 1;

total_iterations = 200;

Capacities = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H = (1/sqrt(2))*(randn(N, M) + i*randn(N,M));
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    Q = (P/M)*Identity_M;
    
    trace_Q = trace(Q);
    
    Capacities(j) = log2(det(Identity_N + H*Q*ctranspose(H)));
end

cdfplot(Capacities)
xlabel("Capacity")
ylabel("Cumulative of capacity")
title('K = 0, SISO: M = 1, N = 1')


% EX 2

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
xlabel("Capacity")
ylabel("Cumulative of capacity")
title('K = 0, SIMO: M = 1, N = 3')


% EX 3

% Compute and plot separately the CDF of capacity for the case of MISO fading
% (K = 0) channel at 10 dB SNR, with M = 3 transmit antennas with no CSIT.

% MISO fading
N = 1;
M = 3;

total_iterations = 200;

Capacities = zeros(1, total_iterations);

for j = 1:total_iterations
    % H must be complex
    H = (1/sqrt(2)) * (randn(N, M) + 1i*randn(N, M));
        
    Identity_N = eye(N);
    Identity_M = eye(M);
    
    P = 10;
    Q = (P/M)*Identity_M;
    
    trace_Q = trace(Q);
    
    Capacities(j) = log2(det(Identity_N + H*Q*ctranspose(H)));
end

cdfplot(Capacities)
xlabel("Capacity")
ylabel("Cumulative of capacity")
title('K = 0, MISO: M = 3, N = 1')


% EX 4

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

        S = S .^ 2;

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
title('K = 0, MIMO: M = N = 3')

% ES 5

% Compute and plot separately the CDF of capacity for the case of MIMO Ricean
% (K = 20) channel at 10 dB SNR, with M = 3 transmit antennas and N = 3
% receive antennas.
% First: with no CSIT.
% Second: With CSIT (in the same plot).

% MIMO fading
N = 3;
M = 3;

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

    while trace_Q < P

        prev_Q = Q;

        [U,S,V] = svd(H);

        S = S .^ 2;

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



