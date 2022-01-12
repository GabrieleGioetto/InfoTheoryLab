
max_avg_no_feed = 0
max_avg_feed = 0
max_avg_no_feed_N = 0
max_avg_feed_N = 0

for N = 1:11
    M = 12 - N
    [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(N,M,20,200);

    N
    M
    avg_throughput_no_feedback = median(Capacities_no_feedback)
    avg_throughput_feedback = median(Capacities_feedback)

    if avg_throughput_no_feedback > max_avg_no_feed
        max_avg_no_feed = avg_throughput_no_feedback
        max_avg_no_feed_N = N
    end

    if avg_throughput_feedback > max_avg_feed
        max_avg_feed = avg_throughput_feedback
        max_avg_feed_N = N
    end
end

fprintf("---------------------------")

max_avg_no_feed
max_avg_feed
max_avg_no_feed_N
max_avg_feed_N




function [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(N,M,K, total_iterations)
    
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
end

