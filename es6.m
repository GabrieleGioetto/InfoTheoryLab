
max_avg_no_feed = 0;
max_avg_feed = 0;
max_avg_no_feed_N = 0;
max_avg_feed_N = 0;

max_10_LOS = 0;
max_10_LOS_N = 0;

for N = 1:11
    M = 12 - N;
    K = 50;

    [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(N,M,0,100);
    [Capacities_no_feedback_LOS, Capacities_feedback_LOS] = calculate_capacities(N,M,K,100);

    avg_throughput_no_feedback = median(real(Capacities_no_feedback));
    avg_throughput_feedback = median(real(Capacities_feedback));
    
    ten_throughput_LOS = prctile(real(Capacities_no_feedback_LOS),10);

    if avg_throughput_no_feedback > max_avg_no_feed
        max_avg_no_feed = avg_throughput_no_feedback;
        max_avg_no_feed_N = N;
    end

    if avg_throughput_feedback > max_avg_feed
        max_avg_feed = avg_throughput_feedback;
        max_avg_feed_N = N;
    end

    if ten_throughput_LOS > max_10_LOS
        max_10_LOS = ten_throughput_LOS;
        max_10_LOS_N = N;
    end

end

fprintf("---------------------------")

max_avg_no_feed
max_avg_feed

fprintf("---------------------------")

max_avg_no_feed_N
max_avg_feed_N

fprintf("---------------------------")

max_10_LOS
max_10_LOS_N


% 
% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(3,3,20,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)
% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(3,20,20,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)
% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(20,3,20,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)
% 
% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(20,20,20,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)
% 

% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(5,5,0,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)
% [Capacities_no_feedback, Capacities_feedback] = calculate_capacities(5,5,50,200);
% median(Capacities_no_feedback)
% median(Capacities_feedback)


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


%     cdfplot(Capacities_feedback)
%     hold on;
%     cdfplot(Capacities_no_feedback)
%     legend("CSIT", "No CSIT")
%     xlabel("Capacity")
%     ylabel("Cumulative of capacity")
%     title('K = 50/ K = 0, MIMO: M = N = 5')
% 
end

