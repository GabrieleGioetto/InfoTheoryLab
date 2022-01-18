
K = 0

for N = 1:11
    M = 12 - N;
    
    total_iterations = 200;
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

    subplot(3,4,N)
    cdfplot(Capacities_feedback)
    xlabel("Capacity")
    ylabel("Cumulative of capacity")
    title('Subplot 2: sin(2x)')
    a = strcat('M =', num2str(M));
    a = strcat(a, 'x');
    b = strcat('N =', num2str(N));
    title(strcat(a, b))
end