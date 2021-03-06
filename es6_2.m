% ES 6.2
colors = ["r","b","g","k"];

clear h;

[h,Capacities_no_feedback, Capacities_feedback] = calculate_capacities(3,3,20,200, colors(1),1);
h_tot(1,1) = h(1,1);
h_tot(1,2) = h(1,2);
median(Capacities_no_feedback)
median(Capacities_feedback)

[h,Capacities_no_feedback, Capacities_feedback] = calculate_capacities(3,20,20,200, colors(2),2);
h_tot(2,1) = h(1,1);
h_tot(2,2) = h(1,2);
median(Capacities_no_feedback)
median(Capacities_feedback)

[h,Capacities_no_feedback, Capacities_feedback] = calculate_capacities(20,3,20,200, colors(3),3);
h_tot(3,1) = h(1,1);
h_tot(3,2) = h(1,2);
median(Capacities_no_feedback)
median(Capacities_feedback)

[h,Capacities_no_feedback, Capacities_feedback] = calculate_capacities(20,20,20,200, colors(4),4);
h_tot(4,1) = h(1,1);
h_tot(4,2) = h(1,2);
median(Capacities_no_feedback)
median(Capacities_feedback)

h_tot;

legend(h_tot(:,1), {'N = 3, M = 3' 'N = 3, M = 20' 'N = 20, M = 3' 'N = 20, M = 20'})



function [h, Capacities_no_feedback, Capacities_feedback] = calculate_capacities(N,M,K, total_iterations, color,i)

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


    h(1,1) = cdfplot(Capacities_feedback)
    hold on;
    h(1,2) = cdfplot(Capacities_no_feedback)
%     legend_text = sprintf('N = %d, M = %d',N, M);
%     legend(legend_text)
    set(h(1,1),'Color',color)
    set(h(1,2),'Color',color)


    xlabel("Capacity [bit/s]")
    ylabel("Cumulative of capacity")
% 
end



