% SISO
N = 1;
M = 1;

H_los = ones(N,M);

SNRdB = 10;
SNR = 10^(SNRdB/10);

P = SNR*1;

%% TASK 1
% SISO AWGN - Case 1
I_m = eye(M);
Q = (P/M) .* I_m;
if trace(Q) == P
    disp('Check OK');
end

I_n = eye(N);
% select MR as number of Monte Carlo runs
MR = 200;
C_runs = zeros(MR,1);
for i = 1:MR
    H = H_los; %+ randn(1); adding gaussian noise
    % define capacity
    C_runs(i) = log2(det(I_n + H * Q* ctranspose(H)));
end
cdfplot(C_runs)
title('CDF of capacity - SISO AWGN')
figure(2)

% SISO Rayleigh fading channel

C_runs = zeros(MR,1);
for i = 1:MR
    H_r = (1/sqrt(2))*(randn(N,M) + 1i*randn(N,M));
    H = H_r + randn(1); 
    % define capacity
    C_runs(i) = log2(det(I_n + H * Q* ctranspose(H)));
end
cdfplot(C_runs)
title('CDF of capacity - SISO Rayleigh')