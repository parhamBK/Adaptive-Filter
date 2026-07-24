%% part a ----------------------------
clear;
clc;

rng(1) ;

h = [1, 1.8, 0.81] ;
N = 100 ;                 % length of u (input)
u_n = randn(N,1) ;
desired_signal = filter(h, 1, u_n) ;

% -------------
M = length(h) ;
K = N - M + 1 ;
U = zeros(K, M) ;
for i = 1:K
    U(i,:) = u_n(i+M-1:-1:i).' ;
end

d0 = desired_signal(M:end) ;

R_for_Jmin = (U.'*U) / K ;
p0 = (U.'*d0) / K ;
sigma_d2 = (d0.'*d0) / K ;
wopt = R_for_Jmin \ p0 ;
Jmin = sigma_d2 - p0.' * wopt ; % Since Jmin ≈ 0, the misadjustment is 0/0
                               % therefore, in this example the LMS
                               % algorithm converges
% -------------

R = (u_n) * (u_n') / (u_n' * u_n) ; % normalize matrix R with ||u(n)||^2
lambda = eig(R) ;
lambda_max = max(lambda) ;
mu_in_theory = linspace(0, 0.05, 2/lambda_max) ; % mu = μ

% -------but we can try in other way

W_opt = h ;
M = length(W_opt) ;
iter = 10 ;
E = zeros(M, N);

mu_test = 0:0.1:2 ;
means_mu_test_nlms = zeros(M,length(mu_test));

mu_min = 0.1;
mu_max = 0.4; % base on HW5

p = 0.1;
method = 1;

for i = 1:length(mu_test)
    mu = mu_test(i);

    E(:,:) = 0;
    for r = 1:iter
        u = randn(N,1);
        d = filter(h, 1, u);
        w = zeros(M,1);

        for n = M:N
            x = u(n:-1:n-M+1) ;
            y = w' * x ;
            e = d(n) - y ;
            w = w + (mu / (1e-6 + (x' * x))) * x * conj(e) ;
            E(:, n) = E(:, n) + (w - W_opt(:)) ;
        end
    end
    E = E / iter ;
    means_mu_test_nlms(:,i) = E(:,N) ;

end

x = mu_test(:) ;

figure;
Y = means_mu_test_nlms.' ;
bar(x, Y, 'grouped') ;
ylim([-5 5]) ;
xlabel('\mu') ;
ylabel('E(\epsilon)') ;
legend('Error w1','Error w2','Error w3','Location','best') ;
grid on ;
title('NLMS - mean convergence test');

%% part b ----------------------------
mu = 0.5;
taps = 4; % N = taps
[Wopt_for_tap4_nlms, ~ , ~] = NLMS(u_n, desired_signal, taps, mu);
[Wopt_for_tap4_vlms, ~ , ~] = VLMS(u_n, desired_signal, taps, mu, mu_min, mu_max, p);

%% part c ----------------------------

K = 5 ;
mu = 0.5; % watch mu = 0.2 too

W_opt_new = [W_opt , 0 ] .' ;

ek_nlms = zeros(taps, N - taps + 1) ;
ek_vlms = zeros(taps, N - taps + 1) ;

for i = 1:K

    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u);

    [~, W_steps_nlms, ~ ] = NLMS(u, d, taps, mu);
    L1 = size(W_steps_nlms,2);
    ek_nlms(:,1:L1) = ek_nlms(:,1:L1) + (W_steps_nlms - W_opt_new) ;

    [~, W_steps_vlms, ~ ] = VLMS(u, d, taps, mu, mu_min, mu_max, p);
    L2 = min(size(W_steps_vlms,2), size(ek_vlms,2));
    ek_vlms(:,1:L2) = ek_vlms(:,1:L2) + (W_steps_vlms(:,1:L2) - W_opt_new) ;

end

mean_e_nlms = (1/K) * ek_nlms ;
mean_e_vlms = (1/K) * ek_vlms ;

iters1 = 1:size(mean_e_nlms,2);
iters2 = 1:size(mean_e_vlms,2);

figure;
subplot(2,1,1);
plot(iters1, mean_e_nlms.') ;
grid on ;
xlabel('Iteration n');
ylabel('Mean Error weight value');
title('mean convergence - NLMS');
k = 1:taps;
labels = compose("w_%d", k);
legend(labels, 'Location', 'best');

subplot(2,1,2);
plot(iters2, mean_e_vlms.') ;
grid on ;
xlabel('Iteration n');
ylabel('Mean Error weight value');
title('mean convergence - VLMS');
k = 1:taps;
labels = compose("w_%d", k);
legend(labels, 'Location', 'best');

%% part d ----------------------------

ek2_nlms = zeros(1, N - taps + 1) ;
ek2_vlms = zeros(1, N - taps + 1) ;

for i = 1:K
    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u);

    [~, ~, error_nlms] = NLMS(u, d, taps, mu) ;
    temp = error_nlms.^2 ;
    ek2_nlms = ek2_nlms + temp(1:length(ek2_nlms)).' ;

    [~, ~, error_vlms] = VLMS(u, d, taps, mu, mu_min, mu_max, p) ;
    temp = error_vlms.^2 ;
    ek2_vlms = ek2_vlms + temp(1:length(ek2_vlms)).' ;
end

MSE_e_nlms = (1/K) * ek2_nlms ;
MSE_e_vlms = (1/K) * ek2_vlms ;

iters = 1:size(MSE_e_nlms,2) ;

figure ;
subplot(2,1,1);
plot(iters, MSE_e_nlms) ;
grid on ;
xlabel('Iteration n') ;
ylabel('MSE Error (output)') ;
title('MSE convergence - NLMS') ;

subplot(2,1,2);
plot(iters, MSE_e_vlms) ;
grid on ;
xlabel('Iteration n') ;
ylabel('MSE Error (output)') ;
title('MSE convergence - VLMS') ;

%% part e ----------------------------
mu = [0.1 , 0.05];
taps = 4; % N = taps
W_opt_new = [W_opt , 0 ] .' ; % for taps = 4
K = 5 ;

figure;

for m = 1:length(mu)

    mu_m = mu(m);
    ek_nlms = zeros(taps, N - taps + 1);
    ek_vlms = zeros(taps, N - taps + 1);

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u);

        [~, W_steps_nlms, ~] = NLMS(u, d, taps, mu_m);
        L1 = size(W_steps_nlms,2);
        ek_nlms(:,1:L1) = ek_nlms(:,1:L1) + (W_steps_nlms - W_opt_new);

        [~, W_steps_vlms, ~] = VLMS(u, d, taps, mu_m, mu_min, mu_max, p);
        L2 = min(size(W_steps_vlms,2), size(ek_vlms,2));
        ek_vlms(:,1:L2) = ek_vlms(:,1:L2) + (W_steps_vlms(:,1:L2) - W_opt_new);
    end

    mean_e_nlms = (1/K) * ek_nlms;
    mean_e_vlms = (1/K) * ek_vlms;

    iters1 = 1:size(mean_e_nlms,2);
    iters2 = 1:size(mean_e_vlms,2);

    subplot(length(mu), 2, (m-1)*2 + 1);
    plot(iters1, mean_e_nlms.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('NLMS - part e (\\mu = %.3f)', mu_m));
    k = 1:taps;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

    subplot(length(mu), 2, (m-1)*2 + 2);
    plot(iters2, mean_e_vlms.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('VLMS - part e (\\mu = %.3f)', mu_m));
    k = 1:taps;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

end

%% part f ----------------------------

mu = 0.3 ;
taps = [4 3 2] ; % N = taps
K = 5 ;

figure;

for t = 1:length(taps)

    taps_t = taps(t);

    temp = length(W_opt);
    if taps_t <= temp
        W_opt_new = W_opt(1:taps_t).';
    else
        W_opt_new = [W_opt , zeros(1, taps_t-temp)].' ;
    end

    ek_nlms = zeros(taps_t, N - taps_t + 1);
    ek_vlms = zeros(taps_t, N - taps_t + 1);

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u);

        [~, W_steps_nlms, ~] = NLMS(u, d, taps_t, mu);
        L1 = size(W_steps_nlms,2);
        ek_nlms(:,1:L1) = ek_nlms(:,1:L1) + (W_steps_nlms - W_opt_new);

        [~, W_steps_vlms, ~] = VLMS(u, d, taps_t, mu, mu_min, mu_max, p);
        L2 = min(size(W_steps_vlms,2), size(ek_vlms,2));
        ek_vlms(:,1:L2) = ek_vlms(:,1:L2) + (W_steps_vlms(:,1:L2) - W_opt_new);
    end

    mean_e_nlms = (1/K) * ek_nlms;
    mean_e_vlms = (1/K) * ek_vlms;

    iters1 = 1:size(mean_e_nlms,2);
    iters2 = 1:size(mean_e_vlms,2);

    subplot(length(taps), 2, (t-1)*2 + 1);
    plot(iters1, mean_e_nlms.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('NLMS - part f (tap = %.3f)', taps_t));
    k = 1:taps_t;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

    subplot(length(taps), 2, (t-1)*2 + 2);
    plot(iters2, mean_e_vlms.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('VLMS - part f (tap = %.3f)', taps_t));
    k = 1:taps_t;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

end
