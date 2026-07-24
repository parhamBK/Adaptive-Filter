%% part a ----------------------------
clear;
clc;

h = [1, 1.8, 0.81] ;
N = 100 ;                 % length of u (input)

W_opt = h ;
M = length(W_opt) ;
iter = 10 ;
E = zeros(M, N);

lambda_test= 0.1:0.05:1.1 ;
means_lambda_test_RLS = zeros(M,length(lambda_test));

for i = 1:length(lambda_test)
    lambda = lambda_test(i);

    E(:,:) = 0;
    for r = 1:iter
        u = randn(N,1);
        d = filter(h, 1, u);
        w = zeros(M,1);
        [Wopt_RLS, ~ , ~] = RLS(u, d, M, lambda);
        E(:,i) = Wopt_RLS - W_opt.' ;
    end
    E = E / iter ;
    means_lambda_test_RLS(:,i) = E(:,i) ;
end

x = lambda_test(:) ;

figure;
Y = means_lambda_test_RLS.' ;
bar(x, Y, 'grouped') ;
ylim([-0.5 0.5]) ;
xlabel('\lambda') ;
ylabel('Error W') ;
legend('Error w1','Error w2','Error w3','Location','best') ;
grid on ;
title('RLS - lambda convergence test');

%% part b ----------------------------
h = [1, 1.8, 0.81] ;
N = 100 ;                 % length of u (input)
u_n = randn(N,1) ;
desired_signal = filter(h, 1, u_n) ;

lambda = 0.2;
taps = 4; % N = taps
[Wopt_for_tap4_RLS, ~ , ~] = RLS(u_n, desired_signal, taps, lambda);

%% part c ----------------------------

K = 5 ;
lambda = 0.5; % i change lambda into lambda and i will use same value for lambda
taps = 4 ;
W_opt_new = [W_opt , 0 ] .' ;

ek_RLS = zeros(taps, N - taps + 1) ;
for i = 1:K

    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u);

    [~, W_steps_RLS, E ] = RLS(u, d, taps, lambda);
    L1 = size(W_steps_RLS,2);
    ek_RLS(:,1:L1) = ek_RLS(:,1:L1) + (W_steps_RLS - W_opt_new) ;

end

mean_e_RLS= (1/K) * ek_RLS ;
mean_Error_RLS = (1/K) * E ;
iters1 = 1:size(mean_e_RLS,2);
iters2 = 1:size(mean_Error_RLS,1);
figure;
subplot(2,1,1);
plot(iters1, mean_e_RLS.') ;
grid on ;
xlabel('Iteration n');
ylabel('Mean Error weight value');
title('mean convergence with weight - RLS');
k = 1:taps;
labels = compose("w_%d", k);
legend(labels, 'Location', 'best');

subplot(2,1,2);
plot(iters2, mean_Error_RLS.') ;
grid on ;
xlabel('Iteration n');
ylabel('Mean Error');
title('mean convergence - RLS');
%% part d ----------------------------

ek2_RLS = zeros(1, N - taps + 1) ;

for i = 1:K
    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u);

    [~, ~, error_RLS] = RLS(u, d, taps, lambda) ;
    temp = error_RLS.^2 ;
    ek2_RLS = ek2_RLS + temp(1:length(ek2_RLS)).' ;
end

MSE_e_RLS = (1/K) * ek2_RLS ;

iters = 1:size(MSE_e_RLS,2) ;

figure ;
plot(iters, MSE_e_RLS) ;
grid on ;
xlabel('Iteration n') ;
ylabel('MSE Error (output)') ;
title('MSE convergence - RLS') ;
%% part e ----------------------------
lambda = [0.05 , 0.1 , 0.2 , 0.5, 0.6, 0.7 ];
taps = 4; % N = taps
W_opt_new = [W_opt , 0 ] .' ; % for taps = 4
K = 5 ;

figure;

for m = 1:length(lambda)

    lambda_m = lambda(m);
    ek_RLS = zeros(taps, N - taps + 1);

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u);

        [~, W_steps_RLS, ~ ] = RLS(u, d, taps, lambda_m);
        L1 = size(W_steps_RLS,2);
        ek_RLS(:,1:L1) = ek_RLS(:,1:L1) + (W_steps_RLS - W_opt_new);

    end

    mean_e_RLS = (1/K) * ek_RLS;

    iters1 = 1:size(mean_e_RLS,2);

    subplot(2*length(lambda), 1, m);
    plot(iters1, mean_e_RLS.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('RLS - part e (\\lambda = %.3f)', lambda_m));
    k = 1:taps;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

end

%% part f ----------------------------

lambda = 0.3 ;
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

    ek_RLS = zeros(taps_t, N - taps_t + 1);

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u);

        [~, W_steps_RLS, ~] = RLS(u, d, taps_t, lambda);
        L1 = size(W_steps_RLS,2);
        ek_RLS(:,1:L1) = ek_RLS(:,1:L1) + (W_steps_RLS - W_opt_new);
    end

    mean_e_RLS = (1/K) * ek_RLS;

    iters1 = 1:size(mean_e_RLS,2);

    subplot(length(taps), 1, t);
    plot(iters1, mean_e_RLS.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('RLS - part f (tap = %.3f)', taps_t));
    k = 1:taps_t;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

end
