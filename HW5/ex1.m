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
p = (U.'*d0) / K ;
sigma_d2 = (d0.'*d0) / K ;
wopt = R_for_Jmin \ p ;
Jmin = sigma_d2 - p.' * wopt ; % Since Jmin ≈ 0, the misadjustment is 0/0
                               % therefore, in this example the LMS algorithm converges
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
means_mu_test = zeros(M,length(mu_test));

for i = 1:length(mu_test)
    mu = mu_test(i);
    for r = 1:iter
        u = randn(N,1);                      
        d = filter(h, 1, u);   % We will have a new inputs(u , d) at every iteration to compute(almost) the mean
                               % this is necessary for compute mean
        w = zeros(M,1) ;                       
    
        for n = M:N
            x = u(n:-1:n-M+1) ;                  
            y = w' * x ;
            e = d(n) - y ;
            w = w + mu * x * e ;  
            E(:, n) = E(:, n) + (w - W_opt(:)) ; % We sum over all iterations (for same n (time) )
        end
    end
    % and now we have sum of all result of random input 
    E = E / iter ; % tada! , obtain the mean of E = W(n) - W_opt (approximately)
    means_mu_test(:,i) = E(:,N) ;
end

x = mu_test(:) ;
Y = means_mu_test.' ;
% for mu = 0.5 is a little random some times converge and sometimes not
bar(x, Y, 'grouped') ;
ylim([-5 5]) ;
xlabel('\mu') ;
ylabel('E(\epsilon)') ;
legend('Error w1','Error w2','Error w3','Location','best') ;
grid on ;

% This is different from theory in the following way:
% in mean convergence, we assume mu is very small.
% With this assumption, we can assume independence and obtain simpler math.
% However, if μ becomes larger, the independence assumption is no longer valid

%% part b ----------------------------
mu = 0.5;
taps = 4; % N = taps
[Wopt_for_tap4, ~ , ~] = LMS(u_n, desired_signal, taps, mu);

%% part c ----------------------------

% hmmmm consider this e = ε and ek = εk so we have :

K = 5 ;
mu = 0.5; % watch mu = 0.2 too
ek = zeros(taps, N - taps + 1) ; % N is for length of u_n
W_opt_new = [W_opt , 0 ] .' ;

for i = 1:K

    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u); % every iter must have new random u

    [Wopt_c, W_steps, ~ ] = LMS(u, d, taps, mu);
    ek = ek + (W_steps - W_opt_new ) ;

end
mean_e = (1/K) * ek ;

iters = 1:size(mean_e,2);

figure;
plot(iters, mean_e.') ;
grid on ;
xlabel('Iteration n');
ylabel('Mean Error weight value');
title('mean convergence');

k = 1:taps;
labels = compose("w_%d", k);
legend(labels, 'Location', 'best');

% i tried with mu = 0.2, its better
%for mu = 0.5 is a random
% i dont know its random sometimes 30 iter sometimes more

%% part d ----------------------------

ek2 = zeros(1, N - taps + 1) ; % N is for length of u_n

for i = 1:K
    u = ( randn(N,1) ) ;
    d = filter(W_opt, 1, u); % every iter must have new random u

    [Wopt_c, W_steps, error] = LMS(u, d, taps, mu) ;
    temp = error.^2 ;
    ek2 = ek2 + temp ;
end

MSE_e = (1/K) * ek2 ;

iters = 1:size(MSE_e,1) ;

figure ;
plot(iters, MSE_e.') ;
grid on ;
xlabel('Iteration n') ;
ylabel('MSE Error (output)') ;
title('MSE convergence') ;


%% part e ----------------------------
mu = [0.1 , 0.05];
taps = 4; % N = taps
W_opt_new = [W_opt , 0 ] .' ; % for taps = 4
K = 5 ;

figure;

for m = 1:length(mu)

    mu_m = mu(m);

    [Wopt_for_tap4, ~, ~ ] = LMS(u_n, desired_signal, taps, mu_m);  % part b

    ek = zeros(taps, N - taps + 1); % N is for length of u_n

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u); % every iter must have new random u
        [Wopt_c, W_steps, ~] = LMS(u, d, taps, mu_m);
        ek = ek + (W_steps - W_opt_new);
    end

    mean_e = (1/K) * ek;
    iters = 1:size(mean_e,2);

    subplot(length(mu), 1, m);
    plot(iters, mean_e.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error **weight value**');
    title(sprintf('mean convergence - part e (\\mu = %.3f)', mu_m));

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

    %these line is for update W_opt with considering taps(number of rows)
    temp = length(W_opt);
    if taps_t <= temp
        W_opt_new = W_opt(1:taps_t).';
    else
        W_opt_new = [W_opt , zeros(1, taps_t-temp)].' ; 
    end




    [Wopt_for_tap4, ~, ~ ] = LMS(u_n, desired_signal, taps_t, mu);  % part b

    ek = zeros(taps_t, N - taps_t + 1); % N is for length of u_n

    for i = 1:K
        u = ( randn(N,1) );
        d = filter(W_opt, 1, u); % every iter must have new random u
        [Wopt_c, W_steps, ~] = LMS(u, d, taps_t, mu);
        ek = ek + (W_steps - W_opt_new);
    end

    mean_e = (1/K) * ek;
    iters = 1:size(mean_e,2);

    subplot(length(taps), 1, t);
    plot(iters, mean_e.');
    grid on;
    xlabel('Iteration n');
    ylabel('Mean Error weight value');
    title(sprintf('mean convergence - part f (tap = %.3f)', taps_t));

    k = 1:taps_t;
    labels = compose("w_%d", k);
    legend(labels, 'Location', 'best');

end







