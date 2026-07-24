%% part a ----------------------------
clear;
clc;

h = [1, 1.8, 0.81] ;
N =100 ;
u = randn(N,1) ;
v = randn(N,1) ;
desired_signal = filter(h, 1, u) ;
Gamma = 1 ;
d = desired_signal + Gamma * v ;
mu = 0.5 ;

taps = 4 ;

mu_min = 1e-4;
mu_max = 0.4;
p = 0.01;
method = 1;

[Wopt_a_nlms, ~, ~] = NLMS(u, d, taps, mu);
[Wopt_a_vlms, ~, ~] = VLMS(u, d, taps, mu, mu_min, mu_max, p);


%% part b ----------------------------

Gamma = 0.1;
d = desired_signal + Gamma * v;
taps = 4;
mu_list = 0:0.01:0.5;
W_opt_target = [h , 0 ] .' ;
y_nlms = zeros(taps, numel(mu_list));
y_vlms = zeros(taps, numel(mu_list));

figure;
hold on;
grid on;

for k = 1:numel(mu_list)
    mu = mu_list(k);

    [Wopt_b_nlms, ~, ~] = NLMS(u, d, taps, mu);
    y_nlms(:, k) = Wopt_b_nlms - W_opt_target ;

    [Wopt_b_vlms, ~, ~] = VLMS(u, d, taps, mu, mu_min, mu_max, p);
    y_vlms(:, k) = Wopt_b_vlms - W_opt_target ;
end

subplot(2,1,1);
plot(mu_list, y_nlms.', 'LineWidth', 2);
ylim([-5 5]);
xlabel('\mu');
ylabel('Error');
title('NLMS - Error weight value');
legend("w 1","w 2","w 3","w 4",'Location','best');
grid on;

subplot(2,1,2);
plot(mu_list, y_vlms.', 'LineWidth', 2);
ylim([-5 5]);
xlabel('\mu');
ylabel('Error');
title('VLMS - Error weight value');
legend("w 1","w 2","w 3","w 4",'Location','best');
grid on;
