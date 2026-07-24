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
lambda = 0.5 ;

taps = 4 ;

[Wopt_a_RLS, ~, ~] = RLS(u, d, taps, lambda);


%% part b ----------------------------

Gamma = 0.1;
d = desired_signal + Gamma * v;
taps = 4;
lambda_list = 0:0.01:1;
W_opt_target = [h , 0 ] .' ;
y_RLS = zeros(taps, numel(lambda_list));

figure;
hold on;
grid on;

for k = 1:numel(lambda_list)
    lambda = lambda_list(k);

    [Wopt_b_RLS, ~, ~] = RLS(u, d, taps, lambda);
    y_RLS(:, k) = Wopt_b_RLS - W_opt_target ;

end

plot(lambda_list, y_RLS.', 'LineWidth', 2);
ylim([-5 5]);
xlabel('\lambda');
ylabel('Error');
title('RLS - Error weight value');
legend("w 1","w 2","w 3","w 4",'Location','best');
grid on;

