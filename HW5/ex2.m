%% part a ----------------------------
clear;
clc;

h = [1, 1.8, 0.81] ;
N =100 ; %lenght of u (input)
u = randn(N,1) ;
v = randn(N,1) ;
desired_signal = filter(h, 1, u) ;
Gamma = 1 ;
d = desired_signal + Gamma * v ;
mu = 0.5 ; % try mu = 0.2

taps = 4 ;
%taps = 5 ;

[Wopt_a, ~, ~] = LMS(u, d, taps, mu);


%% part b ----------------------------

Gamma = 0.1;
d = desired_signal + Gamma * v;
taps = 4;
mu_list = 0:0.01:0.5;
W_opt_target = [h , 0 ] .' ;   % for tap = 4
y = zeros(taps, numel(mu_list));   % store Wopt_b for each mu

figure;
hold on;
grid on;

for k = 1:numel(mu_list)
    mu = mu_list(k);

    [Wopt_b, ~, ~] = LMS(u, d, taps, mu);

    y(:, k) = Wopt_b - W_opt_target ;
end


plot(mu_list, y.', 'LineWidth', 2);
ylim([-5 5]);
xlabel('\mu');
ylabel('Error');
title('Error weight value');
legend("w 1","w 2","w 3","w 4",'Location','best');



