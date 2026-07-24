%% part a ----------------------------
clear;
clc;
b = [1 0.5];
a = [1 -0.9];
N = 300;
mu = 0.3 ;
taps = 4 ;

u = Genarate_N_from_rand(N) ;

d = filter(b, a, u) ;

M = taps ;
iter = 30 ;

E_nlms = zeros(M, N - M + 1);
E_vlms = zeros(M, N - M + 1);

error_nlms = zeros(N, 1) ;
error_vlms = zeros(N, 1) ;

sum_Wopt_a_nlms = zeros(M, 1) ;
sum_Wopt_a_vlms = zeros(M, 1) ;

mu_min = 1e-4;
mu_max = 0.4;
p = 0.01;
method = 1;

for r = 1:iter
    u = Genarate_N_from_rand(N) ;
    d = filter(b, a, u) ;

    [Wopt_a_nlms, W_steps_nlms , e_nlms] = NLMS(u, d, taps, mu);
    [Wopt_a_vlms, W_steps_vlms , e_vlms] = VLMS(u, d, taps, mu, mu_min, mu_max, p);

    error_nlms = error_nlms + e_nlms;
    error_vlms = error_vlms + e_vlms;

    L1 = size(W_steps_nlms,2);
    L2 = size(W_steps_vlms,2);

    E_nlms(:,1:L1) = E_nlms(:,1:L1) + W_steps_nlms;
    E_vlms(:,1:min(L2,size(E_vlms,2))) = E_vlms(:,1:min(L2,size(E_vlms,2))) + W_steps_vlms(:,1:min(L2,size(E_vlms,2)));

    sum_Wopt_a_nlms = sum_Wopt_a_nlms + Wopt_a_nlms ;
    sum_Wopt_a_vlms = sum_Wopt_a_vlms + Wopt_a_vlms ;
end

E_nlms = E_nlms / iter ;
E_vlms = E_vlms / iter ;

mean_Wopt_a_nlms = sum_Wopt_a_nlms / iter ;
mean_Wopt_a_vlms = sum_Wopt_a_vlms / iter ;

error_nlms = error_nlms / iter ;
error_vlms = error_vlms / iter ;

Yw1 = E_nlms.';
nW1 = 0:size(Yw1,1)-1;
fig1 = figure;
clf;
plot(nW1, Yw1, 'LineWidth', 1.2);
grid on;
ylim([-5 5]);
xlabel('n');
ylabel('W');
legend('w1','w2','w3','w4','Location','best');
title('Taps - NLMS');

Yw2 = E_vlms.';
nW2 = 0:size(Yw2,1)-1;
fig2 = figure;
clf;
plot(nW2, Yw2, 'LineWidth', 1.2);
grid on;
ylim([-5 5]);
xlabel('n');
ylabel('W');
legend('w1','w2','w3','w4','Location','best');
title('Taps - VLMS');

Ye1 = error_nlms ;
nE1 = 0:numel(Ye1)-1;

fig3 = figure;
clf;
plot(nE1, Ye1, 'LineWidth', 1.2);
grid on;
xlabel('n');
ylabel('error');
title('Error signal - part a (NLMS)');

Ye2 = error_vlms ;
nE2 = 0:numel(Ye2)-1;

fig4 = figure;
clf;
plot(nE2, Ye2, 'LineWidth', 1.2);
grid on;
xlabel('n');
ylabel('error');
title('Error signal - part a (VLMS)');

%% part b ----------------------------
taps_list = [2 4 8 16 32];

for ti = 1:numel(taps_list)

    taps = taps_list(ti);

    error_nlms = zeros(N, 1);
    error_vlms = zeros(N, 1);

    for r = 1:iter

        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);

        [~, ~, e_nlms] = NLMS(u, d, taps, mu);
        [~, ~, e_vlms] = VLMS(u, d, taps, mu, mu_min, mu_max, p);

        error_nlms = error_nlms + e_nlms;
        error_vlms = error_vlms + e_vlms;

    end

    error_nlms = error_nlms / iter;
    error_vlms = error_vlms / iter;

    Ye1 = error_nlms;
    nE1 = 0:numel(Ye1)-1;

    fig = figure('Name', sprintf('Error NLMS - %d taps', taps));
    clf;
    plot(nE1, Ye1, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal (NLMS, taps = %d)', taps)) ;

    Ye2 = error_vlms;
    nE2 = 0:numel(Ye2)-1;

    fig = figure('Name', sprintf('Error VLMS - %d taps', taps));
    clf;
    plot(nE2, Ye2, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal (VLMS, taps = %d)', taps)) ;
end

%% part c ----------------------------
A_list = [0.1 0.3 1];
taps = 4 ;
M = taps;
iter = 30;

for ai = 1:numel(A_list)

    A = A_list(ai);

    error_nlms = zeros(N, 1);
    error_vlms = zeros(N, 1);

    for r = 1:iter
        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);
        v = randn(N,1) ;
        d = d + A * v ;

        [~, ~, e_nlms] = NLMS(u, d, taps, mu);
        [~, ~, e_vlms] = VLMS(u, d, taps, mu, mu_min, mu_max, p);

        error_nlms = error_nlms + e_nlms;
        error_vlms = error_vlms + e_vlms;
    end

    error_nlms = error_nlms / iter;
    error_vlms = error_vlms / iter;

    Ye1 = error_nlms;
    nE1 = 0:numel(Ye1)-1;
    fig = figure('Name', sprintf('Error NLMS - A=%.3g', A));
    clf;
    plot(nE1, Ye1, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal, NLMS, A = %.3g', A));

    Ye2 = error_vlms;
    nE2 = 0:numel(Ye2)-1;
    fig = figure('Name', sprintf('Error VLMS - A=%.3g', A));
    clf;
    plot(nE2, Ye2, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal, VLMS, A = %.3g', A));

end


function u = Genarate_N_from_rand(N)
    u1 = rand(N,1);
    u2 = rand(N,1);
    x  = sqrt(-2*log(u1)) .* cos(2*pi*u2);
    u = x ./ max(abs(x));
end
