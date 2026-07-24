%% part a ----------------------------
clear;
clc;
b = [1 0.5];
a = [1 -0.9];
N = 300;
lambda = 0.3 ;
taps = 4 ;

u = Genarate_N_from_rand(N) ;

d = filter(b, a, u) ;

M = taps ;
iter = 30 ;

E_RLS = zeros(M, N - M + 1);

error_RLS = zeros(N, 1) ;

sum_Wopt_a_RLS = zeros(M, 1) ;

for r = 1:iter
    u = Genarate_N_from_rand(N) ;
    d = filter(b, a, u) ;

    [Wopt_a_RLS, W_steps_RLS , e_RLS] = RLS(u, d, taps, lambda);

    error_RLS = error_RLS + e_RLS;

    L1 = size(W_steps_RLS,2);

    E_RLS(:,1:L1) = E_RLS(:,1:L1) + W_steps_RLS;
    sum_Wopt_a_RLS = sum_Wopt_a_RLS + Wopt_a_RLS ;
end

E_RLS = E_RLS / iter ;

mean_Wopt_a_RLS = sum_Wopt_a_RLS / iter ;

error_RLS = error_RLS / iter ;

Yw1 = E_RLS.';
nW1 = 0:size(Yw1,1)-1;
fig1 = figure;
clf;
plot(nW1, Yw1, 'LineWidth', 1.2);
grid on;
ylim([-5 5]);
xlabel('n');
ylabel('W');
legend('w1','w2','w3','w4','Location','best');
title('Taps - RLS');


Ye1 = error_RLS ;
nE1 = 0:numel(Ye1)-1;

fig3 = figure;
clf;
plot(nE1, Ye1, 'LineWidth', 1.2);
grid on;
xlabel('n');
ylabel('error');
title('Error signal - part a (RLS)');

%% part b ----------------------------
taps_list = [2 4 8 16 32];

for ti = 1:numel(taps_list)

    taps = taps_list(ti);

    error_RLS = zeros(N, 1);

    for r = 1:iter

        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);

        [~, ~, e_RLS] = RLS(u, d, taps, lambda);

        error_RLS = error_RLS + e_RLS;

    end

    error_RLS = error_RLS / iter;

    Ye1 = error_RLS;
    nE1 = 0:numel(Ye1)-1;

    fig = figure('Name', sprintf('Error RLS - %d taps', taps));
    clf;
    plot(nE1, Ye1, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal (RLS, taps = %d)', taps)) ;

end

%% part c ----------------------------
A_list = [0.1 0.3 1];
taps = 4 ;
M = taps;
iter = 30;

for ai = 1:numel(A_list)

    A = A_list(ai);

    error_RLS = zeros(N, 1);

    for r = 1:iter
        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);
        v = randn(N,1) ;
        d = d + A * v ;

        [~, ~, e_RLS] = RLS(u, d, taps, lambda);

        error_RLS = error_RLS + e_RLS;
    end

    error_RLS = error_RLS / iter;

    Ye1 = error_RLS;
    nE1 = 0:numel(Ye1)-1;
    fig = figure('Name', sprintf('Error RLS - A=%.3g', A));
    clf;
    plot(nE1, Ye1, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal, RLS, A = %.3g', A));


end


function u = Genarate_N_from_rand(N)
    u1 = rand(N,1);
    u2 = rand(N,1);
    x  = sqrt(-2*log(u1)) .* cos(2*pi*u2);
    u = x ./ max(abs(x));
end
