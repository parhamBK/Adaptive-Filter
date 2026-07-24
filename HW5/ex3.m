%% part a ----------------------------
clear;
clc;
b = [1 0.5];
a = [1 -0.9];
N = 300;
mu = 0.3 ;
taps = 4 ;
 
u = Genarate_N_from_rand(N) ;  % Attached file 1

d = filter(b, a, u) ;

M = taps ;
iter = 30 ;
E = zeros(M, N - M + 1);
error = zeros(N, 1) ;
sum_Wopt_a = zeros(M, 1) ;

for r = 1:iter
    u = Genarate_N_from_rand(N) ;
    d = filter(b, a, u) ; % We will have a new inputs(u , d) at every iteration to compute the mean
                          % this is necessary for compute mean                   
    [ Wopt_a, W_steps , e] = LMS(u, d, taps, mu);
    error = error + (e);
    E = E + W_steps ;
    sum_Wopt_a = sum_Wopt_a + Wopt_a ;
end

E = E / iter ; 
mean_Wopt_a = sum_Wopt_a / iter ;
error = error / iter ;

Yw = E.';
nW = 0:size(Yw,1)-1;
fig1 = figure;
clf;
h1 = plot(nW, Yw, 'LineWidth', 1.2);
grid on;
ylim([-5 5]);
xlabel('n');
ylabel('W');
legend('w1','w2','w3','w4','Location','best');
title('Taps');


Ye = error ;
nE = 0:numel(Ye)-1;

fig2 = figure;
clf;
h2 = plot(nE, Ye, 'LineWidth', 1.2);
grid on;
xlabel('n');
ylabel('error');
title('Error signal - part a');

%% part b ----------------------------
taps_list = [2 4 8 16 32];

for ti = 1:numel(taps_list)

    taps = taps_list(ti);
    error = zeros(N, 1);

    for r = 1:iter

        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);

        [Wopt_a, W_steps, e] = LMS(u, d, taps, mu);
        error = error + (e);

    end

    error = error / iter;

   
    Ye = error;
    nE = 0:numel(Ye)-1;

    fig = figure('Name', sprintf('Error - %d taps', taps));
    clf;
    plot(nE, Ye, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal (taps = %d)', taps)) ;
end
%% part c ----------------------------
A_list = [0.1 0.3 1];
taps = 4 ;
M = taps;
iter = 30;

for ai = 1:numel(A_list)

    A = A_list(ai);
    error = zeros(N, 1);

    for r = 1:iter
        u = Genarate_N_from_rand(N);
        d = filter(b, a, u);
        v = randn(N,1) ;
        d = d + A * v ;
        [Wopt_a, W_steps, e] = LMS(u, d, taps, mu);

        error = error + (e);
    end

    error = error / iter;

    Ye = error;
    nE = 0:numel(Ye)-1;
    fig = figure('Name', sprintf('Error - A=%.3g', A));
    clf;
    plot(nE, Ye, 'LineWidth', 1.2);
    grid on;
    xlabel('n');
    ylabel('error');
    title(sprintf('Error signal, A = %.3g', A));

end








function u = Genarate_N_from_rand(N)
    u1 = rand(N,1);
    u2 = rand(N,1);
    x  = sqrt(-2*log(u1)) .* cos(2*pi*u2);  % ------Attached file 1
    u = x ./ max(abs(x));
end












