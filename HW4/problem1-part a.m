clc; clear; close all;
%% Example
R = [4 1; 
     1 3];
p = [1; 0.5];

Wo = R \ p;

w0 = [-2; 3];

[Q, Lambda] = eig(R);
lambda = diag(Lambda);
lambda_min = min(lambda);
lambda_max = max(lambda);

alpha_opt = 2/(lambda_min + lambda_max);

alphas = [0.05, 0.15, alpha_opt];

N = 40;

v0 = Q' * (w0 - Wo);

%% cost function
[w1, w2] = meshgrid(linspace(-3,4,300), linspace(-3,4,300));
J = zeros(size(w1));
for i = 1:numel(w1)
    wv = [w1(i); w2(i)];
    e = p - R*wv;
    J(i) = e'*e;
end

%% Plot contour
figure;
hold on;
grid on;
contour(w1, w2, J, 35, 'LineWidth',1);
plot(Wo(1), Wo(2), 'rp', 'MarkerSize',14, 'MarkerFaceColor','r');
text(Wo(1), Wo(2), 'Wo (optimum)');

%% Plot trajectories and first steps
colors = lines(length(alphas));
labels = cell(1,length(alphas));

for a = 1:length(alphas)
    alpha = alphas(a);
    path = zeros(2,N);  % 2 for R 2*2
    path(:,1) = w0;

    for n = 1:N
        M = (eye(2)-alpha*Lambda)^n; % M = (I−αΛ)
        w_n = Wo + Q*(M*v0); % w(n)= w0​ + Q (M)^n v(0)
        path(:,n+1) = w_n;
    end
    
    h(a) = plot(path(1,:), path(2,:), 'Color', colors(a,:), 'LineWidth',2);
    hold on
    plot(path(1,1:2), path(2,1:2), 'o', 'Color', colors(a,:), ...
        'MarkerFaceColor', colors(a,:), 'MarkerSize',7);
    labels{a} = sprintf('\\alpha = %.3f', alpha);
end

legend(h, labels);
xlabel('w_1');
ylabel('w_2');

%% in other way
errors = zeros(length(alphas), N);

for a = 1:length(alphas)
    alpha = alphas(a);
    for n = 1:N
        M = (eye(2)-alpha*Lambda)^n;
        w_n = Wo + Q*(M*v0);
        errors(a,n) = norm(w_n - Wo);
    end
end

figure;
hold on;
grid on;
for a = 1:length(alphas)
    plot(1:N, errors(a,:), 'LineWidth',2);
end
xlabel('Iteration'); ylabel('||w(n)-Wo||');
legend(arrayfun(@(x) sprintf('\\alpha = %.3f', x), alphas, 'UniformOutput', false));
