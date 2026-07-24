clc ;
clear ;
close all ;

rng(1);

[x, Fs] = audioread('sound.mp3') ;
if size(x,2) > 1
    x = mean(x,2) ;
end

Fs_target = 8000 ;
D = Fs / Fs_target ;
D = round(D) ;
S = x(1:D:end);
fs = Fs / D ;
S = S(320000:440000,1) ;
Duration = length(S) / fs ;
S = S -  mean(S) ;
S = S / std(S) ;

L = length(S) ;
n0 = randn(L,1) ;
t = (0:L-1).' / fs;

N = 5 ;
h = rand(N,1) ;
n = filter(h, 1, n0) ;
x = S + n ;

taps = 10 ;
mu = 0.1 ;

mu_min = 1e-12 ;
mu_max = 0.4 ;
p = 1.05 ;
method =  3 ;

[~, W_nlms, e_nlms] = NLMS(n0, x, taps, mu) ;
[~, W_vlms, e_vlms] = VLMS(n0, x, taps, mu , mu_min, mu_max, p, method) ;

figure;

subplot(4,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(4,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(4,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_nlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('NLMS output vs S');
legend('S (clean sound)','e NLMS');
xlim([5 5.2]);

subplot(4,1,4);
plot(t, S);
hold on;
grid on;
plot(t, e_vlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('VLMS output vs S');
legend('S (clean sound)','e VLMS');
xlim([5 5.2]);

figure;

subplot(2,1,1);
grid on;
plot(W_nlms.');
xlabel('Iteration');
ylabel('Weight value');
title('NLMS weights');
legend([ arrayfun(@(k)sprintf('NLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');

subplot(2,1,2);
grid on;
plot(W_vlms.');
xlabel('Iteration');
ylabel('Weight value');
title('VLMS weights');
legend([ arrayfun(@(k)sprintf('VLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');
xlim([1000 120000]);

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_nlms(1:length(S)), fs) ;
pause(Duration + 0.5) ;
sound(e_vlms(1:length(S)), fs) ;
pause(Duration + 1) ;

%% --------------- part v

var = 10 ;
new_n0 = sqrt(var) * n0 ;

n = filter(h, 1, new_n0) ;
x = S + n ;

[~, W_nlms, e_nlms] = NLMS(new_n0, x, taps, mu) ;
[~, W_vlms, e_vlms] = VLMS(new_n0, x, taps, mu, mu_min, mu_max, p,method) ;

figure;

subplot(4,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(4,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(4,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_nlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('NLMS output vs S');
legend('S (clean sound)','e NLMS');
xlim([5 5.2]);

subplot(4,1,4);
plot(t, S);
hold on;
grid on;
plot(t, e_vlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('VLMS output vs S');
legend('S (clean sound)','e VLMS');
xlim([5 5.2]);

figure;

subplot(2,1,1);
grid on;
plot(W_nlms.');
xlabel('Iteration');
ylabel('Weight value');
title('NLMS weights');
legend([ arrayfun(@(k)sprintf('NLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');

subplot(2,1,2);
grid on;
plot(W_vlms.');
xlabel('Iteration');
ylabel('Weight value');
title('VLMS weights');
legend([ arrayfun(@(k)sprintf('VLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');
xlim([5000 120000]);

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_nlms(1:length(S)), fs) ;
pause(Duration + 0.5) ;
sound(e_vlms(1:length(S)), fs) ;
pause(Duration + 1) ;

%% --------------- part z

b = [1 0.5];
a = [1 -0.9];
n = filter(b, a, n0) ;
x = S + n ;

taps = 20 ; % we need more taps !!!!
            % doesnt work with 10 taps

[~, W_nlms, e_nlms] = NLMS(n0, x, taps, mu) ;
[~, W_vlms, e_vlms] = VLMS(n0, x, taps, mu, mu_min, mu_max, p, method) ;

figure;

subplot(4,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(4,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(4,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_nlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('NLMS output vs S');
legend('S (clean sound)','e NLMS');
xlim([5 5.2]);

subplot(4,1,4);
plot(t, S);
hold on;
grid on;
plot(t, e_vlms(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('VLMS output vs S');
legend('S (clean sound)','e VLMS');
xlim([5 5.2]);

figure;

subplot(2,1,1);
grid on;
plot(W_nlms.');
xlabel('Iteration');
ylabel('Weight value');
title('NLMS weights');
legend([ arrayfun(@(k)sprintf('NLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');

subplot(2,1,2);
grid on;
plot(W_vlms.');
xlabel('Iteration');
ylabel('Weight value');
title('VLMS weights');
legend([ arrayfun(@(k)sprintf('VLMS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');
xlim([5000 120000]);

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_nlms(1:length(S)), fs) ;
pause(Duration + 0.5) ;
sound(e_vlms(1:length(S)), fs) ;
pause(Duration + 1) ;

