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
lambda = 0.999 ;


[~, W_RLS, e_RLS] = RLS(n0, x, taps, lambda) ;

figure;

subplot(3,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_RLS(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('RLS output vs S');
legend('S (clean sound)','e RLS');
xlim([5 5.2]);

figure;

grid on;
plot(W_RLS.');
xlabel('Iteration');
ylabel('Weight value');
title('RLS weights');
legend([ arrayfun(@(k)sprintf('RLS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');
xlim([0 3000]);


sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_RLS(1:length(S)), fs) ;
pause(Duration + 0.5) ;

%% --------------- part v

var = 10 ;
new_n0 = sqrt(var) * n0 ;

n = filter(h, 1, new_n0) ;
x = S + n ;

[~, W_RLS, e_RLS] = RLS(new_n0, x, taps, lambda) ;

figure;

subplot(3,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_RLS(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('RLS output vs S');
legend('S (clean sound)','e RLS');
xlim([5 5.2]);

figure;

grid on;
plot(W_RLS.');
xlabel('Iteration');
ylabel('Weight value');
title('RLS weights');
legend([ arrayfun(@(k)sprintf('RLS w_%d',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');


sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_RLS(1:length(S)), fs) ;
pause(Duration + 0.5) ;

%% --------------- part z

b = [1 0.5];
a = [1 -0.9];
n = filter(b, a, n0) ;
x = S + n ;

taps = 25 ; % we need more taps !!!!
            % doesnt work well with 10 taps

[~, W_RLS, e_RLS] = RLS(n0, x, taps, lambda) ;

figure;

subplot(3,1,1);
plot(t, S);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('S (clean sound)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, x);
grid on;
xlabel('Time');
ylabel('Amplitude');
title('x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,3);
plot(t, S);
hold on;
grid on;
plot(t, e_RLS(1:length(S)));
xlabel('Time');
ylabel('Amplitude');
title('RLS output vs S');
legend('S (clean sound)','e RLS');
xlim([5 5.2]);

figure;

grid on;
plot(W_RLS.');
xlabel('Iteration');
ylabel('Weight value');
title('RLS weights');
legend([ arrayfun(@(k)sprintf('RLS w_{%d}',k),0:taps-1,'UniformOutput',false)],'Location','bestoutside');


sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e_RLS(1:length(S)), fs) ;
pause(Duration + 0.5) ;

