clc ;
clear ;
close all ;

 rng(1);

[x, Fs] = audioread('sound.mp3') ;  %hmmmmm, i dont know how sample signal with fs = 8 KHz
                                    % and this function sample base on file
                                    % so i have this idea :                                  
if size(x,2) > 1    % is number of columns greater then 1 ? 
    x = mean(x,2) ; % if yes we make this vector
end

Fs_target = 8000 ;
D = Fs / Fs_target ;
D = round(D) ; % saftey ;)
S = x(1:D:end); % sample ~ 8KHz
fs = Fs / D ;

S = S(320000:440000,1) ; % not all song :)
Duration = length(S) / fs ;

S = S -  mean(S) ;
S = S / std(S) ;

L = length(S) ;
n0 = randn(L,1) ; % ----------ref input
t = (0:L-1).' / fs;

N = 5 ;
h = rand(N,1) ;
n = filter(h, 1, n0) ;
x = S + n ; % ----------primary input

taps = 10 ;
mu = 0.001 ;
[~, W, e] = problem2_LMS(n0, x, taps, mu) ;

% Compare signals
figure;

subplot(3,1,1);
plot(t, S);
hold on;
grid on;
plot(t, x);
xlabel('Time');
ylabel('Amplitude');
title('Before');
legend('S (clean sound)','x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, S);
hold on;
grid on;
plot(t, e);
xlabel('Time');
ylabel('Amplitude');
title('After');
legend('S (clean sound)','e after filtering');
xlim([5 5.2]);

subplot(3,1,3);
grid on;
plot(W.');
xlabel('Iteration');
ylabel('Weight value');
title('weights');
legend(arrayfun(@(k)sprintf('w_%d',k),0:taps-1,'UniformOutput',false),'Location','bestoutside');

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e, fs) ;
pause(Duration + 1) ;

%% --------------- part v

% we need just change n0 n and x

var = 10 ;
new_n0 = sqrt(var) * n0 ; % ----------ref input

n = filter(h, 1, new_n0) ;
x = S + n ; % ----------primary input

[~, W, e] = problem2_LMS(new_n0, x, taps, mu) ;

% Compare signals
figure;

subplot(3,1,1);
plot(t, S);
hold on;
grid on;
plot(t, x);
xlabel('Time');
ylabel('Amplitude');
title('Before (var = 10)');
legend('S (clean sound)','x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, S);
hold on;
grid on;
plot(t, e);
xlabel('Time');
ylabel('Amplitude');
title('After (var = 10)');
legend('S (clean sound)','e after filtering');
xlim([5 5.2]);

subplot(3,1,3);
grid on;
plot(W.');
xlabel('Iteration');
ylabel('Weight value');
title('weights (var = 10) ');
legend(arrayfun(@(k)sprintf('w_%d',k),0:taps-1,'UniformOutput',false),'Location','bestoutside');

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e, fs) ;
pause(Duration + 1) ;

%% --------------- part z

% we use n0 in first palace(std = 1)

b = [1 0.5];
a = [1 -0.9];
n = filter(b, a, n0) ;
x = S + n ; % ----------primary input

[~, W, e] = problem2_LMS(n0, x, taps, mu) ;

% Compare signals
figure;

subplot(3,1,1);
plot(t, S);
hold on;
grid on;
plot(t, x);
xlabel('Time');
ylabel('Amplitude');
title('Before (IIR)');
legend('S (clean sound)','x (sound with noise)');
xlim([5 5.2]);

subplot(3,1,2);
plot(t, S);
hold on;
grid on;
plot(t, e);
xlabel('Time');
ylabel('Amplitude');
title('After (IIR)');
legend('S (clean sound)','e after filtering');
xlim([5 5.2]);

subplot(3,1,3);
grid on;
plot(W.');
xlabel('Iteration');
ylabel('Weight value');
title('weights (IIR) ');
legend(arrayfun(@(k)sprintf('w_%d',k),0:taps-1,'UniformOutput',false),'Location','bestoutside');

sound(x, fs) ;
pause(Duration + 0.5) ;
sound(e, fs) ;
pause(Duration + 1) ;








