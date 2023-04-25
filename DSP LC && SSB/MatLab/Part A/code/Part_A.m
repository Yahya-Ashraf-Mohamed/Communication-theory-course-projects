%% Sawtooth Waveform Generation
clc; clear all; close all;


%============================================
% PART A: DSB LC USING MATLAB
%============================================


% [1] Generate the message signal m1(t) and plot it.
% Define the parameters
sample_time = 1/100000;
steps = 1000000;

% Create the time vector
% t = 0:sample_time:(steps*sample_time);
t = sample_time*(0:steps);
t = t';
% Generate the sawtooth signal using the sawtooth function
sawtooth_signal = sawtooth((-2*pi*t+1.26*2.5));

% Plot the sawtooth signal
figure;
plot(t, sawtooth_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Sawtooth Signal');


% [2] Generate the modulated signal, s(t), using a carrier wave of 1 Volt amplitude and 10 KHz frequency, Ka = 0.5.
% s(t) = Ac [1 + kam1(t)] cos (2πfct)

% Define the parameters
Ac = 1; % Carrier amplitude
fc = 10000; % Carrier frequency in Hz
Ka = 0.5; % Modulation index

% Generate the modulated signal
st = Ac * (1 + Ka * sawtooth_signal) .* cos(2*pi*fc*t);

% Plot the modulated signal
figure;
plot(t, st);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal with Ka=0.5');

%=============================================

Ka = 1; % Modulation index

% Generate the modulated signal
st = Ac * (1 + Ka * sawtooth_signal) .* cos(2*pi*fc*t);

% Plot the modulated signal with Ka = 1
figure;
plot(t, st);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal with Ka=1');

%=============================================

Ka = 2; % Modulation index

% Generate the modulated signal
st = Ac * (1 + Ka * sawtooth_signal) .* cos(2*pi*fc*t);

% Plot the modulated signal with Ka = 2
figure;
plot(t, st);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal with Ka=2');

%============================================
% PART B: DSB LC USING SIMULINK
%============================================

% Create a Simulink model to visualize the sawtooth signal
simin.time = t;
simin.signals.values = sawtooth_signal;
