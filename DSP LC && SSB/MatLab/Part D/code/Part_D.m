%% Sawtooth Waveform Generation
clc; clear all; close all;

%============================================
% PART D: SSB USING MATLAB
%============================================

%=================================================%
% [1] Generate and plot message signal m2(t):
%=================================================%

% Define the parameters
B = 1000; % Bandwidth in Hz
t = linspace(-0.005,0.005,2000); % Time vector

% Generate the message signal
m2 = sinc(B*t);

% Plot the message signal
figure;
plot(t,m2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Message signal m2(t)');
%%
%=================================================%
% [2] Generate the modulated signal, s2(t)
%=================================================%

% Define the carrier frequency
fc = 10000; % Hz

% Generate the Hilbert transform of the message signal
h_m2 = hilbert(m2);

% Generate the modulated signal
s2 = real(m2 .* cos(2*pi*fc*t) - imag(h_m2) .* sin(2*pi*fc*t));

% Plot the modulated signal
figure;
plot(t,s2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal s2(t)');

%%
%=================================================%
% [3] Plot the USB output
%=================================================%

% Generate the USB output
usb = s2 .* exp(1j*2*pi*fc*t);

% Plot the USB output
figure;
plot(t,real(usb));
xlabel('Time (s)');
ylabel('Amplitude');
title('USB Output');

%%
%=================================================%
% [4] Plot the LSB output
%=================================================%

% Generate the LSB output
lsb = s2 .* exp(-1j*2*pi*fc*t);

% Plot the LSB output
figure;
plot(t,real(lsb));
xlabel('Time (s)');
ylabel('Amplitude');
title('LSB Output');

%%
%=================================================%
% [5] Plot the spectrum of the modulated signal in both cases
%=================================================%
dt = t(2)-t(1);     %dt = 0.005/1999
fs = 1/dt;          %fs = 200 kHz

% Compute the Fourier transform of the modulated signal
S2 = fft(s2);

% Compute the frequency vector
f = linspace(-fs/2,fs/2,length(S2));

% Plot the magnitude spectrum of the modulated signal
figure;
plot(f,abs(fftshift(S2)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Modulated Signal s2(t)');

% Compute the Fourier transform of the USB output
USB = fft(usb);

% Plot the magnitude spectrum of the USB output
figure;
plot(f,abs(fftshift(USB)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of USB Output');

% Compute the Fourier transform of the LSB output
LSB = fft(lsb);

% Plot the magnitude spectrum of the LSB output
figure;
plot(f,abs(fftshift(LSB)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of LSB Output');



%%
%===================================================================%
% [6] Implement a suitable demodulator to extract m2(t) from s2(t)
%===================================================================%
% Generate the local oscillator at the carrier frequency
lo = cos(2*pi*fc*t);

% Multiply the modulated signal by the local oscillator
y = s2 .* lo;

% Design a low-pass filter with cutoff frequency of 1 kHz
f_cutoff = 1000;
[b,a] = butter(6,f_cutoff/(fs/2));

% Filter the demodulated signal
m2_demod = filtfilt(b,a,y);

% Plot the demodulated message signal
figure;
plot(t,m2_demod);
xlabel('Time (s)');
ylabel('Amplitude');
title('Demodulated Message Signal m2(t)');


%%
%===================================================================%
% [7] Investigate the output of the previous steps if the generator carrier
% wasn’t perfectly synchronized
%===================================================================%

% Case 1: Local carrier frequency at the receiver is f1 = fc + 0.1B.

% Define the local oscillator frequency with a frequency offset of 0.1*B
f1 = fc + 0.1*B;
lo1 = cos(2*pi*f1*t);

% Multiply the modulated signal by the local oscillator
y1 = s2 .* lo1;

% Filter the demodulated signal
m2_demod1 = filtfilt(b,a,y1);

% Plot the demodulated message signal with frequency offset
figure;
plot(t,m2_demod1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Demodulated Message Signal with Frequency Offset');


% Case 2: Local carrier frequency at the receiver is f1 = fc - 0.1B.
% Define the local oscillator frequency with a frequency offset of -0.1*B
f2 = fc - 0.1*B;
lo2 = cos(2*pi*f2*t);

% Multiply the modulated signal by the local oscillator
y2 = s2 .* lo2;

% Filter the demodulated signal
m2_demod2 = filtfilt(b,a,y2);

% Plot the demodulated message signal with frequency offset
figure;
plot(t,m2_demod2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Demodulated Message Signal with Frequency Offset');






% Create a Simulink model to visualize the sawtooth signal
t = transpose(t);
m2 = transpose(m2);
plot(t,m2);
% Create a Simulink model to visualize the sawtooth signal
%simin.time = t;
%simin.signals.values = m2;

%m2 = transpose(m2);
message = [t.*1000.+5 m2];





