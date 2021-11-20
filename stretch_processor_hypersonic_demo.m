% Stretch Processor Demo
% 2, 6000m/s hypersonic targets, separated by 0.75m of range (5nS)
% example output from a stretch processor with slope correction applied
clear;
close all;

tau_p = 100e-6; % uncompressed pulse width
BW = 500e6; % LFM pulse bandwidth
c=3e8;
Ps = 1; % relative power of target #1
Ps2 = 0.5; % relative power of target #2
tau_m = 0;
tau_h = 104e-6; % % range extent of processor (predetermined by previous analysis)
Nfft = 2 ^ 19; 
sample_rate = 20e6; % ADC sample rate

scale_ = (tau_p/(0.25e-9)) / (1e9/sample_rate);

alpha = BW / tau_p; % LFM slope

tau_r1 = 0;     %initial relative target positions (seconds)
tau_r2 = 5e-9;

Rdot = -6000; %range rate, m/s
% 6km/s ~= 13400 miles/hour, 
% this velocity is on the upper end of current estimated 
% hypersonic missile capabilities

R1_0 = tau_r1 * c / 2;
R2_0 = tau_r2 * c / 2;

t = (-0.25 * Nfft / 2) * 1e-9: 0.25e-9 : 0.25 * (Nfft/2 - 1) * 1e-9;

R1 = R1_0 + t * c / 2;
R2 = R2_0 + t * c / 2;

alpha_r = alpha * (1 - 2 * Rdot / c) ^ 2;

tau_RRdot1 = (2*R1_0/c) / (1-2*Rdot/c);
tau_RRdot2 = (2*R2_0/c) / (1-2*Rdot/c);

phi_1 = pi * alpha_r * (t - tau_RRdot1) .^ 2;
phi_2 = pi * alpha_r * (t - tau_RRdot2) .^ 2;

% voltages out of receiver, after carrier is removed, before adc
v1 = sqrt(Ps) * exp(1j * phi_1).*rectpuls((t-tau_r1) / tau_p); 
v2 = sqrt(Ps2) * exp(1j * phi_2) .* rectpuls((t-tau_r1) / tau_p);

% Matched filter with slope correction applied
% To see output without slope correction, change alpha_r to alpha below
h = exp(1j*pi*alpha_r*t.^2).*rectpuls((t) / tau_h); 

r1 = v1;
r2 = v1 + v2;

vo1 = h.*conj(r1);
vo2 = h.*conj(r2);

sample_time = 1 / sample_rate; %sample time in seconds
sample_time_nS = sample_time * 1e9;

% t is a time array in nS
v_o_sampled1 = downsample(vo1, sample_time_nS);
v_o_sampled2 = downsample(vo2, sample_time_nS);

Nfft_downsampled = 2^17;
V_o_f1 = fftshift(abs(fft(v_o_sampled1, Nfft_downsampled)));
V_o_f2 = fftshift(abs(fft(v_o_sampled2, Nfft_downsampled)));

% FFT bin to frequency conversion
dt = sample_time_nS * 0.25e-9;
Dt = dt * Nfft_downsampled;
df = 1/Dt;
freqs = [-Nfft_downsampled/2:Nfft_downsampled/2-1]*df;
tau_r_calc = freqs/alpha;

figure1 = figure;
plot(tau_r_calc*1e9, V_o_f1/scale_);
xlim([-30 30]);
grid on;
xlabel('\tau (nS)');
titlestr = sprintf('1 Target, velocity=%.2f m/s', Rdot);
title(titlestr);
ylabel('Normalized Amplitude');
%movegui(figure2, [1800 700]);

figure2 = figure;
plot(tau_r_calc*1e9, V_o_f2/scale_);
titlestr = sprintf('2 Targets, velocity=%.2f m/s, separation=%.2f nS', Rdot, (tau_r2-tau_r1)*1e9);
title(titlestr);
xlim([-30 30]);
grid on;
xlabel('\tau (nS)');
ylabel('Normalized Amplitude');
%movegui(figure2, [1800 700]);




