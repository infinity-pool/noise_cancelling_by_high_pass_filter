%% Read audio file
% Get audio data and sample rate
file_name = input("Type audio name : ", 's');
[y, Fs] = audioread("audio/"+file_name+".wav");
t = 0:length(y)-1;
save(file_name+".mat",'y','Fs','t');

%% Load mat file and initialize
clear;
close all;

mat_name = input("Type mat file name : ", 's');
load(mat_name+".mat", 'y', 'Fs', 't');

%% DFT of audio
% massive size vector to matrix
m = floor(sqrt(length(y)));
y = y(1:m*m);
y_mat = reshape(y, m, m);

% Show original audio
figure(1);
movegui("northwest");
plot(y);
title("Original signal : " + mat_name);

% DFT
F = [ ];
for i = 1:length(y_mat)
    F(:, end+1) = dft1(y_mat(:, i)');
end
F_vec = reshape(F', m*m, 1);

pause(0.1);
% Plot siganl in freq.-domain
figure(2);
movegui("north");
subplot(211);
stem(real(F_vec));
title("Original signal in Frequency Domain");

% High-pass Filter
filter_ratio = 0.137;
F_vec_filtered = F_vec;
len = length(F_vec);
mid = len / 2;
rect = ones(size(F_vec_filtered));
rect(floor(-filter_ratio * len / 2 + mid):floor(filter_ratio * len / 2 + mid)) = 0;
F_vec_filtered = rect.*F_vec_filtered;

% Plot filtered signal in freq.-domain
subplot(212);
stem(real(F_vec_filtered));
title("Filtered signal in Frequency Domain");

% reshape massive size vector to matrix
F_filtered = reshape(F_vec_filtered, m, m)';

% iDFT(original)
f_vec = [];
for i = 1:length(F)
    f_vec(:, end+1) = idft1(F(:, i)');
end
f_vec = reshape(f_vec, m*m, 1);

% iDFT(filtered)
f_filtered = [];
for i = 1:length(F_filtered)
    f_filtered(:, end+1) = idft1(F_filtered(:, i)');
end
f_vec_filtered = reshape(f_filtered, m*m, 1);

pause(0.1);
% Plot siganl in time-domain
figure(3);
movegui("northeast");
subplot(311);
plot(f_vec);
subplot(312);
plot(f_vec_filtered);
subplot(313);
audio_wo_noise = f_vec - f_vec_filtered;
plot(audio_wo_noise);

% Play original sound and filtered sound
sound(f_vec, Fs);
pause(3);
sound(y - f_vec_filtered, Fs);
audiowrite(mat_name + "_res.wav", audio_wo_noise, Fs);

%% Discrete Fourier Transform 1D
function [F] = dft1(time_domain_value)
M = length(time_domain_value);
x = -M/2:M/2-1;
n = -M/2:M/2-1;
power = -1i*2*pi/M*n'*x;
expn = exp(power);
F = (time_domain_value*expn)';
end

%% Inverse Discrete Fourier Transform 1D
function [F] = idft1(freq_domain_value)
M = length(freq_domain_value);
x = -M/2:M/2-1;
n = -M/2:M/2-1;
power = 1i*2*pi/M*x'*n;
expn = exp(power);
F = real(freq_domain_value*expn/M);
end