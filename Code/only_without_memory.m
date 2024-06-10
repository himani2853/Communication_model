clc;
clear all;
[y, Fs] = audioread('project.wav');
y_mono = mean(y,2);
y_mono = y_mono(1:100);

% Fs = 48000;
% y_mono = [119/128,-9/128,10/128,115/128]; %............................x(t)

% sound(y_mono,Fs);

%% ADC conversion
M = 4;
levels = 2^8;
levels_found = linspace(min(y_mono), max(y_mono), levels); 
y_quantised = quantize_signal(y_mono,levels_found);
y_normalized = (y_quantised - min(y_quantised)) / (max(y_quantised) - min(y_quantised));
y_scaled = round((y_normalized)* 255);
y_check = de2bi(y_scaled);
y_binary = de2bi(y_scaled, 8);
y_plot = reshape(y_binary',1,[]);

figure; %............................x1(t)
stem(y_plot);
title('x1(t) - Analog to digital Conversion');
xlabel('sample index');
ylabel('sample amplitude');
grid on;
%% DAC Conversion

y_decimal = bi2de(y_binary);

%% Encoding

decimal_values = zeros(size(y_binary,1)*4, 1);

for i = 1:size(y_binary, 1)
    for j = 1:4
        bits = y_binary(i, (j-1)*2+1 : j*2);
        decimal_values((i-1)*4 + j) = bi2de(bits, 'left-msb');
    end
end

gray_code = [0, 1, 3, 2];

theta_0 = 0;
M = 4;
theta_m = zeros(4,1);

for i = 1:M
    m = i;
    theta_m(i) = theta_0 + ((2*pi)/M)*(m-1);
end

j = sqrt(-1);

encoded_angle = zeros(length(y_decimal),1);
encoded_signal = zeros(length(y_decimal),1);

for i= 1:length(decimal_values)
    gray_index = find(gray_code == decimal_values(i));
    encoded_angle(i) = theta_m(gray_index);
    encoded_signal(i) = exp(1j*encoded_angle(i));
end



%% line coding
%% a) rect pulse

T0 = (1/Fs)*2;
N = 19;
rect_pulse = ones(N,1);

inphase = zeros(length(encoded_signal),1);
quadphase = zeros(length(encoded_signal),1);

for i = 1: length(encoded_signal)
    inphase(i) = round(real(encoded_signal(i)));
    quadphase(i) = round(imag(encoded_signal(i)));
end


figure;
subplot(2,1,1); %........................x2(t)
stem(inphase);
title('x2(t) - Encoded inphase signal');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(quadphase);
title('x2(t) - Encoded quadphase signal');
xlabel('sample index');
ylabel('Amplitude');
grid on;

inphase_upsample = upsample(inphase, N);
quadphase_upsample = upsample(quadphase, N);


inphase_linecoded_signal = conv(inphase_upsample,rect_pulse);
quadphase_linecoded_signal = conv(quadphase_upsample,rect_pulse);

inphase_linecoded_signal = inphase_linecoded_signal(1:end-N);
quadphase_linecoded_signal = quadphase_linecoded_signal(1:end-N);

figure; %........................x3(t)
subplot(2,1,1);
plot(inphase_linecoded_signal);
title('x3(t) - Linecoded inphase signal (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(quadphase_linecoded_signal);
title('x3(t) - Linecoded inphase signal (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

figure;
plot(inphase_linecoded_signal,quadphase_linecoded_signal,'o');
title('Input constellation (Rect pulse)');
xlabel('inphase');
ylabel('quadphase');
grid on;

%% b) raised cosine
a = 1;
m = 9;
length_rc = 1;
[transmit_filter, ~] = raised_cosine(a,m,length_rc);

inphase_rc_linecoded_signal = conv(inphase_upsample,transmit_filter);
quadphase_rc_linecoded_signal = conv(quadphase_upsample,transmit_filter);

inphase_rc_linecoded_signal = inphase_rc_linecoded_signal(1:end-N);
quadphase_rc_linecoded_signal = quadphase_rc_linecoded_signal(1:end-N);


figure; %........................x3(t)
subplot(2,1,1);
plot(inphase_rc_linecoded_signal);
title('x3(t) - Linecoded inphase signal (Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(quadphase_rc_linecoded_signal);
title('x3(t) - Linecoded quadphase signal (Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

figure;
plot(inphase_rc_linecoded_signal,quadphase_rc_linecoded_signal,'o');
title('Input constellation (Raised cosine)');
xlabel('inphase');
ylabel('quadphase');
grid on;

%% modulation
%rect pulse

fc = 1e6;
t_rect = 0:1/(10*fc):(length(inphase_linecoded_signal)-1)/(10*fc);
A = 5;
Tb = 19;
factor = sqrt(2/Tb);
in_cos_rect = cos(2*pi*fc*t_rect);
quad_sin_rect = sin(2*pi*fc*t_rect);

modulated_rect_inphase = A*inphase_linecoded_signal.*in_cos_rect';
modulated_rect_quadphase = A*quadphase_linecoded_signal.*quad_sin_rect';

%raised cosine

t_rc = 0:1/Fs:(length(inphase_rc_linecoded_signal)-1)/Fs;
in_cos_rc = cos(2*pi*fc*t_rc);
quad_sin_rc = sin(2*pi*fc*t_rc);

modulated_rc_inphase = A*inphase_rc_linecoded_signal.*in_cos_rc';
modulated_rc_quadphase = A*quadphase_rc_linecoded_signal.*quad_sin_rc';


figure; %........................x4(t)
subplot(2,1,1);
plot(modulated_rect_inphase);
title('x4(t) - modulated inphase signal (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(modulated_rect_quadphase);
title('x4(t) - modulated quadphase signal (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

figure;
subplot(2,1,1);
plot(modulated_rc_inphase);
title('x4(t) - modulated quadphase signal (Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(modulated_rc_quadphase);
title('x4(t) - modulated quadphase signal (Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

final_modulated_rect = modulated_rect_inphase + (1i)*modulated_rect_quadphase;
final_modulated_rc = modulated_rc_inphase + (1i)*modulated_rc_quadphase;



%% Memoryless AWGN Channel

sigma_square = 1;
j = sqrt(-1);
noise_rect = sqrt(sigma_square)*(randn((size(final_modulated_rect)))+ j*(randn(size(final_modulated_rect))));
received_signal_rect = final_modulated_rect + noise_rect;

noise_rc = sqrt(sigma_square)*(randn((size(final_modulated_rc)))+ j*(randn(size(final_modulated_rc))));
received_signal_rc = final_modulated_rc + noise_rc;

figure;%.....................y4(t)
subplot(2,1,1);
plot(real(received_signal_rect));
title('y4(t) - modulated inphase signal with noise ml(Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(imag(received_signal_rect));
title('y4(t) - modulated quadphase signal with noise ml(Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;


figure;
subplot(2,1,1);
plot(real(received_signal_rc));
title('y4(t) - modulated inphase signal with noise ml(Raised cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(imag(received_signal_rc));
title('y4(t) - modulated quadphase signal with noise ml(Raised cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

%% Demodulation

%Memoryless Channel
%Rect pulse

received_signal_rect_cos = 2*received_signal_rect.*cos(2*pi*fc*t_rect)';
received_signal_rect_sine = 2*received_signal_rect.*sin(2*pi*fc*t_rect)';

received_inphase = lowpass(real(received_signal_rect_cos),5500,Fs);
received_quadphase = lowpass(imag(received_signal_rect_sine),5500,Fs);

figure;%........................y3(t)
subplot(2,1,1);
plot(received_inphase);
title('y3(t) - demodulated quadphase signal with noise ml(Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(received_quadphase);
title('y3(t) - demodulated quadphase signal with noise ml(Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;


figure;
plot(received_inphase,received_quadphase,'o');
title('Output constellation with memoryless noise (Rect pulse)');
xlabel('inphase');
ylabel('quadphase');
grid on;


%Raised Cosine pulse

received_signal_rc_cos = 2*received_signal_rc.*cos(2*pi*fc*t_rc)';
received_signal_rc_sine = 2*received_signal_rc.*sin(2*pi*fc*t_rc)';

received_rc_inphase = lowpass(real(received_signal_rc_cos),5500,Fs);
received_rc_quadphase = lowpass(imag(received_signal_rc_sine),5500,Fs);

figure;
subplot(2,1,1);
plot(received_rc_inphase);
title('y3(t) - demodulated quadphase signal with noise ml(Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(received_rc_quadphase);
title('y3(t) - demodulated quadphase signal with noise ml(Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;


figure;
plot(received_rc_inphase,received_rc_quadphase,'o');
title('Output constellation with memoryless noise (Raised cosine)');
xlabel('inphase');
ylabel('quadphase');
grid on;


%% line decode

% Memoryless
% Rect pulse and raised cosine

line_decoded_inphase = zeros(length(inphase),1);
line_decoded_quadphase = zeros(length(quadphase),1);
t_seq_rect = 10:19:length(received_quadphase)-1;
k=1;
for i = t_seq_rect
    line_decoded_inphase(k) = received_inphase(i);
    line_decoded_quadphase(k) = received_quadphase(i);
    k = k+1;
end

line_decoded_inphase_rc = zeros(length(inphase),1);
line_decoded_quadphase_rc = zeros(length(quadphase),1);
t_seq_rc = 10:19:length(received_rc_inphase)-1;
k=1;
for i = t_seq_rc
    line_decoded_inphase_rc(k) = received_rc_inphase(i);
    line_decoded_quadphase_rc(k) = received_rc_quadphase(i);
    k = k+1;
end

figure;%........................y2(t)
subplot(2,1,1);
stem(line_decoded_inphase);
title('y2(t) - Linedecoded inphase signal with noise ml (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(line_decoded_quadphase);
title('y2(t) - Linedecoded quadphase signal with noise ml(Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

figure;
subplot(2,1,1);
stem(line_decoded_inphase_rc);
title('y2(t) - Linedecoded inphase signal with noise ml(Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(line_decoded_quadphase_rc);
title('y2(t) - Linedecoded quadphase signal with noise ml (Raised Cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;


%% Decoder
%Memoryless
%Rect and raised cosine
decision_rect = zeros(length(line_decoded_inphase),1);
decision_rc = zeros(length(line_decoded_inphase_rc),1);

for i=1:length(line_decoded_inphase)
    decision_rect(i) = line_decoded_inphase(i) + (1i)*(line_decoded_quadphase(i));
end

for i=1:length(line_decoded_inphase_rc)
    decision_rc(i) = line_decoded_inphase_rc(i) + (1i)*(line_decoded_quadphase_rc(i));
end

final_decision_rect = decision_making(decision_rect);
final_decision_rc = decision_making(decision_rc);



%% Decimal Decoding
% Memoryless Channel

decimal_rect = decoding(final_decision_rect);
decimal_rc = decoding(final_decision_rc);


%% Decimal to Binary Conversion
% Memoryless Channel

binary_rect = dec2bin(decimal_rect,2) - '0';
binary_rc = dec2bin(decimal_rc,2) - '0';


%% DAC Conversion
% Memoryless Channel

final_binary_rect = reshape(binary_rect.',[],1);
final_binary_rc = reshape(binary_rc.',[],1);
figure;%....................y1(t)
subplot(2,1,1);
stem(final_binary_rect);
title('y1(t) - Digital to analog Conversion with noise ml (Rect pulse)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(final_binary_rc);
title('y1(t) - Digital to analog Conversion with noise ml (Raised cosine)');
xlabel('sample index');
ylabel('Amplitude');
grid on;

% Reshaping

reshaped_final_binary_rect = (reshape(final_binary_rect,8,[]))';
reshaped_final_binary_rc = (reshape(final_binary_rc,8,[]))';

% analog signal

analog_binary_rect = bi2de(reshaped_final_binary_rect)/255;
analog_binary_rc = bi2de(reshaped_final_binary_rc)/255;

sound(analog_binary_rect,Fs);

% [psd_rc_linecoded,~] = pwelch(rc_linecoded_signal,[],[],[],'centered',Fs);
% figure;
% semilogy(psd_rc_linecoded);
% xlabel('Frequency');
% ylabel('PSD');
% title('PSD of linecoded signal (raised cosine)');

function decision = decision_making(input_signal)
    decision = zeros(length(input_signal),1);
    for i = 1:length(input_signal)
        check = angle(input_signal(i));
        if(check > -pi/4 && check <= pi/4)
        decision(i) = 1 + (1i)*0;
        elseif(check > pi/4 && check <= 3*pi/4)
        decision(i) = 0 + (1i)*1;
        elseif(check > 3*pi/4 && check <= pi)
        decision(i) = -1 + (1i)*0;
        elseif(check > -pi && check <= -3*pi/4)
        decision(i) = -1 + (1i)*0;
        elseif(check > -3*pi/4 && check <= -pi/4)
        decision(i) = 0 + (1i)*(-1);
        end
    end
end

function bit_decision = decoding(input_signal)
         bit_decision = zeros(length(input_signal),1);
         for i = 1 : length(input_signal)
             check = input_signal(i);
             if(real(check) == 1 && imag(check) == 0)
                 bit_decision(i) = 0;
             elseif(real(check) == 0 && imag(check) == 1)
                 bit_decision(i) = 1;
             elseif(real(check) == -1 && imag(check) == 0) 
                 bit_decision(i) = 3;
             elseif(real(check) == 0 && imag(check) == -1)
                 bit_decision(i) = 2;
             end
         end
end



