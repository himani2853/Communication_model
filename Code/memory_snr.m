lc;
clear;

[y, Fs] = audioread('project.wav');
y_mono = mean(y,2);

%% ADC conversion
M = 4;
levels = 2^8;
levels_found = linspace(min(y_mono), max(y_mono), levels); 
y_quantised = quantize_signal(y_mono,levels_found);
y_normalized = (y_quantised - min(y_quantised)) / (max(y_quantised) - min(y_quantised));
y_scaled = round((y_normalized)* 255);
y_check = de2bi(y_scaled);
y_binary = de2bi(y_scaled, 8);
y_plot = reshape(y_binary.',1,[]);

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

inphase_upsample = upsample(inphase, N);
quadphase_upsample = upsample(quadphase, N);


inphase_linecoded_signal = conv(inphase_upsample,rect_pulse);
quadphase_linecoded_signal = conv(quadphase_upsample,rect_pulse);

inphase_linecoded_signal = inphase_linecoded_signal(1:end-N);
quadphase_linecoded_signal = quadphase_linecoded_signal(1:end-N);


%% b) raised cosine
a = 1;
m = 9;
length_rc = 1;
[transmit_filter, ~] = raised_cosine(a,m,length_rc);

inphase_rc_linecoded_signal = conv(inphase_upsample,transmit_filter);
quadphase_rc_linecoded_signal = conv(quadphase_upsample,transmit_filter);

inphase_rc_linecoded_signal = inphase_rc_linecoded_signal(1:end-N);
quadphase_rc_linecoded_signal = quadphase_rc_linecoded_signal(1:end-N);


%% modulation
%rect pulse

fc = 1e6;
t_rect = 0:1/(10*fc):(length(inphase_linecoded_signal)-1)/(10*fc);
A = 1;
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

final_modulated_rect = modulated_rect_inphase + (1i)*modulated_rect_quadphase;
final_modulated_rc = modulated_rc_inphase + (1i)*modulated_rc_quadphase;


%% Memoryless AWGN Channel
EbN0_dB = 0:1:10;
EbN0_linear = 10.^(EbN0_dB / 10); 
Es = sum(transmit_filter.^2);

Eb = Es/2;
N0 = Eb./EbN0_linear;
j = sqrt(-1);

sigma_square = N0/2;

BER = zeros(size(EbN0_dB));

for i = 1 : length(EbN0_dB)
j = sqrt(-1);
noise_rect = sqrt(sigma_square(i))*(randn((size(final_modulated_rect)))+ j*(randn(size(final_modulated_rect)))); 
noise_rc = sqrt(sigma_square(i))*(randn((size(final_modulated_rc)))+ j*(randn(size(final_modulated_rc))));


%% AWGN Channel with Memory
del = zeros(length(rect_pulse)+1,1);
a = 0.8;
del(1) = a;
del(end) = (1-a);
rect_mem = conv(final_modulated_rect,del);
received_signal_rect_mem = rect_mem(1:end-length(rect_pulse)) + noise_rect;

rc_mem = conv(final_modulated_rc,del);
received_signal_rc_mem = rc_mem(1:end-length(transmit_filter)) + noise_rc;

%% Demodulation
%% Channel with Memory
%Rect pulse

received_signal_rect_cos_mem = 2*received_signal_rect_mem.*cos(2*pi*fc*t_rect)';
received_signal_rect_sine_mem = 2*received_signal_rect_mem.*sin(2*pi*fc*t_rect)';

received_inphase_mem = lowpass(real(received_signal_rect_cos_mem),5500,Fs);
received_quadphase_mem = lowpass(imag(received_signal_rect_sine_mem),5500,Fs);

received_signal_rc_cos_mem = 2*received_signal_rc_mem.*cos(2*pi*fc*t_rc)';
received_signal_rc_sine_mem = 2*received_signal_rc_mem.*sin(2*pi*fc*t_rc)';

received_inphase_rc_mem = lowpass(real(received_signal_rc_cos_mem),5500,Fs);
received_quadphase_rc_mem = lowpass(imag(received_signal_rc_sine_mem),5500,Fs);


%% line decode
% With Memory
% Rect pulse and Raised Cosine

line_decoded_inphase_mem = zeros(length(inphase),1);
line_decoded_quadphase_mem = zeros(length(quadphase),1);

t_seq_rect = 10:19:length(received_quadphase_mem)-1;
k=1;
for j = t_seq_rect
    line_decoded_inphase_mem(k) = received_inphase_mem(j);
    line_decoded_quadphase_mem(k) = received_quadphase_mem(j);
    k = k+1;
end

line_decoded_inphase_rc_mem = zeros(length(inphase),1);
line_decoded_quadphase_rc_mem = zeros(length(quadphase),1);
t_seq_rc = 10:19:length(received_inphase_rc_mem)-1;
k=1;
for j = t_seq_rc
    line_decoded_inphase_rc_mem(k) = received_inphase_rc_mem(j);
    line_decoded_quadphase_rc_mem(k) = received_quadphase_rc_mem(j);
    k = k+1;
end

%% Decoder
% with  memory
% rect and raised cosine

decision_rect_mem = zeros(length(line_decoded_inphase_mem),1);
decision_rc_mem = zeros(length(line_decoded_inphase_rc_mem),1);

for j=1:length(line_decoded_inphase_mem)
    decision_rect_mem(j) = line_decoded_inphase_mem(j) + (1i)*(line_decoded_quadphase_mem(j));
end

for j=1:length(line_decoded_inphase_rc_mem)
    decision_rc_mem(j) = line_decoded_inphase_rc_mem(j) + (1i)*(line_decoded_quadphase_rc_mem(j));
end

final_decision_rect_mem = decision_making(decision_rect_mem);
final_decision_rc_mem = decision_making(decision_rc_mem);

%% Decimal Decoding
% With Memory

decimal_rect_mem = decoding(final_decision_rect_mem);
decimal_rc_mem = decoding(final_decision_rc_mem);


%% Decimal to Binary Conversion
% with memory

binary_rect_mem = dec2bin(decimal_rect_mem,2) - '0';
binary_rc_mem = dec2bin(decimal_rc_mem,2) - '0';

%% DAC Conversion
% with memory

final_binary_rect_mem = reshape(binary_rect_mem.',[],1);
final_binary_rc_mem = reshape(binary_rc_mem.',[],1);

% Reshaping

reshaped_final_binary_rect_mem = (reshape(final_binary_rect_mem,8,[]))';
reshaped_final_binary_rc_mem = (reshape(final_binary_rc_mem,8,[]))';

% analog signal

analog_binary_rect_mem = bi2de(reshaped_final_binary_rect_mem)/255;
analog_binary_rc_mem = bi2de(reshaped_final_binary_rc_mem)/255;

error = sum(abs(y_plot(:) - final_binary_rc_mem(:)));
disp(error);
BER(i) = (error / numel(y_plot));
disp(BER(i));

end

figure;
semilogy(EbN0_dB, (BER), '-o');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate (BER) vs Eb/N0 for raised cosine');

grid on;

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



