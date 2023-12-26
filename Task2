Task 2: For System 1, by considering M = 4 and a = 0.25,
develop a simulation in MATLAB to plot the SEP vs. SNR for
both UE1 and UE2, where SNR = Et/No. Verify the simula-
tion results through the SEP expressions for M-PAM provided
in [2].

% -Initialize the variables:
clc, clear, close all;

M = 4;
N = 10000; % Number of symbols
a = 0.25; % Average Energy of S1
Eg1 = sqrt(a/5);
Eg2 = sqrt((1 - a)/5);
snr_dB = -10:0.5:20;
snr = 10.^(snr_dB/10);
num1_Err = zeros(1, length(snr));
num2_Err = zeros(1, length(snr));

% -Generate the symbols:
pam1 = [-3*Eg1, -1*Eg1, 1*Eg1, 3*Eg1]; % 4 - PAM1
pam2 = [-3*Eg2, -1*Eg2, 1*Eg2, 3*Eg2]; % 4 - PAM2
s = zeros(4, 4);

for i = 1:length(pam1)
    for j = 1:length(pam2)
        s(i, j) = pam1(i) + 1i * pam2(j);
    end
end

qam = reshape(s, 1, 16); % 16 - QAM
symbols_sent = randsrc(1, N, qam);
Et = norm(var(symbols_sent) - (mean(symbols_sent))^2);
s_real = real(symbols_sent);
s_imag = imag(symbols_sent);

% -Add Complex Gaussian Noise and compute the SEP:
for i = 1:length(snr)
    std1 = sqrt(Et/(2*snr(i))); % Standard Deviation 1
    std2 = sqrt(0.5*Et/(2*snr(i))); % Standard Deviation 2
    cgnoise1 = std1 * (randn(1, N) + 1i * randn(1, N));
    cgnoise2 = std2 * (randn(1, N) + 1i * randn(1, N));
    r1 = symbols_sent + cgnoise1;
    r2 = symbols_sent + cgnoise2;
    y1 = real(r1);
    y2 = imag(r2);

    % Decision regions
    reg1 = zeros(1, N);
    reg1(find(y1 < (-2 * Eg1))) = -3 * Eg1;
    reg1(find(y1 >= (2 * Eg1))) = 3 * Eg1;
    reg1(find(y1 >= (-2 * Eg1) & y1 < 0)) = -1 * Eg1;
    reg1(find(y1 >= 0 & y1 < (2 * Eg1))) = 1 * Eg1;

    reg2 = zeros(1, N);
    reg2(find(y2 < (-2 * Eg2))) = -3 * Eg2;
    reg2(find(y2 >= (2 * Eg2))) = 3 * Eg2;
    reg2(find(y2 >= (-2 * Eg2) & y2 < 0)) = -1 * Eg2;
    reg2(find(y2 >= 0 & y2 < (2 * Eg2))) = 1 * Eg2;

    num1_Err(i) = nnz(s_real - reg1); % number of errors for user 1
    num2_Err(i) = nnz(s_imag - reg2); % number of errors for user 1
end

SEP1 = num1_Err / N;
SEP2 = num2_Err / N;

% -Plot the SEP vs SNR for the 2 users:
figure
SEP1_theory = (M - 1) / M * erfc(sqrt(3 / (M^2 - 1) / 2 * snr));
SEP2_theory = (M - 1) / M * erfc(sqrt(3 / (M^2 - 1) / 2 * snr * (1 - a) / a));
semilogy(snr_dB, SEP1_theory, 'r', snr_dB, SEP2_theory, 'b')
hold on
grid on
semilogy(snr_dB, SEP1, 'y', snr_dB, SEP2, 'g');
xlabel('SNR(db)');
ylabel('SEP ');
legend('SEP1theory', 'SEP2theory', 'SEP1', 'SEP2');
title('SEP vs SNR(dB)');
