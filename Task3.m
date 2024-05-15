3 Task 3: Also, find at every SNR the value of α such that the
maximum SEP between UE1 and UE2 is minimized.
*The solution below, that we suggest, is with arrays. However, we know that it’s not the optimal
because of the delay and the inaccuracy of the code(because of the delay we had to decrease the
number of values of a so that the code run at a normal speed), but for some technical reasons
we couldn’t use the optimization toolbox to solve the task with optimization.

% -Define the parameters:
clc, clear, close all;

M = 4; % number of symbols in each constellation
N = 1000; % number of symbols to transmit
snr_dB = -10:0.5:20;
snr = 10.^(snr_dB/10);
a_values = linspace(0.001, 0.999, 100); % range of ’a’ values to search
max_SEP = zeros(length(snr_dB), length(a_values));
a_opt = zeros(1, length(snr_dB));

% -Find the optimal value of a:
for i = 1:length(snr_dB)
    for j = 1:length(a_values)
        a = a_values(j);

        Es1 = a;
        Es2 = (1 - a);
        Et = Es1 + Es2;
        Eg1 = sqrt(a/5);
        Eg2 = sqrt((1 - a)/5);

        s1 = randi([0 M-1], 1, N); % symbols for UE1
        s2 = randi([0 M-1], 1, N); % symbols for UE2

        x = Eg1 * pammod(s1, M) + Eg2 * pammod(s2, M) .* 1j;

        % Calculate noise standard deviation
        sigma_1 = sqrt(Et/(2*snr(i)));
        sigma_2 = sqrt(0.5*Et/(2*snr(i)));

        % Generate complex Gaussian noise
        noise1 = sigma_1 * (randn(1, N) + 1i * randn(1, N));
        noise2 = sigma_2 * (randn(1, N) + 1i * randn(1, N));

        % Add noise to the signals and demodulate
        y1 = x + noise1; % received signal at UE1
        y2 = x + noise2; % received signal at UE2
        r1 = real(y1); % demodulated symbols at UE1
        r2 = imag(y2); % demodulated symbols at UE2

        % Maximum Likelihood Detection
        reg1 = zeros(1, N);
        reg1(find(r1 < (-2 * Eg1))) = -3 * Eg1;
        reg1(find(r1 >= (2 * Eg1))) = 3 * Eg1;
        reg1(find(r1 >= (-2 * Eg1) & r1 < 0)) = -1 * Eg1;
        reg1(find(r1 >= 0 & r1 < (2 * Eg1))) = 1 * Eg1;

        reg2 = zeros(1, N);
        reg2(find(r2 < (-2 * Eg2))) = -3 * Eg2;
        reg2(find(r2 >= (2 * Eg2))) = 3 * Eg2;
        reg2(find(r2 >= (-2 * Eg2) & r2 < 0)) = -1 * Eg2;
        reg2(find(r2 >= 0 & r2 < (2 * Eg2))) = 1 * Eg2;

        % Calculate the Symbol Error Probabilities (SEP) for UE1 and UE2
        SEP1 = nnz(reg1 - real(x)) / N;
        SEP2 = nnz(reg2 - imag(x)) / N;

        % Store the maximum SEP between UE1 and UE2
        max_SEP(i, j) = max(SEP1, SEP2);
    end

    l = find(max_SEP(i, :) == min(max_SEP(i, :)), 1, 'first');
    a_opt(i) = a_values(l);
end

% -Plot ’optimal a’ vs SNR:
figure;
plot(snr_dB, a_opt);
xlabel('SNR (dB)');
ylabel('Optimal a');
title('Optimal a vs. SNR');
