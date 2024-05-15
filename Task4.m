Task 4: Consider a 4-QAM constellation for System 2. Also,
find at every SNR the value of θ such that the that the maxi-
mum SEP between UE1 and UE2 is minimized (i.e. the value
of θ that achieves user fairness) and plot SEP of UE1 and UE2
vs. SNR, where SNR = Et/No. Compare the performance of
System 1 and System 2.

*The solution below, that we suggest, is with arrays. However, we know that it’s not the optimal
because of the delay and the inaccuracy of the code(because of the delay we had to decrease
the number of values of theta so that the code run at a normal speed), but for some technical
reasons we couldn’t use the optimization toolbox to solve the task with optimization.

% -Initialize the parameters
clc, clear, close all;

N = 1000; % number of symbols sent
M = 2; % number of symbols in each constellation
a = 0.25; % average constellation energy of 2 -PAM
theta = 0.01 * pi:pi/100:pi; % angle of the 4 - QAM constellation

Eg1 = sqrt(3 * a);
Eg2 = sqrt(3 * (1 - a));
Et = 1; % average constellation energy for the 4 - QAM

% SNR
SNR_dB = -10:0.5:20;
SNR = 10.^(SNR_dB/10);

max_SEP = zeros(length(SNR), length(theta));
qam = zeros(M, M);
theta_opt = zeros(1, length(SNR));

% -Find the optimal value of theta for every SNR:
for i = 1:length(SNR)
    SEP1 = zeros(1, length(theta)); % matrix for SEP for UE1
    SEP2 = zeros(1, length(theta)); % matrix for SEP for UE2
    for j = 1:length(theta)

        % Generate the 4 -QAM constellation
        pam1 = pammod([0, M-1], M, theta(j));
        pam2 = pammod([0, M-1], M, theta(j));
        PAM1 = Eg1 * pam1;
        PAM2 = Eg2 * pam2;
        for k = 1:M
            for l = 1:M
                qam(k, l) = PAM1(k) + 1i * PAM2(l);
            end
        end
        QAM = reshape(qam, 1, 4);
        real_q = real(QAM);
        imag_q = imag(QAM);

        % The received signal before the noise
        x = randsrc(1, N, PAM1) + 1i * randsrc(1, N, PAM2);

        % Calculate noise standard deviation
        sigma_1 = sqrt(Et/(2 * SNR(i)));
        sigma_2 = sqrt(0.5 * Et/(2 * SNR(k)));

        % Generate complex Gaussian noise
        noise_1 = sigma_1 * (randn(1, N) + 1i * randn(1, N));
        noise_2 = sigma_2 * (randn(1, N) + 1i * randn(1, N));

        % Add noise to the signals and demodulate
        r1 = x + noise_1; % received signal UE1
        r2 = x + noise_2; % received signal at UE2
        y1 = real(r1); % demodulated symbols at UE1
        y2 = imag(r2); % demodulated symbols at UE2

        % Maximum Likelihood Detection
        d1 = norm(abs(real_q(4) - real_q(2))); % Euclidean Distance
        d2 = norm(abs(real_q(1) - real_q(2))); % Euclidean Distance
        reg1 = zeros(1, length(theta));
        if (real_q(4) > 0 & real_q(2) > 0)
            reg1(find(y1 > real_q(4) + d1/2)) = real_q(2);
            reg1(find(y1 > 0 & y1 < real_q(4) + d1/2)) = real_q(4);
            reg1(find(y1 < -(real_q(4) + d1/2))) = real_q(3);
            reg1(find(y1 < 0 & y1 > -(real_q(4) + d1/2))) = real_q(1);
        elseif real_q(4) < 0 & real_q(2) > 0
            reg1(find(y1 > min(real_q(1), real_q(2)) + d2/2)) = max(real_q(1), real_q(2));
            reg1(find(y1 < -(min(real_q(1), real_q(2)) + d2/2))) = -max(real_q(1), real_q(2));
            reg1(find(y1 > 0 & y1 < min(real_q(1), real_q(2)) + d2/2)) = min(real_q(1), real_q(2));
            reg1(find(y1 < 0 & y1 > -(min(real_q(1), real_q(2)) + d2/2))) = -min(real_q(1), real_q(2));
        else
            reg1(find(y1 > min(real_q(1), real_q(3)) + d1/2)) = max(real_q(1), real_q(3));
            reg1(find(y1 > 0 & y1 < min(real_q(1), real_q(3)) + d1/2)) = min(real_q(1), real_q(3));
            reg1(find(y1 < -(min(real_q(1), real_q(3)) + d1/2))) = -max(real_q(1), real_q(3));
            reg1(find(y1 < 0 & y1 > -(min(real_q(1), real_q(3)) + d1/2))) = -min(real_q(1), real_q(3));
        end
        d3 = norm(abs(imag_q(3) - imag_q(4))); % Euclidean Distance
        d4 = norm(abs(imag_q(1) - imag_q(3))); % Euclidean Distance
        reg2 = zeros(1, length(theta));
        if imag_q(3) > 0
            reg2(find(y2 > imag_q(3) + d3/2)) = imag_q(4);
            reg2(find(y2 > 0 & y2 < imag_q(3) + d3/2)) = imag_q(3);
            reg2(find(y2 < -(real_q(3) + d3/2))) = imag_q(1);
            reg2(find(y2 < 0 & y2 > -(imag_q(3) + d3/2))) = imag_q(2);
        elseif imag_q(3) < 0 & imag_q(4) > 0
            reg2(find(y2 > min(imag_q(2), imag_q(4)) + d4/2)) = max(imag_q(2), imag_q(4));
            reg2(find(y2 < -(min(imag_q(2), imag_q(4)) + d4/2))) = -max(imag_q(2), imag_q(4));
            reg2(find(y2 > 0 & y2 < min(imag_q(2), imag_q(4)) + d4/2)) = min(imag_q(2), imag_q(4));
            reg2(find(y2 < 0 & y2 > -(min(imag_q(2), imag_q(4)) + d4/2))) = -min(imag_q(2), imag_q(4));
        else
            reg2(find(y2 > min(imag_q(2), imag_q(1)) + d3/2)) = max(imag_q(1), imag_q(2));
            reg2(find(y2 > 0 & y2 < min(imag_q(2), imag_q(1)) + d3/2)) = min(imag_q(1), imag_q(2));
            reg2(find(y2 < -(min(imag_q(2), imag_q(1)) + d3/2))) = -max(imag_q(1), imag_q(2));
            reg2(find(y2 < 0 & y2 > -(min(imag_q(2), imag_q(1)) + d3/2))) = -min(imag_q(1), imag_q(2));
        end

        % Calculate the Symbol Error Probabilities for UE1 and UE2
        SEP1(j) = nnz(reg1 - real(x)) / N;
        SEP2(j) = nnz(reg2 - imag(x)) / N;

        % Store the maximum SEP between UE1 and UE2
        max_SEP(i, j) = max(SEP1(j), SEP2(j));
    end

    [~, l] = find(max_SEP(i, :) == min(max_SEP(i, :)), 1, 'first');
    theta_opt(i) = theta(l);
end

% -Plot ’optimal theta(rad)’ vs SNR:
figure;
plot(SNR_dB, theta_opt, '-o');
xlabel('SNR (dB)');
ylabel('Optimal theta');
title('Optimal theta vs. SNR');
