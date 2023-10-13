clear;

data = load('ASP_Final_Data');
theta_s_noisy = data.theta_s_noisy;
theta_i_noisy = data.theta_i_noisy;
matX = data.matX;
L = numel(theta_s_noisy);  % Time length
t = [1:L];                 % Time vector
N = numel(matX(:,1));      % Number of isotropic antennas
d = 0:(N-1);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

%% Denoise Stage
M = 10;                    % Filter length
w_s = zeros(L-M+1,M);
w_i = zeros(L-M+1,M);
theta_s_hat = zeros(1,L);
theta_s_noisy = [zeros(1,M), theta_s_noisy];
theta_i_hat = zeros(1,L);
theta_i_noisy = [zeros(1,M), theta_i_noisy];

theta = linspace(-pi/2, pi/2, 500);
theta_deg = rad2deg(theta);
a = exp(j*pi*sin(theta).*d.');

mu = 0.0002;
for i = 1 : L
    theta_s_hat(:,i) = theta_s_noisy(i+M-1:-1:i)*w_s(i,:).';
    if(i == 1)
        w_s(i+1,:) = w_s(i,:);
    else
        w_s(i+1,:) = w_s(i,:) + mu * theta_s_noisy(i+M-1:-1:i)*(theta_s_noisy(i+M) - theta_s_hat(:,i));
    end
end

mu = 0.00003;
for i = 1 : L
    theta_i_hat(:,i) = theta_i_noisy(i+M-1:-1:i)*w_i(i,:).';
    if(i == 1)
        w_i(i+1,:) = w_i(i,:);
    else
        w_i(i+1,:) = w_i(i,:) + mu * theta_i_noisy(i+M-1:-1:i)*(theta_i_noisy(i+M) - theta_i_hat(:,i));
    end
end

%% Beamforming Stage

[y_array_steering, w_array_steering, B_array_steering] = array_steering_beamformer(matX, N, theta_s_hat, a, d, L);
[s_t_hat, w, B] = proposed_beamformer(matX, N, theta_s_hat, theta_i_hat, a, d, w_array_steering, L);

figure;
plot(t,theta_s_noisy(M+1:end), t,theta_i_noisy(M+1:end));
ylim([-10 20]);
legend('Signal $\tilde{\theta}_{s}(t)$','Interference $\tilde{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
title('DOA Without Denoised')

figure;
plot(t,theta_s_hat,t,theta_i_hat);
ylim([-10 20]);
legend('Signal $\hat{\theta}_{s}(t)$','Interference $\hat{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
title('DOA After Denoised')

figure;
subplot(2,1,1)
plot(t,real(s_t_hat))
% axis tight
ylim([-5.5 5.5])
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(s_t_hat))
% axis tight
ylim([-5.5 6.5])
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)

%%



%% RLS Method

% clear;
% 
% data = load('ASP_Final_Data');
% theta_s_noisy = data.theta_s_noisy.';
% theta_i_noisy = data.theta_i_noisy.';
% matX = data.matX;
% L = numel(theta_s_noisy);  % Time length
% t = [1:L];                 % Time vector
% N = numel(matX(:,1));      % Number of isotropic antennas
% d = 0:(N-1);
% 
% M = 10;                    % Filter length
% mu_s = 0.18;
% mu_i = 0.128;
% w_s = zeros(M,L-M+1);
% w_i = zeros(M,L-M+1);
% theta_s_hat = zeros(1,L);
% theta_i_hat = zeros(1,L);
% 
% theta = linspace(-pi/2, pi/2, 500);
% theta_deg = rad2deg(theta);
% a = exp(j*pi*sin(theta).*d.');
% 
% w = zeros(M,L-M+1);
% y = zeros(1,L);
% delta = 0.001;
% lambda = 0.992;
% P = eye(M)/delta;
% for i = 1 : L-M
%     k = (P*theta_s_noisy(i+M-1:-1:i))/(lambda + theta_s_noisy(i+M-1:-1:i).'*P*theta_s_noisy(i+M-1:-1:i));
%     xi = theta_s_noisy(i+M) - w(:,i).'*theta_s_noisy(i+M-1:-1:i);
%     w(:,i+1) = w(:,i) + k*xi;
%     theta_s_hat(:,i) = w(:,i+1).'*theta_s_noisy(i+M-1:-1:i);
%     P = P/lambda - k*theta_s_noisy(i+M-1:-1:i).'*P/lambda;
% end
% 
% for i = 1 : L-M
%     k = (P*theta_i_noisy(i+M-1:-1:i))/(lambda + theta_i_noisy(i+M-1:-1:i).'*P*theta_i_noisy(i+M-1:-1:i));
%     xi = theta_i_noisy(i+M) - w(:,i).'*theta_i_noisy(i+M-1:-1:i);
%     w(:,i+1) = w(:,i) + k*xi;
%     theta_i_hat(:,i) = w(:,i+1).'*theta_i_noisy(i+M-1:-1:i);
%     P = P/lambda - k*theta_i_noisy(i+M-1:-1:i).'*P/lambda;
% end
% 
% figure();
% hold on
% plot(theta_s_noisy(M+1:end))
% plot(theta_s_hat)
% figure();
% hold on
% plot(theta_i_noisy(M+1:end))
% plot(theta_i_hat)
% 
% 
% theta_s_noisy = data.theta_s_noisy;
% theta_i_noisy = data.theta_i_noisy;
% 
% smooth_lowess_s = smoothdata(theta_s_noisy,"lowess");
% noise_noisy_s = theta_s_noisy - smooth_lowess_s;
% n_power_noisy_s = var(noise_noisy_s);
% 
% s_power_s = sum(smooth_lowess_s.^2)/2000;
% SNR_noisy_s = 10*log10(s_power_s/n_power_noisy_s)
% noise_result_s = theta_s_hat - smooth_lowess_s;
% n_power_result_s = var(noise_result_s);
% SNR_result_s = 10*log10(s_power_s/n_power_result_s)
% smooth_rloess_i = smoothdata(theta_i_noisy,"rloess");
% noise_noisy_i = theta_i_noisy - smooth_rloess_i;
% n_power_noisy_i = var(noise_noisy_i);
% 
% s_power_i = sum(smooth_rloess_i.^2)/2000;
% SNR_noisy_i = 10*log10(s_power_i/n_power_noisy_i)
% 
% noise_result_i = theta_i_hat - smooth_rloess_i;
% n_power_result_i = var(noise_result_i);
% 
% SNR_result_i = 10*log10(s_power_i/n_power_result_i)