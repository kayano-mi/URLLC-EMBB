
close all;
clear all;
clc;
% 理论值
% 参数
BS_number = 4;
P_t_dBm = 1;
P_t = 10.^(0.1 * P_t_dBm)/1000
noise_power_dBm =-97;
noise_power = 10.^(0.1 * noise_power_dBm)/1000


beta1 = 1;
c1 = 102000;
d1 = 73.22;
N1 = 4;
lambda1 = -8.2;
sigama_j = 0.00001;
epsilon = [10^-8 10^-7 10^-6 10^-5 10^-4];


gj = noise_power.*beta1./(noise_power.*beta1 + P_t.*sigama_j^2)
sj = noise_power.*beta1./(noise_power.*beta1 + 2.*P_t.*sigama_j^2)

gj_sum = prod(gj)
sj_sum = prod(sj)

a_B = gj_sum.* (gj_sum .* N1 - sj_sum - (N1 - 1) .* gj_sum .* gj_sum)./(sj_sum - gj_sum .* gj_sum)
b_B = (1 - gj_sum) .* a_B

b_x = exp(-lambda1/beta1)

B = beta(a_B, b_B)
B_x = betainc(b_x, a_B, b_B)

Achievability_AC_JT = B_x./B

semilogx(epsilon,Achievability_AC_JT);
legend('P_T=5dBm')
title('Low interference setting(JT)');
% axis([10^-8 10^-4 0 1.1]);
xlabel('Error target(\epsilon)');
ylabel('Achievability')
grid on;

