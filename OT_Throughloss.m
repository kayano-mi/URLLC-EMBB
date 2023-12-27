close all;
clear all;
clc;

% % % Monte carmello仿真 
% 场景参数
simulation_num = 1000;
BS = [0 0; 0 300; 300 0; 300 300];
BS_number = 4;
URLLC_number = 1;
subband_number = 4;

% 功率、频率参数
P_t = [10 15 20 25 30]; %发射功率：dBm
% P_t = [10];
omega2 = -97;%子带噪声功率：dBm 
fc = 4e9; % carrier frequency in Hz

% 信道参数
shadow_mu = 1; % 对数正态lnX的均值
shadow_std = 7.8; % 对数正态lnX的标准差
shadow_mu_N = log((shadow_mu^2)/sqrt(shadow_std+shadow_mu^2)); % 对数正态对应的正态分布的均值
shadow_std_N = sqrt(log(shadow_std/(shadow_mu^2)+1)); % 对数正态对应的正态分布的标准差
Rayleigh_mu = 0.001;  % 瑞利信道造成的功率衰减值所服从的指数分布的均值

% MCS指标参数
lambda_m = -8.2;
c_m = 102000;
d_m = 73.22;   
MCS_parameter = struct;
MCS_parameter.lambda_m = [-8.2 -6.1 -4.64 -1.01 0.7 2.48 4.97 6.9 8.71 10.3 12.36 14.44];
MCS_parameter.c_m = [102000 197000 702000 313000 49700 522000 45000 46500 53400 15600 8770 4090 ];
MCS_parameter.d_m = [73.22 47.07 38.96 16.27 9.47 7.42 3.4 2.19 1.46 0.9 0.54 0.29];
MCS_parameter.beta_m = [1 1.4 1.4 1.48 1.5 1.62 3.1 4.32 5.37 7.71 15.5 19.6];
MCS_parameter.L_m = [3.078 2.052 1.714 1.539 1.368 1.197 1.026 0.855 0.684 0.52 0.342 0.171];
MCS_parameter_num = 1:length(MCS_parameter.beta_m);
MCS_parameter_num_flag = length(MCS_parameter.beta_m);


epsilon = [10^(-8) 10^(-7) 10^(-6) 10^(-5) 10^(-4) 10^(-3)];
distance_between_UE_URLLC = zeros(URLLC_number, BS_number);
pathloss_db = zeros(URLLC_number, BS_number); 
BLER = zeros(subband_number, 1);
A_success = zeros(length(epsilon),length(P_t));
L_m_sum = zeros(length(epsilon),length(P_t));

for P_t_flag = 1:length(P_t)
    for epsilon_flag = 1:length(epsilon)
        for simulation_num_flag = 1:simulation_num                
            %% 大尺度衰落      
            % pathloss 
            for car_number_flag = 1:URLLC_number
                for BS_number_flag = 1:BS_number
                   distance_between_UE_URLLC(car_number_flag,BS_number_flag) = pdist([[randi([0,300]),randi([0,300])];BS(BS_number_flag,:)],'euclidean');
                   pathloss_db(car_number_flag,BS_number_flag) = -(32.4+20.*log10(fc./1000000)+30.*log10(distance_between_UE_URLLC(car_number_flag,BS_number_flag)./1000));
                end
            end
            % 阴影衰落  
            shadow = lognrnd(shadow_mu_N, shadow_std_N, [URLLC_number, BS_number]);
            shadow_db = 10*log10(shadow);
           
            % 小尺度衰落  
            Rayleigh = exprnd(Rayleigh_mu, [subband_number, BS_number]);
            Rayleigh_db = 10*log10(Rayleigh);
    
            % 子带信噪比
%             Gamma_db = P_t(P_t_flag) + pathloss_db + shadow_db + Rayleigh_db - omega2; 
%             Gamma = 10.^(0.1*Gamma_db);
%             Gamma_BS = sum(Gamma,2);
%             Gamma_sub = exp(-Gamma_BS);
%             gamma = -log(0.25.*sum(Gamma_sub,1));
%             gamma_db = 10*log10(gamma)
            
            Gamma_db = P_t(P_t_flag) + pathloss_db + shadow_db + Rayleigh_db - omega2; 
            Gamma = 10.^(0.1*Gamma_db);
            Gamma = exp(-Gamma);
            Gamma_sub = sum(Gamma,1);
            gamma = -log(0.25.*Gamma_sub);
            gamma_db = 10*log10(gamma);
            

            % MCS选择
            for MCS_parameter_flag = flip(MCS_parameter_num)
                BLER = MCS_parameter.c_m(MCS_parameter_flag).*exp(-MCS_parameter.d_m(MCS_parameter_flag).*gamma);
                % 误块率计算     
                if (prod(BLER) <= epsilon(epsilon_flag)) || (MCS_parameter_flag == 1)
                    MCS_parameter_num_flag = MCS_parameter_flag;
                    L_m_sum(epsilon_flag, P_t_flag)  = L_m_sum(epsilon_flag, P_t_flag) + MCS_parameter.L_m(MCS_parameter_num_flag);
                    MCS_parameter_num_flag = length(MCS_parameter.beta_m);                
                    break;
                end
            end
        end
    end
end


L_m_sum = L_m_sum.';
L_m_sum = L_m_sum./simulation_num*60

figure(1);
semilogx(epsilon, L_m_sum(1,:),'-p');
legend('FSS and MCMSA:JT P_t=10dBm');
title('Low interference setting(OT)');
axis([10^-8 10^-3 20 140]);
xlabel('Error target(\epsilon)');
ylabel('Throughput Loss per PRB in kbps');
grid on;

figure(2);
plot(P_t, L_m_sum(:,4),'-p');
legend('FSS and MCMSA:JT \epsilon=10^{-5}');
title('Low interference setting(OT)');
axis([10 30 0 120]);
xlabel('Transmit power(dBm)');
ylabel('Throughput Loss per PRB in kbps');
grid on;




            
    