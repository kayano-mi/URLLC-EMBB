close all;
clear all;
clc;

% % % Monte carmello仿真 
% 场景参数
simulation_num = 50000;
BS = [0 0; 0 300; 300 0; 300 300];
BS_number = 4;
URLLC_number = 1;
subband_number = 4;

% 功率、频率参数
P_t = [-10 -5 0 5 10 15 20]; %发射功率：dBm
omega2 = -97;%子带噪声功率：dBm 
fc = 4e9; % carrier frequency in Hz

% 信道参数
shadow_mu = 1; % 对数正态lnX的均值
shadow_std = 7.8; % 对数正态lnX的标准差
shadow_mu_N = log((shadow_mu^2)/sqrt(shadow_std+shadow_mu^2)); % 对数正态对应的正态分布的均值
shadow_std_N = sqrt(log(shadow_std/(shadow_mu^2)+1)); % 对数正态对应的正态分布的标准差

Rayleigh_mu = 0.001;  % 瑞利信道造成的功率衰减值所服从的指数分布的均值

% 评价指标参数
lambda_m = -8.2;
c_m = 102000;
d_m = 73.22;    

epsilon = [10^(-8) 10^(-7) 10^(-6) 10^(-5) 10^(-4) 10^(-3)];
distance_between_UE_URLLC = zeros(URLLC_number, BS_number);
pathloss_db = zeros(URLLC_number, BS_number); 
BLER = zeros(subband_number, 1);
A_success = zeros(length(epsilon),length(P_t));

for P_t_flag = 1:length(P_t)
    for epsilon_flag = 1:length(epsilon)
        for simulation_num_flag = 1:simulation_num                
            %% 大尺度衰落      
            % pathloss 
            for car_number_flag = 1:URLLC_number
                for BS_number_flag = 1:BS_number
                   distance_between_UE_URLLC(car_number_flag,BS_number_flag) = pdist([randi([0,300]),[randi([0,300])];BS(BS_number_flag,:)],'euclidean');
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
            Gamma_db = P_t(P_t_flag) + pathloss_db + shadow_db + Rayleigh_db - omega2; 
            Gamma = 10.^(0.1*Gamma_db);
            Gamma = exp(-Gamma);
            Gamma_sub = sum(Gamma,1);
            gamma = -log(0.25.*Gamma_sub);
            gamma_db = 10*log10(gamma);
            
            % 误块率        
            for BLER_flag = 1:subband_number
                if gamma_db(BLER_flag) <= lambda_m
                    BLER(BLER_flag) = 1;
                else
                    BLER(BLER_flag) = c_m.*exp(-d_m.*gamma(BLER_flag));
                end
            end       
            % 可到达率
            A = prod(BLER);
            if A <= epsilon(epsilon_flag)
                A_success(epsilon_flag, P_t_flag) = A_success(epsilon_flag, P_t_flag) + 1;
            else
                A_success(epsilon_flag, P_t_flag) = A_success(epsilon_flag, P_t_flag);
            end    
        end
    end
end
A_success = A_success.';
A_success = A_success./simulation_num

figure(1);
semilogx(epsilon, A_success(1,:),'p', epsilon, A_success(2,:),'s', epsilon, A_success(3,:),'*', epsilon, A_success(4,:),'p', epsilon, A_success(5,:),'O', epsilon, A_success(6,:),'^', epsilon, A_success(7,:),'v');
legend('P_T=-10dBm', 'P_T=-5dBm', 'P_T=0dBm', 'P_T=5dBm', 'P_T=10dBm', 'P_T=15dBm', 'P_T=20dBm');
title('Low interference setting(OT)');
axis([10^-8 10^-3 0 1.1]);
xlabel('Error target(\epsilon)');
ylabel('Achievability');
grid on;

figure(2);
plot(P_t, A_success(:,1),'p', P_t, A_success(:,2),'s', P_t, A_success(:,3),'*', P_t, A_success(:,4),'p', P_t, A_success(:,5),'O', P_t, A_success(:,6),'v');
legend('\epsilon=10^{-8}', '\epsilon=10^{-7}', '\epsilon=10^{-6}', '\epsilon=10^{-5}', '\epsilon=10^{-4}', '\epsilon=10^{-3}');
title('Low interference setting(OT)');
axis([-10 20 0 1]);
xlabel('Transmit power(dBm)');
ylabel('Achievability');
grid on;

