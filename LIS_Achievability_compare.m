figure;
P_t = [-10 -5 0 5 10 15 20];
epsilon = [10^(-8) 10^(-7) 10^(-6) 10^(-5) 10^(-4) 10^(-3)];
A_success_JT=[
   0.033220000000000   0.034500000000000   0.035420000000000   0.038620000000000   0.041560000000000   0.044920000000000
   0.068300000000000   0.074120000000000   0.078360000000000   0.082720000000000   0.090140000000000   0.098800000000000
   0.153640000000000   0.160360000000000   0.174720000000000   0.184460000000000   0.198140000000000   0.220180000000000
   0.328080000000000   0.343520000000000   0.364900000000000   0.393240000000000   0.413500000000000   0.444680000000000
   0.616280000000000   0.642520000000000   0.662200000000000   0.691620000000000   0.719920000000000   0.756120000000000
   0.887880000000000   0.900860000000000   0.914640000000000   0.927120000000000   0.940200000000000   0.951680000000000
   0.988820000000000   0.989560000000000   0.992120000000000   0.995040000000000   0.995820000000000   0.997080000000000
];
A_success_OT=[
   0.029820000000000   0.032460000000000   0.034680000000000   0.036640000000000   0.040840000000000   0.044740000000000
   0.066880000000000   0.069340000000000   0.074600000000000   0.079160000000000   0.084260000000000   0.094140000000000
   0.135740000000000   0.142700000000000   0.151180000000000   0.160620000000000   0.172960000000000   0.187180000000000
   0.273240000000000   0.288960000000000   0.310280000000000   0.323260000000000   0.341480000000000   0.367680000000000
   0.510620000000000   0.529260000000000   0.547540000000000   0.575800000000000   0.598800000000000   0.626940000000000
   0.789580000000000   0.809100000000000   0.822980000000000   0.838780000000000   0.853900000000000   0.872300000000000
   0.957600000000000   0.963340000000000   0.968940000000000   0.974080000000000   0.977560000000000   0.982660000000000
];

subplot(2,1,1);
semilogx(epsilon, A_success_JT(4,:),'p', epsilon, A_success_JT(6,:),'p', epsilon, A_success_OT(4,:),'s', epsilon, A_success_OT(6,:),'s');
legend('P_T=5dBm, JT', 'P_T=15dBm, JT', 'P_T=5dBm, OT', 'P_T=15dBm, OT');
axis([10^-8 10^-4 0 1.1]);
xlabel('Error target(\epsilon)');
ylabel('Achievability');
grid on;

subplot(2,1,2);
plot(P_t, A_success_JT(:,3),'p', P_t, A_success_JT(:,6),'p', P_t, A_success_OT(:,3),'s', P_t, A_success_OT(:,6),'v');
legend('\epsilon=10^{-6}, JT', '\epsilon=10^{-3}, JT', '\epsilon=10^{-6}, OT', '\epsilon=10^{-3}, OT');
axis([-10 20 0 1]);
xlabel('Transmit power(dBm)');
ylabel('Achievability');
grid on;

sgtitle('Low interference setting')