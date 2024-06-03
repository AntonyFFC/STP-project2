T_T0 = [1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
K_K0_PID = [1.585 1.515 1.45 1.40 1.34 1.30 1.255 1.215 1.174 1.137 1.102];
K_K0_DMC = [1.95 1.93 1.80 0.87 0.35 0.12 0.05 0.05 0.04 0.04 0.04];

figure;
plot(T_T0, K_K0_PID, '-o')
grid on;
xlabel('T/T_0');
ylabel('K/K_0');
title('PID')
print ('pics/PID_ob_stab', '-dpng', '-r400')

figure;
plot(T_T0, K_K0_DMC, '-o')
grid on;
xlabel('T/T_0'); 
ylabel('K/K_0');
title('DMC')
print ('pics/DMC_ob_stab', '-dpng', '-r400')