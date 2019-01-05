clc;
close all;


M = 10;
f = 1;
k = 9;

r = 10;



lambda0 = 1;
p0 = 1;

kp = 1;
km = 1;
sgn = sign(kp/km);

Gamma = [1 0 0 0;
         0 2 0 0;
         0 0 2 0;
         0 0 0 3];
     
 %% gains for r = 10
 %% gains for r = 10
lambda0 = 2;
p0 = 1;

kp = 1;
km = 1;
sgn = sign(kp/km);

Gamma = [1 0 0 0;
         0 1 0 0;
         0 0 1 0;
         0 0 0 1];
     
     
     %% plots
figure(1);
plot(y_p.time, y_p.signals.values)
title('comparison of process output and reference model output when r = 10');
hold on;
plot(y_m.time, y_m.signals.values);
legend('y_p', 'y_m');
xlabel('time');
grid on;

figure(2);
plot(y_p_2.time, y_p_2.signals.values)
title('comparison of process output and reference model output when r = 2sin(3t)+5sin(t)');
hold on;
plot(y_m_2.time, y_m_2.signals.values);
legend('y_p', 'y_m');
xlabel('time');
grid on;
