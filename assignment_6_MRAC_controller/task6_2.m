

%% 

a = 0.02;
b = 1.3;
d = 10;


r = 35;



gamma1 = 5;
gamma2 = 1;
gamma3 = 1;

%% 
figure(1);
hold on;
grid on;
plot(x.time(1:2801),x.signals.values(1:2801));
plot(x_m.time(1:2801),x_m.signals.values(1:2801));
legend('model output', 'reference model output');