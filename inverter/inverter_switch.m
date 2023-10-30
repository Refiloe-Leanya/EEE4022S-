% SVPWM 
f_switching = 1e3;  % Switching frequency in Hz
f_fundamental = 50;  % Fundamental frequency in Hz
dt = 1 / f_switching;
t = 0:dt:1/f_fundamental;

% Modulation Index
m = 0.8;

% SVPWM Generation
u_alpha = m * sin(2 * pi * f_fundamental * t);
u_beta = m * cos(2 * pi * f_fundamental * t);

u_max = max(u_alpha);
u_min = min(u_alpha);

T1 = 0.5 * (1 + u_alpha / u_max);
T2 = 0.5 * (1 + u_beta / u_max);
T0 = 1 - T1 - T2;

% Display results
figure;
subplot(3,1,1);
plot(t, u_alpha, t, u_beta);
xlabel('Time');
ylabel('Voltage');
title('Control Voltages');

subplot(3,1,2);
plot(t, T1, t, T2, t, T0);
xlabel('Time');
ylabel('Duty Cycle');
title('SVPWM Switching Signals');

subplot(3,1,3);
stairs(t, 1*(T1 == 1), 'r');
hold on;
stairs(t, 2*(T2 == 1), 'g');
stairs(t, 3*(T0 == 1), 'b');
xlabel('Time');
ylabel('Switching State');
title('Switching Sequence');
