% Find the Thevenin Model Parameters
clear
clc

%% Variables
% Current Pulse Condition (A)
Ipulse = 30;

% Vertical Voltage rise (V)
Ur0 = 0.45;

% Total Voltage Rise (V)
Utotal = 0.6;

%Time to reach 98% of its final value (s)
Ts = 30;

%% Calculating Parameters
R0 = Ur0/Ipulse;

R1 = (Utotal/Ipulse)-R0;

C1 = Ts/(4*R1);

%% Print Outputs
msg1 = ['R0 = ', num2str(R0), 'Ω'];
msg2 = ['R1 = ', num2str(R1), 'Ω'];
msg3 = ['C1 = ', num2str(C1), 'F'];

disp(msg1);
disp(msg2);
disp(msg3);