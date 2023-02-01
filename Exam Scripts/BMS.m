% BMS Calculations
clear
clc
%% Variables
% Battery Voltage (V)
Vbatmax = 450;

% DC Link Capacitor (F)
Cdclink = 500e-6;

% DC Link Capacitor dV/dt limit (V/s)
Vlim = 100e6;

% Discharge Time (s)
Td = 5;

% Discharge Amount - percentage of maximum voltage
DC = 0.9;

%% Calculations
% Precharge resistor size
Rprecharge = Vbatmax/(Cdclink*Vlim);

% Discharge Resistor Size
Rdischarge = -Td/(log(1-DC)*Cdclink);

%% Outputs
msg1 = ['The minimum precharge resistor value is: ',num2str(Rprecharge), 'Ω'];
msg2 = ['The maximum discharge resistor value is: ', num2str(Rdischarge), 'Ω'];

disp(msg1)
disp(msg2)