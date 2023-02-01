clear
clc
%% Variables

% Power (W)
P = 200;

% Ambient Temperature (C)
Ta = 40;

% Thermal Resistances (K/W)
Rthsa = 0.1;
Rthcs = 0.05;
Rthjc = 0.15;

%% Temperature Calculation
Ts = Rthsa * P + Ta;
Tc = Rthcs * P + Ts;
Tj = Rthjc * P + Tc; % For Si should be less than 175c
