clear
clc
%% Variables
% DC Voltage for Inverter
Vdc = 400;

% Target Line Voltage (V RMS)
Vline_rms = 220;

%% Modulation SPWM
Vline_peak = Vline_rms * 2^0.5;

Vphase_peak = Vline_peak/3^0.5;

Mspwm = Vphase_peak/(Vdc/2);

%% Modulation SVM

Msvm = Vline_peak/Vdc;

%% Difference
msg = ['It is clear that MSVM (',num2str(Msvm), ') < MSPWM (',(num2str(Mspwm)),') ,therefore using SVM (min-max injection) saves modulation ratio.'];

disp(msg);