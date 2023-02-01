% Cell is discharged, find hte final temperature
clear
clc
%% Variables
% Cell Dimensions (eg 18650: D = 18mm, L = 65mm)
D = 0.018;
L = 0.065;

% Cell Capacity (Ah)
Ecell = 2.5;

% Cell Internal Resistance (Ω)
Rint = 0.06;

% Discharge Rate (C)
Crate = 2;

% Cooling rate (W/m^2K)
h = 12;

% Ambient Temperature
Ta = 25;

% Maximum Battery Temperature
Tmax = 80;

%% Calculations
% Current (A)
I = Ecell * Crate;

% Cylindrical Cell Area
A = 2*pi*(D/2)*L;

% Final Battery Temperature
Tbat = ((I^2*Rint)/(h*A)) + Ta;

% Increase in battery temperature
Tinc = Tbat - Ta;

%% Maintain temperature below Tmax
C = ((((Tmax-Ta)*h*A)/Rint)^0.5)/Ecell;

%% Output
msg = ['The final cell temperature is: ',num2str(Tbat)];
disp(msg);

msg1 = ['To keep the cell temperature below ',num2str(Tmax),', the maximum C rate is: ', num2str(C)];
disp(msg1)
