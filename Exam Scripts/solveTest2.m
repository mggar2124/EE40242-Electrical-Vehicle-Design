clear
clc

%% Variables

syms Ld Lq Rs np RPM P Ke Iphase_rmsq 

%% Values
%Stator Inductance
Ld = 80e-6;
Lq = 200e-6;

% Stator Resistance
Rs = 0;

% Pole pairs
np = 1;

% Max Speed
%RPM = 6000;

% Power at Max Speed
P = 200e3;

% Back EMF constant
Ke = 0.036;

% Stator Phase Current (Given in Q)
Iphase_rmsq = 400;

%% Stator Phase Current & Line Voltage at max speed and power

% Max Speed in rad/s
Wmax = RPM * ((2*pi)/60);

% Torque at max speed
TeWmax = P/Wmax;

% Back EMF constant in Vs/rad
lambdaF = Ke*(60/(2*pi));

iq = ((2*TeWmax)/(3*np*lambdaF));

Iphase_rms = iq/2^0.5;

vd = -(Wmax * Lq * iq);

vq = Wmax * lambdaF;

Vphase_peak = ((vd^2)+(vq^2))^0.5 == 229.5063;

Vline_rms = (3/2)^0.5 * Vphase_peak ;


%% Solve
variableToBeFound = RPM;

output = vpa(solve(Vline_rms, variableToBeFound));

termToBeFound = double(output(2));






