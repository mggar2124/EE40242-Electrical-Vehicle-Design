clear
clc
%% Variables

syms Ld Lq Rs np RPM P Ke Iphase_rmsq 

% Stator Inductance
Ld = 80e-6;
%Lq = 200e-6;

% Stator Resistance
Rs = 0;

% Pole pairs
np = 1;

% Max Speed
RPM = 6000;

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

vd = -Wmax * Lq * iq;

vq = Wmax * lambdaF;

Vphase_peak = ((vd^2)+(vq^2))^0.5;

Vline_rms = (3/2)^0.5 * Vphase_peak;

%% ZDC Torque

% Max Speed in rad/s
iqZDC = Iphase_rmsq * 2^0.5;

TeZDC = 3/2 * np * lambdaF * iqZDC;

%%  MTPA Torque

ipeakMTPA = Iphase_rmsq * 2^0.5;

deltaI = lambdaF/((Ld-Lq));

K = deltaI/4;

idMTPA = -((((ipeakMTPA^2)/2)+K^2)^0.5) - K;

iqMTPA = ((ipeakMTPA^2)-(idMTPA^2))^0.5;

TeMTPA = 3/2 * np * (Ld-Lq) * idMTPA * iqMTPA + 3/2 * np * lambdaF * iqMTPA;

%%  ZDC RMS Motor Line Voltage

vdZDC = -Wmax * Lq * iqZDC;

vqZDC = Wmax * lambdaF;

VpeakZDC = ((vdZDC^2)+(vqZDC^2))^0.5;

Vline_rmsZDC = (3/2)^0.5 * VpeakZDC;

%% MTPA RMS Motor Line Voltage

vdMTPA = (idMTPA*Rs) - (Wmax*Lq*iqMTPA);

vqMTPA = (iqMTPA*Rs) + (Wmax*Ld*idMTPA) + (Wmax*lambdaF);

VpeakMTPA = ((vdMTPA^2)+(vqMTPA^2))^0.5;

Vline_rmsMTPA = (3/2)^0.5 * VpeakMTPA;

%% Print Outputs
msg1 = ['The Stator Phase Current (RMS) at Max Speed & Power: ', num2str(Iphase_rms), 'A'];
disp(msg1);

msg2 = ['The Line Voltage (RMS) at Max Speed & Power: ', num2str(Vline_rms), 'V'];
disp(msg2);

msg3 = ['The Torque under MTPA control at a Stator Phase Current of ', num2str(Iphase_rmsq), 'A (RMS):  ', num2str(TeMTPA), 'Nm'];
disp(msg3);

msg4 = ['The Torque under ZDC control at a Stator Phase Current of ', num2str(Iphase_rmsq), 'A (RMS):  ', num2str(TeZDC), 'Nm'];
disp(msg4);

msg4 = ['The Line Voltage (RMS) under MTPA control at a Stator Phase Current of ', num2str(Iphase_rmsq), 'A (RMS):  ', num2str(Vline_rmsMTPA), 'V'];
disp(msg4);

msg5 = ['The Line Voltage (RMS) under ZDC control at a Stator Phase Current of ', num2str(Iphase_rmsq), 'A (RMS):  ', num2str(Vline_rmsZDC), 'V'];
disp(msg5);

%% Solve
%{
variableToBeFound = RPM;

output = vpa(solve(Vline_rms, variableToBeFound));

termToBeFound = double(output(2));
%}
