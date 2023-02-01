% Find the peak charging power and the OCV at which it occurs
clc
clear
%% Variables
% Pack Config
Ns = 96;
Np = 1;

% Cell Parameters
Vcellmax = 4.2;
Rint = 190e-3;
Icellmax = 200;

%% Pack Calculations
% Pack maximum voltage
Vpackmax = Ns * Vcellmax;

% Pack maximum current
Ipackmax = Np * Icellmax;

% Pack maximum power
Ppackmax = Vpackmax * Ipackmax;

% OCV it occurs at
Vdelta = Ipackmax * Rint;

OCV = Vpackmax - Vdelta;

%% Outputs
msg1 = ['The peak charging power is: ', num2str(Ppackmax), 'W'];
msg2 = ['This occurs at an OCV of  : ', num2str(OCV),'V'];

disp(msg1)
disp(msg2)
