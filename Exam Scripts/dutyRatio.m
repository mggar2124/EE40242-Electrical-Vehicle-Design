% Calculate the battery current and power
clear
clc
%% Variables
% Duty Rate (%)
DR = 100;

% Battery Voltage (V)
Vbat = 528;

% Armature Resistance (Ω)
R = 7.5;

% Motor Speed (RPM)
RPM = 2521.6;

% Back EMF (V)
Vemf = 240.6;

% Desired Power (W)
Pdesire = 15000;

%% Calculations
% Effective Motor Voltage (V)
Vmotor = Vbat - Vemf;

% Battery Current (A)
Ibat = Vmotor/R;

% Battery Power (W)
Pbat = Vbat * Ibat;


%% Calculating Duty Rate

A = 1/R;
B = -Vemf/R;
C = -Pdesire;

Answer1 = ((-B)+((B^2-4*A*C))^0.5)/(2*A);
Answer2 = ((-B)-((B^2-4*A*C))^0.5)/(2*A);

if Answer1 >= 0
    Umotor = Answer1;
else
    Umotor = Answer2;
end

DRcalc = Umotor/Vbat;
DRcalcpercent = DRcalc * 100;

%% Outputs

msg1 = ['For a 100% DR the battery current is: ',num2str(Ibat),'A'];
msg2 = ['For a 100% DR the battery power is  : ', num2str(Pbat),'W'];
msg3 = ['To achieve a battery power of ',num2str(Pdesire),'W a ',num2str(DRcalcpercent),'% DR is required.'];

disp(msg1)
disp(msg2)
disp(msg3)

