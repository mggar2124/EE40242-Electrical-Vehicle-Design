clear
clc
%% Variables

Vstart = 3.2;

Vstop = 3.4;

Ipulse = 30;

%% Charge or Discharge

if Vstart - Vstop >= 0
    disp('The battery has been discharged.');
else
    disp('The battery has been charged.');
end

%% Internal Resistance of Battery

Rintt = abs((Vstart-Vstop)/Ipulse);

msg = ['The internal resistance of the battery is: ', num2str(Rintt) , 'Ω'];
disp(msg);