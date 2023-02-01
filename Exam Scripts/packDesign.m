% Design a battery pack to meet certain requirements
clear
clc
%% Variables

% Pack Energy (Wh)
ENreq = 50000;

% Pack Voltage Range (V)
Vpackminreq = 320;
Vpackmaxreq = 550;

% Vehicle Power (kW)
Pvehiclereq = 75000;

% Cell Energy (Wh)
Ecell = 200;

% Cell Voltage Range (V)
Vcellmin = 2.5;
Vcellmax = 4.2;

% Cell max discharge (A)
Imaxcell = 80;

%% Calculating Pack Configuration

% Power Constraints
Ns = ceil(ENreq/Ecell);

% Voltage Constraints
Nsmin = ceil(Vpackminreq/Vcellmin);
Nsmax = ceil(Vpackmaxreq/Vcellmax);

% Check Nsmin gives enough voltage to meet pack requirement

while Nsmin <= Nsmax
    if Vcellmin * Nsmin < Vpackminreq
        Nsmin = Nsmin + 1;
    else
        break
    end
end

% Number of parallel cells required
Np = 1;

while Np < 100
    if Nsmin * Np < Ns
        Np = Np + 1;
    else
        break
    end
end


% Second Higher Np

Np2 = 1;

while Np2 < 100
    if Vcellmin * Np2 * Imaxcell < Pvehiclereq
        Np2 = Np2 + 1;
    else
        break
    end
end



%% Calculating Pack Current, Voltage and Power

Ipack = Imaxcell * Np;

Vpackmin = Vcellmin * Nsmin;
Vpackmax = Vcellmax * Nsmin;

Ppackmin = Ipack * Vpackmin;
Ppackmax = Ipack * Vpackmax;

Epack = Nsmin * Np * Ecell ;

%% Outputs
msg = ['Pack Configuration   = ',num2str(Nsmin),'s',num2str(Np),'p'];
disp(msg)

msg1 = ['Pack Current         = ',num2str(Ipack)];
msg2 = ['Minimum Pack Voltage = ',num2str(Vpackmin)];
msg3 = ['Maximum Pack Voltage = ',num2str(Vpackmax)];
msg4 = ['Minimum Pack Power   = ',num2str(Ppackmin)];
msg5 = ['Maximum Pack Power   = ',num2str(Ppackmax)];

disp(msg1)
disp(msg2)
disp(msg3)
disp(msg4)
disp(msg5)

