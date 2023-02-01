% ======================================================
% Description:   Modified Track Data Analysis Function
% Author :       Peter Wilson, Jacob Wyborn
% Creation date: 08/12/2022
% Name:          trackData.m
% TASK 1 TEST CODE
% ======================================================

% ======================================================
% Step 0: define variables
% ======================================================

% Define the Track Data Name
trName = 'Silverstone.csv';

sim("test_car_2021a.slx",120)
% Define some constants relating to the car
mass = 300.0;



% grip between 0.0 (no grip) and 1.0 (max grip)
grip = 0.25;
% The greater kslip is the higher the sensitivity to slip angle
kslip = 2;

% define some analysis constants
% limit the slope change
dslopeMax = 5;

% The distance step between points on the track
dx = 5;

% ======================================================
% Step 2: Load in the track data into two arrays x and y
% ======================================================

% Load in the Data from the csv file as x,y data
trXY = readtable(trName);

% define columns of data for X and Y
x = table2array(trXY(:,1));
y = table2array(trXY(:,2));

% ======================================================
% Step 3: Plot the Track Layout (Racing Line)
% ======================================================

% Remember, each [x,y] point represents a point on track 5 m away from
% those next to it.
% plot(x,y)


% ======================================================
% Step 4: Analyse racing line data
% ======================================================

% calculate the gradient of the line
slope=diff(y)./diff(x);
% calculate the angle of travel
theta = atan(slope);

% calculate the difference in angle (rad)
slip = abs(diff(theta));
%ss = size(slip)
%ns = ss(1)

% ======================================================
% Step 5: Calculate the reduction in speed due to slip
% ======================================================

% The greater the slip, the greater the reduction in speed
slipMag = slip * kslip;
speedReduction = 1.0-min(slipMag,1.0-grip); %1

%duplicate the array of speedReduction
speedReductionConcat = speedReduction; 
speedReduction = cat(1,speedReduction, speedReductionConcat); 
speedReduction = cat(1,speedReduction, speedReductionConcat); 
speedReduction = cat(1,speedReduction, speedReductionConcat); 
speedReduction = cat(1,speedReduction, speedReductionConcat); 

%figure
%plot(speedReduction)

%Average the Speeds every 5 meters
dxAvSpeedArray = [];

count1 = 5; 
samples = 0; 
dxSpeedTotal = 0; 

for i=1:length(vehiclePos)
    if (count1 + 5) > 28980 
        break
    end

    if vehiclePos(i,2) <= count1 
        samples = samples + 1;
        dxSpeedTotal = dxSpeedTotal + vehicleSpeed(i,2);
    end
    
    if vehiclePos(i,2) > count1
            dxSpeedAverage = dxSpeedTotal/samples;
            dxAvSpeedArray = [dxAvSpeedArray, dxSpeedAverage];
            dxSpeedTotal = vehicleSpeed(i,2);
            count1 = count1 + 5;
            samples = 1;
    end
end



% ======================================================
% Step 6: Calculate the target speed based on maximum
% ======================================================
dxAvSpeedArray = transpose(dxAvSpeedArray);
LapSpeed = .*speedReduction;

figure;
plot(LapSpeed);
title('Vehicle Speed against Track Distance')
xlabel('Number of 5 Meter Steps') 
ylabel('Vehicle Speed (m/s)') 
writematrix(LapSpeed,'lapspeed.txt');


% ======================================================
% Step 7: Calculate Lap Times
% ======================================================

timeforStep = rdivide( dx , LapSpeed ); 
lapTime5 = sum(timeforStep, 'all'); 
disp('Lap time for 1 Lap = '); disp(lapTime5/5);
disp('Lap time for 5 Lap = '); disp(lapTime5);

% ======================================================
% Step 8: Calculate Energy Efficiency
% ======================================================

firstE = energy(1,2); %Initial energy of the battery system.

%Find the amount of energy used after 5 Laps
val = 0; 
for i=1:length(energy) 
    if energy(i,1) >= (lapTime5) && val == 0 
        lastE = energy(i,2); 
        val = 1;
    end
end

%Calculate the Energy Efficiency
Ediff  = firstE - lastE; %difference in start and finish energies
joulesPerkWh = 3.6e+06; % Joules converted to kWh
kWhPer5Laps = Ediff/joulesPerkWh; 
endDist = 29.455; %5 laps in km
KMPerkWh= endDist/kWhPer5Laps; %Efficiency in km/kWh
disp('Energy Effiency km/kWh:'); disp(KMPerkWh);

% ======================================================
% Step 9: Calculate Acceleration to 75 m
% ======================================================

val = 0; 
for i=1:length(vehicleSpeed)
    if vehiclePos(i,2) >= 75 && val == 0 
        timeat75 = vehiclePos(i,1); %Record that time for later analysis.
        speedat75 = vehicleSpeed(i,2);
        val = 1;
    end
    if val == 1 %Break out of loop when value is found to save time.
        break
    end
end

%Use the recorded time to calculate acceleration over 75m distance in g's
%and m/s^2.
a = calcAccel(0,75,timeat75); 
disp("Time at 75 meters =  ");  disp(timeat75);
disp("Speed at 75 meters (m/s) = "); disp(speedat75);

% Acceleration Times from 75m event

val = 0;
for i=1:length(vehicleSpeed)
    if vehiclePos(i,2) >= 75 && val == 0
        
    end
end
t = 4.05;
a = calcAccel(0,75,t);

% ======================================================
% Step N: define useful functions
% ======================================================

function ave = calculateAverage(x)
    ave = sum(x(:))/numel(x); 
end

function a = calcAccel(u,s,t)
   a = 2 * ( s - u*t ) / t^2;
end


