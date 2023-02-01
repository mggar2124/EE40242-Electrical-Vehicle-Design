% ======================================================
% Description:   Track Data Analysis Function
% Author :       Matthew garcia
% Creation date: 20/12/2022
% Name:          TestCode3.m
% ======================================================
% Clear workspace and run simulink file

% need simscape, simscape electrical,

% clear
% clc

% Set number of laps
numLaps = 5;

% SIMULATION VARIABLES

% BATTERY Prismatic Pouch Cells
cellsSeries = 125;
CellStrings = 8;
cellVoltage = 3.6;
cellNomV = 3.3;
voltMax = cellVoltage * cellsSeries;
nominalV = cellNomV * cellsSeries;

cellCapacity = 4;
arrayCharge = cellCapacity * CellStrings;
cellR = 13.8e-3;
arrayR = cellsSeries*cellR;
packReistance = (1/arrayR)*CellStrings;
packResistance = 1/packReistance;



mQ = arrayCharge * 0.5;
mV = nominalV * 0.9;

voltageMin = 2.5 * cellsSeries;

% Gearbox
gearRatio = 13; %

% Motor emrax 228 Low Voltage CC
PMax = 75; % for continuous power
TMax = 130; %for continuous toreque
timeConstant = 0.01;
rotorInertia = 0.02521;
rotorDamping = 1e-5;
initRotorSp = 0;

speedVector = [0,500.5500,6500];
torqueVector = [100,140,140,90];

initSOC = 1;
batteryQ = arrayCharge * initSOC;


spEff = 2750;

motEff = 96;

torEff = 130;

% Mass
cellMass = 0.05;
nomMass = 100;
driverMass = 70;

inverterMass = 3;
batteryMass = cellMass*cellsSeries*CellStrings;
motorMass = 13.2;
mass = nomMass + driverMass + batteryMass + motorMass

tyreR = 0.24;

frontalArea = 0.9;
tyreRollingCoeff = 2e-6;
airDragCoeff = 0.85;

maxG = 29.43;

sim("myPMSMFinal.slx",100000) %Run model

Speed = Speed .* (5/18);

Position = tout .* Speed;

interpolatedXY = linspace(0, max(tout), max(tout)*10000);
interpolatedPosition = interp1(tout,Position,interpolatedXY);
interpolatedSpeed = interp1(tout,Speed,interpolatedXY);

% ======================================================
% Step 0: define variables
% ======================================================

% Define the Track Data Name
trackName = 'SilverstoneTrack.csv';
rlName = 'SilverstoneRaceline.csv';


% Define some constants relating to the car
Vmax = max(Speed);
gmax = 1.5;

% grip between 0.0 (no grip) and 1.0 (max grip), 0.25(base)
grip = 0.25;
% The greater kslip is the higher the sensitivity to slip angle
kslip = 2;

% Acceleration Times from 75m event
targetedDistance = 75;
[xmTime, xmSpeed] = xmData(targetedDistance, interpolatedXY, interpolatedPosition, interpolatedSpeed);

a = calcAccel(0,targetedDistance,xmTime);

% define some analysis constants
% limit the slope change
slopeMax = 5;

% The distance step between points on the track
dx = 5;

% ======================================================
% Step 2: Load in the track data into two arrays x and y
% ======================================================

% Load in the raceline data from the csv file as x,y data
rlXY = readtable(rlName);

% define columns of data for X and Y
xRl = table2array(rlXY(:,1));
yRl = table2array(rlXY(:,2));

% Make raceline laps times longer
xRl = repmat(xRl,numLaps,1);
yRl = repmat(yRl,numLaps,1);


%[interpXarray, interpYarray] = interpTrack(xRaceline, yRaceline);

% Do the same for track
trXY = readtable(trackName);
xTr = table2array(trXY(:,1));
yTr = table2array(trXY(:,2));

% ======================================================
% Step 3: Plot the Track Layout (Racing Line)
% ======================================================


% Sum of datapoints along tracklength
[trSeg,trSect] = findTrackLength(xTr, yTr);
trLength = trSeg(length(trSeg));

% Sum of datapoints along raceline
[rlSeg,rlSect] = findTrackLength(xRl, yRl);
rlLength = rlSeg(length(rlSeg));

% ======================================================
% Step 4: Analyse racing line data
% ======================================================

% calculate the gradient of the line
diffX = diff(xRl);
diffY = diff(yRl);

diffX(end+1) = xRl(1) - xRl(end);
diffY(end+1) = yRl(1) - yRl(end);

slope=diffY./diffX;

% calculate the angle of travel
theta = atan(slope);

% calculate the difference in angle (rad)
slip = abs(diff(abs(theta)));

slip(end+1) = abs(theta(1) - theta(end));

%ss = size(slip); %lateral slip
%ns = ss(1); %longitudinal slip

%%%%%%%%%%%%%%%%%
% MAX CORNER SPEED
%%%%%%%%%%%%%%%%%%%

%Finnd radius of curvature for each point in 1161

% Need the size of each segment ~ 5 and slip racelineParts

% size of seg(i) + size of seg(i+1) / slip(i)
% Gives R Radius of curvature
rlR = zeros(1, length(rlSect));
for ii = 1:length(rlSect)-1
    rlR(ii) = rlSect(ii+1)/slip(ii) + rlSect(ii);
end
rlR(end) = rlSect(1)/slip(end) + rlSect(end);
rlR = rlR.';


frictionalCoeff = 1.7;

% max cornner speed = square root of coefficient of friction of tyre and road for
% slip (not rolling) * R * gravity

rlMaxVel = sqrt(9.81*frictionalCoeff*rlR).';

% actual speed = max speed
rollSpeed = rlMaxVel;
% actualspeed(actaulspeed > Vmax) = Vmax;
rollSpeed(rollSpeed > Vmax) = Vmax;

standSpAv = rollSpeed;

initCorners = [];

cornerX = [];
cornerY = [];

steps = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initCorners = [];
cornerX = [];
cornerY = [];

for count1 = 1:length(theta) - steps
    prevTerms = theta(count1:count1+steps-1);
    avPrevterms = calculateAverage(prevTerms);

    if rlMaxVel(count1) <= Vmax
        initCorners = [initCorners, count1];
        cornerX = [cornerX, xRl(count1)];
        cornerY = [cornerY, yRl(count1)];
    end
end

if rlMaxVel(end) <= Vmax
    initCorners = [initCorners, length(theta)];
    cornerX = [cornerX, xRl(end)];
    cornerY = [cornerY, yRl(end)];
end

AccelerationVel =[];
accelerationXY  =[];
breakingXY      =[];
breakingXY      =[initCorners(1)];


for count1=1:length(initCorners)-1
    if initCorners(count1+1) > initCorners(count1)+10

        accelerationXY = [accelerationXY, initCorners(count1)];


        AccelerationVel = [AccelerationVel,rollSpeed(count1)];

    end

end

accelerationXY(end+1) = initCorners(end);
AccelerationVel(end+1) = rollSpeed(accelerationXY(end));

% Stop Point
for count1=2:length(initCorners)-1
    if initCorners(count1) > initCorners(count1-1)+10
        breakingXY = [breakingXY, initCorners(count1)];
    end
end

minSpeedList =[];
cornersMinVel = [];
accelList = [];

% Find min speed corner
minSpeedList = [];
cornersMinVel = [];
for count1=1:length(accelerationXY)
    minSpeedIdx = find(rollSpeed(breakingXY(count1):accelerationXY(count1)) == min(rollSpeed(breakingXY(count1):accelerationXY(count1))), 1);
    minSpeedList = [minSpeedList, (minSpeedIdx + breakingXY(count1)-1)];
    cornersMinVel = [cornersMinVel, rollSpeed(minSpeedList(count1))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ======================================================
% Step 5: Calculate the reduction in speed due to slip
% ======================================================

% The greater the slip, the greater the reduction in speed
slipMag = slip * kslip;
speedReduction = 1.0-min(slipMag,1.0-grip);

% ======================================================
% Step 6: Calculate the target speed based on maximum
% ======================================================

% Find the acceleration distance from start
[accelDist, maxSpeedRed] = accelD(Vmax, 0, Position, Speed);

% Find corresponding index in racelineSegments - gives distance to accelerate to
accelIdxTrack = find(rlSeg >= accelDist, 1, 'first');
accel_1 = find(Position >= accelDist, 1, 'first');

% Calculate the standing start points and their corresponding speed and distance values
standingStartPoints = linspace(1, accel_1, accelIdxTrack);
standingStartList = interp1(1:accel_1, Speed(1:accel_1), standingStartPoints);
standingStartDist = interp1(1:accel_1, Position(1:accel_1), standingStartPoints);


% Replace with standing start list
standSpAv(1:length(standingStartList)) = standingStartList;

% Accelerate out of every corner
for count1=1:length(minSpeedList)  

    [~,currentSpeedIdx] = min(abs(standingStartList-cornersMinVel(count1)));

    % Change it so that the minSpeedIdx isnt the choosing spot
    % Check if the end poinnt of array is at max speed on the max speed
    % graph  

    for j = 1:length(standingStartList)-currentSpeedIdx

        standSpAv( j + accelerationXY(count1) - 1) = standingStartList(j+currentSpeedIdx);
    end
end

% %Break into corners
for count1 = 1:length(minSpeedList)
    brakingDistance = breakDis(Vmax, cornersMinVel(count1), maxG);
    brakingSteps = round(minSpeedList(count1) - breakingXY(count1));
    brakingLine = linspace(standSpAv(breakingXY(count1) - brakingSteps - 1), cornersMinVel(count1), brakingSteps);
    standSpAv(breakingXY(count1) - brakingSteps : breakingXY(count1) - 1) = brakingLine;
end


standSpAv = [standSpAv, standSpAv(1:50)];
rlMaxVel = [rlMaxVel, rlMaxVel(1:50)];

initCorners = unique(initCorners);


for count1=initCorners
standSpAv(count1) = min(rlMaxVel(max(count1-30, 1):min(count1+30, end)));
end

rlMaxVel = rlMaxVel(1:length(xRl));
standSpAv = standSpAv(1:length(xRl));

% ======================================================
% Step 6: Calculate the lap speed
% ======================================================

rollT = findLapTime(rollSpeed,rlLength);
standingT = findLapTime(standSpAv,rlLength);

% Average Speed
rollSpAv = calculateAverage(rollSpeed);
standSpAv = calculateAverage(standSpAv);

% ======================================================
% Step 6: Calculate Efficienncy
% ======================================================

accelRat = rollSpAv / Vmax;

acceleratorT = [rollT(:),accelRat(:)];

% Find the index of racing line length in vehiclePosSIM
[~,rlEndj] = min(abs(Position-rlLength));

E_per_lap = Speed(rlEndj)/1000;
maxEnergy = Speed(1)/1000;
energy = Position(end)/1000;
batefficiency = energy/16.2*CellStrings;
energyLoss = maxEnergy - E_per_lap;
efficiency = (rlLength) / (energyLoss/(3600));

% ======================================================
% Print Some Stuff
% ======================================================
fprintf('Time to 75m = %.2fs\n', xmTime);
fprintf('Lap time for 1 Lap = %.2fs\n', standingT/5);
fprintf('Lap time for %d laps = %.2fs\n', numLaps, standingT);
fprintf('Efficiency of 1 lap = %.2fkm/kWh\n', batefficiency);
% ======================================================
% Step N: define useful functions
% ======================================================

function ave = calculateAverage(x)
    ave = sum(x(:))/numel(x);
end

function a = calcAccel(u,s,t)
   % s = ut + 0.5 a t**2
   % a = 2 * ( s - ut ) / t**2
   a = 2 * ( s - u*t ) / t^2;
end

% Calculate lap length
function [lapLength, segmentLength] = findTrackLength(x, y)
    % Calculate differences between consecutive elements of x and y
    segX = abs(diff(x));
    segY = abs(diff(y));

    % Append difference between first and last elements of x and y
    segX(end+1) = x(1) - x(end);
    segY(end+1) = y(1) - y(end);

    % Initialize segmentLength array
    segmentLength = zeros(size(segX));

    % Iterate over segX and segY
    for i = 1:length(segX)
        % Calculate segment length using Pythagorean theorem
        segmentLength(i) = sqrt((segX(i)^2)+(segY(i)^2));
    end

    % Calculate lap length as cumulative sum of segment lengths
    lapLength = cumsum(segmentLength);
end

% Calculate the laptime
function lapTime = findLapTime(lapSpeed, lapLength)

    lapTime = lapLength / calculateAverage(lapSpeed);

end

% Find time and speed at 75m
function [xmTime, xmSpeed] = xmData(targetDist, interpPoints, interpVehiclePos, interpVehicleSpeed)
    % Find index of element in interpVehiclePos closest to targetDist
    [~,idx] = min(abs(interpVehiclePos-targetDist));

    % Extract time and speed at that index from interpPoints and interpVehicleSpeed
    xmTime = interpPoints(idx);
    xmSpeed = interpVehicleSpeed(idx);
end


function [interpXarray, interpolatedYarray] = interpTrack(trackX, trackY)

    interpolatedYarray = [];
    interpXarray = [];

    for z = 1.0: 2: 1160
       
        x_n = trackX(z:z+1);
        y_n = trackY(z:z+1);
        interpolateedY = linspace(trackY(z),trackY(z+1),1000);
        interpedX = interp1(y_n,x_n,interpolateedY,'makima');
        interpolatedYarray = [interpolatedYarray, interpolateedY];
        interpXarray = [interpXarray, interpedX];

    end
end    

%Distance to accelerate to max speed from the current speed
function [accelDist, maxSpeedReduction] = accelD(maxSpeed, currentSpeed, vehiclePosArray, vehicleSpeedArray)
    maxSpeedReduction = maxSpeed * 0.8;
   
    % Find index of maxSpeed
    [~,maxSpeedIdx] = min( abs( vehicleSpeedArray - maxSpeedReduction));
   
    % Travelling at a current speed
    % Want to know how far it takes to accelerate to max speed
    maxSpeedDist = vehiclePosArray(maxSpeedIdx);

    % Current speed is in lapSpeed
    % Find index of current speed in vehicleSpeedArray
    [~,currentSpeed_1] = min(abs(vehicleSpeedArray-currentSpeed));

    % Distance at which you reach current speed
    currentSpeedDist = vehiclePosArray(currentSpeed_1);
   
    accelDist = maxSpeedDist - currentSpeedDist;
   
end

function breakingD = breakDis(initialV, finalV, maxG)

    breakingD = (finalV^2 - initialV^2)/(2*(maxG));
    breakingD = abs(breakingD);
end

function beginingCList = findStartCorner(theta)
   
    beginingCList = [];
   
    for i = 10:length(theta)
        prevX = theta(i-9:i-1);
        avX = calculateAverage(prevX);
       
        if i-avX >= 100

            beginingCList = [beginingCList,i];

        end
    end
end