% ======================================================
% Description:   Track Data Analysis Function
% Author :       Matthew Garcia
% Creation date: 20/12/2022
% Name:          TestCode2.m
% ======================================================


clear
clc

% Set number of laps
%laps = 1
numLaps = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim("test_car_2021.slx",120)

interpolatedXY = linspace(0,max(tout),max(tout)*10000);
interpolatedPosition = interp1(tout,Position,interpolatedXY);
interpolatedSpeed = interp1(tout,Position,interpolatedXY);

% ======================================================
% Step 0: define variables
% ======================================================

% Define the Track Data Name
rlName = 'SilverstoneRaceline.csv';
trackName = 'SilverstoneTrack.csv';

% Define some constants relating to the car
Vmax = max(Speed);
maxG = 1.5;

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


%%%%%%%%%%%%%%%%%
% MAX CORNER SPEED
%%%%%%%%%%%%%%%%%%%

%Finnd radius of curvature for each point in 1161

% Need the size of each segment ~ 5 and slip racelineParts

% size of seg(i) + size of seg(i+1) / slip(i)
% Gives R Radius of curvature
for count1=1:length(rlSect)-1
    rlR(count1) = (rlSect(count1+1)/slip(count1)+rlSect(count1));
end

rlR(end+1) = (rlSect(1)/slip(end)+rlSect(end));
rlR = rlR.';


frictionalCoeff = 1.75;

% max cornner speed = square root of coefficient of friction of tyre and road for
% slip (not rolling) * R * gravity
rlMaxVel = sqrt(frictionalCoeff*rlR*9.81);
rlMaxVel = rlMaxVel.';

% actual speed = max speed
rollSpeed = rlMaxVel;
% actualspeed(actaulspeed > Vmax) = Vmax;
rollSpeed(rollSpeed > Vmax) = Vmax;

standSpAv = rollSpeed;

initCorners = [];

cornerX = [];
cornerY = [];

steps = 10;


count1 = steps;    
while count1 <= length(theta)
    prevTerms = theta(count1-steps+1:count1-1);
    avPrevterms = calculateAverage(prevTerms);
  

    if rlMaxVel(count1) <= Vmax
        initCorners = [initCorners,count1];
        cornerX = [cornerX, xRl(count1)];
        cornerY = [cornerY, yRl(count1)];       

        count1 = count1 + 1;
    else
        count1 = count1 + 1;

    end
end


accelerationXY=[];
AccelerationVel = [];
breakingXY= [];
breakingXY = [initCorners(1)];


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
for count1=1:length(accelerationXY)

    [~,minSpeedIdx] = min(rollSpeed(breakingXY(count1):accelerationXY(count1)));
    minSpeedList = [minSpeedList, (minSpeedIdx + breakingXY(count1)-1)];
    cornersMinVel = [cornersMinVel, rollSpeed(minSpeedList(count1))]; 

end

% ======================================================
% Step 5: Calculate the reduction in speed due to slip
% ======================================================

% The greater the slip, the greater the reduction in speed
slipMag = slip * kslip;
speedReduction = 1.0-min(slipMag,1.0-grip);

% ======================================================
% Step 6: Calculate the target speed based on maximum
% ======================================================

% Find the acceleratio distance from start
[accelDist , maxSpeedRed] = accelD(Vmax, 0, Position, Speed);

% Find corresponding index in racelineSegmennts - gives distannce to
% accelerate too
[~,accelIdxTrack] = min(abs(rlSeg-accelDist));
[~,accel_1] = min(abs(Position-accelDist));

standingStartPoints = 1:accel_1/accelIdxTrack:accel_1;
standingStartList = interp1(1:accel_1,Speed(1:accel_1),standingStartPoints);
standingStartDist = interp1(1:accel_1,Position(1:accel_1),standingStartPoints);


% Replace with standing start list
for count1 = 1:length(standingStartList)
    standSpAv(count1) = standingStartList(count1);
end

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
for count1=1:length(minSpeedList)

    brakingDistance = breakDis(Vmax,cornersMinVel(count1), maxG);
    brakingSteps = round(minSpeedList(count1)-breakingXY(count1));

    brakingLine = linspace(standSpAv(breakingXY(count1)-brakingSteps-1),cornersMinVel(count1),brakingSteps);

    for j = 1:length(brakingLine)

        standSpAv(breakingXY(count1)-brakingSteps+j-1) = brakingLine(j);

    end

end 


standSpAv = [standSpAv, standSpAv(1:50)];
rlMaxVel = [rlMaxVel, rlMaxVel(1:50)];

for count1=1:length(initCorners)
    
    if initCorners(count1+1) ~= initCorners(count1)+1 && initCorners(count1)+10 >= initCorners(count1+1)
        initCorners = [initCorners(1:count1), initCorners(count1)+1 ,initCorners(count1+1:end)];

    end
end


for count1=initCorners

    standSpAv(count1) = min(rlMaxVel(count1-30:count1+30));
    
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
% rlLength = rlLength;
E_per_lap = Speed(rlEndj)/1000;
maxEnergy = 7943633;
% energy = Position(end)/1000
energyLoss = maxEnergy - E_per_lap;
efficiency = (rlLength) / (energyLoss/(3600));

% ======================================================
% Print Some Stuff
% ======================================================
fprintf('Time to 75m = %.2fs\n', xmTime);
fprintf('Lap time for 1 Lap = %.2fs\n', standingT/5);
fprintf('Lap time for %d laps = %.2fs\n', numLaps, standingT);
fprintf('Efficiency of 1 lap = %.2fkm/kWh\n', efficiency);
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

function [interpXarray, interpolatedYarray] = interpTrack(trackX, trackY) %function for interpolated track values

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


%Distance get to max speed from the current speed
function [accelDist, maxSpeedReduction] = accelD(maxSpeed, currentSpeed, vehiclePosArray, vehicleSpeedArray)
    maxSpeedReduction = maxSpeed * 0.8;
    
    % Find index of maxSpeed
    [~,maxSpeedIdx] = min( abs( vehicleSpeedArray - maxSpeedReduction));

    maxSpeedDist = vehiclePosArray(maxSpeedIdx);

    [~,currentSpeed_1] = min(abs(vehicleSpeedArray-currentSpeed));

    % Distance at which you reach current speed
    currentSpeedDist = vehiclePosArray(currentSpeed_1);
    
    accelDist = maxSpeedDist - currentSpeedDist;
    
end

function breakingD = breakDis(initialV, finalV, maxG) %function to determine breaking operations

    breakingD = (finalV^2 - initialV^2)/(2*(maxG));
    breakingD = abs(breakingD);
end

function beginingCList = findStartCorner(theta) %creates list of corner starting positions
   
    beginingCList = [];
   
    for i = 10:length(theta)
        prevX = theta(i-9:i-1);
        avX = calculateAverage(prevX);
       
        if i-avX >= 100

            beginingCList = [beginingCList,i];

        end
    end
end

