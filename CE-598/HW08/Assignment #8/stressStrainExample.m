function ...
    stressStrainExample(LatRotAngle, numberOfMovieFrames, numTimeSteps)
% peridynamic plane stress example
if(nargin == 0)
  LatRotAngle = 0.;   %lattice rotation angle about Z-axis, in degrees
  numberOfMovieFrames = 1; % number of frames in the movie
  numTimeSteps = 3;
end
strain = [0.001 0.000 0.0];
L = 0.01; % lattice spacing
density = 2323.0; % material density of concrete, kg/m^3
E = 24.86E9; % Young's modulus of concrete
planeStressFlag = 1;  % 1 = plane stress; 2 = plane strain
dampRatio = 0.5; % external damping coefficient
width = 0.05; % span of beam between supports, m
height = 0.05;  % width of beam in Y-direction, m
thickness = 0.03; % thickness of beam in Z-direction, m
writeInterval = floor(numTimeSteps/numberOfMovieFrames); % interval for writing restart file
PtclToMonitor = 1; % particle to monitor for time history output
history = zeros(numTimeSteps, 4);
tol = 0.0001;   % a small tolerance
volPtcl = (sqrt(3.)/2.)*thickness*(L^2); %Eq. 7.5 of Practical Peridynamics
mPtcl = density*volPtcl;
if(planeStressFlag == 1)
    c = sqrt((3.)/2.)*E*L*thickness;  %plane stress
    nu = 1/3;
else(planeStressFlag == 2)
    c = sqrt((3.)/2.)*E*L*thickness;  %plane strain
    nu = 1/4;
end
dT = sqrt(mPtcl/c)/20.;
latticeOrigin = [width/2.; height/2.];
latBmatrix = setupLattice(L, LatRotAngle);
latticeRadius = 2.*sqrt((width/2.)^2 + (height/2.)^2);
nX = floor(latticeRadius/L);
nY = floor(latticeRadius/L);
%
% Insert the particles
%
figure('name', 'unmagnified deformed shape');
hold on
numPtcls = 0;
for i = -nX : nX
    for j = -nY : nY
        xI = [i; j];
        xF = latBmatrix*xI + latticeOrigin;
        if(((xF(1) > 0) && (xF(1) < width) ) &&  ...
                ((xF(2) > tol) && (xF(2) < height - tol) ) )
            numPtcls = numPtcls + 1;
            Xref(numPtcls, :) = xF(:);
            BCcodes(numPtcls, 1:2) = 1; %applied displacement
            BCvals(numPtcls, 1:2) = strain(1)*xF(1);  %value
            BCcodes(numPtcls, 2) = 1; %applied displacement
            BCvals(numPtcls, 2) = strain(2)*xF(2); %value
            if((i == 1) && (j == 0))
                PtclToMonitor = numPtcls;
            end
        end
    end
end
% BCvals = BCvals*loadPerPtcl; 
xCur = Xref;      %initial positions of the particles
velocity = zeros(numPtcls, 2);
bondList = setupBondList(Xref, latBmatrix, numPtcls, L);
oppBond = [2 1 4 3 6 5];
for iStep = 1 : numTimeSteps
    curTime = iStep*dT;
    % Find internal forces acting upon each particle:
    internForce = zeros(numPtcls, 2);
    for iPtcl = 1 : numPtcls
        for iBond = 1 : 6
            if(bondList(iPtcl, iBond) ~= 0)
                dx = xCur(bondList(iPtcl, iBond), :) - xCur(iPtcl, :);
                curBondLength = sqrt(dx(1)^2 + dx(2)^2);
                dc = dx/curBondLength;
                stretch = (curBondLength - L)/L;
                internForce(iPtcl, :) = ...
                                    internForce(iPtcl, :) + dc*c*stretch;
            end
        end
    end
    % Update the particle positions x:
    for iPtcl = 1 : numPtcls
        for iDir = 1 : 2
            xCur(iPtcl, iDir)= Xref(iPtcl, iDir) + BCvals(iPtcl, iDir);
            velocity(iPtcl, iDir) = 0.;
        end
    end
    if(mod(iStep, writeInterval) == 0)
        cla;
        for iPtcl = 1 : numPtcls
            if(BCcodes(iPtcl, 1) + BCcodes(iPtcl, 2) > 0)
                markerColor = 'r';
                markerSize = 6;
            elseif(iPtcl == PtclToMonitor)
                markerColor = 'g';
                markerSize = 6;
            else
                markerColor = 'b';
                markerSize = 1;
            end
            line(Xref(iPtcl, 1), Xref(iPtcl, 2), ...
                'Marker', ...
                'o','MarkerFaceColor', 'k','MarkerEdgeColor', 'none', ...
                'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 2);
            line(xCur(iPtcl, 1), xCur(iPtcl, 2), ...
                'Marker', ...
                'o','MarkerFaceColor', markerColor,'MarkerEdgeColor', 'none', ...
                'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', markerSize);
        end
        plotBonds(xCur, bondList);
        axis equal;
        axis([(min(Xref(:,1))-2*L) (max(Xref(:,1))+2*L) ...
             (min(Xref(:,2)) - width/2.) (max(Xref(:,2))+3*L)]);
        title('unmagnified deformed shape');
        drawnow;
    end
    history(iStep, 1:2) = xCur(PtclToMonitor, :) - Xref(PtclToMonitor, :);
end

[stress strain] = computeStressStrain ...
          (xCur, Xref, L, LatRotAngle, bondList, E, nu, planeStressFlag);

displayTimeHist = 0;
if(displayTimeHist == 1)
    timeSteps = 0 : dT : (numTimeSteps-1)*dT;
    figure('name', 'Time History of monitored particle');
    hold on
    subplot(2,1,1);
    plot(timeSteps, history(:, 1));
    title('horizontal displacement of monitored particle versus time');
    xlabel('time');
    ylabel('horizontal displacement');
    subplot(2,1,2);
    plot(timeSteps, history(:, 2));
    title('vertical displacement of monitored particle versus time');
    xlabel('time');
    ylabel('vertical displacement');
end
stress
strain
return
end

%----------------------------------------------------------------------------------------
function latB = setupLattice(L, LatRotAngle)
%----------------------------------------------------------------------------------------
%
% This function sets up the rotated lattice matrix latB for a hexagonal lattice
%
LatRotAngle = LatRotAngle*pi/180.;  %convert to radians
latB = L*[1 1/2;
          0 sqrt(3)/2];
latRotZ = [cos(LatRotAngle) -sin(LatRotAngle);
           sin(LatRotAngle)  cos(LatRotAngle)];
latB = latRotZ*latB;
return
end

%----------------------------------------------------------------------------------------
function bondList = setupBondList(Xref, latBmatrix, numPtcls, L)
%----------------------------------------------------------------------------------------
%
% This function sets up the lattice bond list for a hexagonal lattice body
%
bondList = zeros(numPtcls, 6);
bond = zeros(6, 2);
tol = .1*L;
bond(1, :) =  (latBmatrix*[1; 0])';
bond(2, :) =  (latBmatrix*[-1; 0])';
bond(3, :) =  (latBmatrix*[0; 1])';
bond(4, :) =  (latBmatrix*[0; -1])';
bond(5, :) =  (latBmatrix*[-1; 1])';
bond(6, :) =  (latBmatrix*[1; -1])';
for iPtcl = 1 : numPtcls
    Xi = Xref(iPtcl, :);
    for jPtcl = 1 : numPtcls
        Xj = Xref(jPtcl, :);
        dX = Xj - Xi;
        bondLength = sqrt(dX(1)^2 + dX(2)^2);
        if(abs(bondLength - L) < tol)
            for iBond = 1 : 6
                dBonds = dX - bond(iBond, :);
                dist = sqrt(dBonds(1)^2 + dBonds(2)^2);
                if(dist < tol)
                    bondList(iPtcl, iBond) = jPtcl;
                end
            end
        end
    end
end
return
end

function plotBonds(x, bondList, bondDamage)
% plot the bonds in a lattice
numPtcls = size(x, 1);
numBonds = size(bondList, 2);
for iPtcl = 1 : numPtcls
    for jBond = 1 : numBonds
        jPtcl = bondList(iPtcl, jBond);
        if(jPtcl ~= 0)
            bondLine = [x(iPtcl, :); x(jPtcl, :)];
            
            jPtcl = bondList(iPtcl, jBond);
            bondLine = [x(iPtcl, :); x(jPtcl, :)];

            line(bondLine(:, 1), bondLine(:, 2), ...
                'Color', 'b', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
end
xlabel 'Link damage is indicated by colors of links'
return
end

function [stress strain] = computeStressStrain ...
            (xCur, Xref, L, LatRotAngle, bondList, E, nu, planeStressFlag)
numPtcls = size(xCur, 1);
numBonds = size(bondList, 2);
oppBond = [2 1 4 3 6 5];

D = (E/((1+nu)*(1-2*nu)))* ...
      [(1-nu)  nu   0          nu     0            0;
       nu    (1-nu) 0          nu     0            0;
       0      0     (1-2*nu)/2 0      0            0;
       nu     nu    0          (1-nu) 0            0;
       0      0     0          0      (1-2*nu)/2   0;
       0      0     0          0      0           (1-2*nu)/2];

stress = zeros(numPtcls, 6);
strain = zeros(numPtcls, 6);

%Fill in code here for the homework problem:   

for iPtcl = 1 : numPtcls
    
    
    
end
return
end
