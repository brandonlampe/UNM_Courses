function percentDifference = ...
    beamExample(LatRotAngle, numberOfMovieFrames, numTimeSteps)
% peridynamic plane stress beam example
if(nargin == 0)
  LatRotAngle = 0.;   %lattice rotation angle about Z-axis, in degrees
  numberOfMovieFrames = 20; % number of frames in the movie
  numTimeSteps = 2000;
end
loadFactor = -5.E4; %(times beam self-weight)
L = 0.01; % lattice spacing
density = 2323.0; % material density of concrete, kg/m^3
E = 24.86E9; % Young's modulus of concrete
dampRatio = 0.5; % external damping coefficient
beamSpan = 0.30; % span of beam between supports, m
beamLength = beamSpan + 4*L; % length of beam in X-direction, m
beamHeight = 0.06;  % width of beam in Y-direction, m
beamThickness = 0.03; % thickness of beam in Z-direction, m
writeInterval = floor(numTimeSteps/numberOfMovieFrames); % interval for writing restart file
PtclToMonitor = 1; % particle to monitor for time history output
history = zeros(numTimeSteps, 4);
tol = 0.0001;   % a small tolerance
volPtcl = (sqrt(3.)/2.)*beamThickness*(L^2); %Eq. 7.5 of Practical Peridynamics
mPtcl = density*volPtcl;
c = sqrt((3.)/2.)*E*L*beamThickness;  %see pg. 210 of Practical Peridynamics
dT = sqrt(mPtcl/c)/20.;
latticeOrigin = [beamSpan/2.; beamHeight/2.];
latBmatrix = setupLattice(L, LatRotAngle);
latticeRadius = 2.*sqrt((beamSpan/2.)^2 + (beamHeight/2.)^2);
nX = floor(latticeRadius/L);
nY = floor(latticeRadius/L);
%
% Insert the particles
%
figure('name', 'unmagnified deformed shape');
hold on
beamOverhang = (beamLength - beamSpan)/2;
numPtcls = 0;
for i = -nX : nX
    for j = -nY : nY
        xI = [i; j];
        xF = latBmatrix*xI + latticeOrigin;
        if(((xF(1) > -beamOverhang) && (xF(1) < beamSpan + beamOverhang) ) &&  ...
                ((xF(2) > tol) && (xF(2) < beamHeight - tol) ) )
            numPtcls = numPtcls + 1;
            Xref(numPtcls, :) = xF(:);
            BCcodes(numPtcls, 1:2) = 0;
            BCvals(numPtcls, 1:2) = 0;
            if((abs(xF(1)) < 0.5*L) && (xF(2) < L))
                BCcodes(numPtcls, 1:2) = 1; %pin the lower left-most particle
                leftSupport = numPtcls;
            elseif((abs(xF(1) - beamSpan) < L/2) && (xF(2) < L))
%                 BCcodes(numPtcls, 2) = 1; %roller at the lower right-most particle
                rightSupport = numPtcls;
            else
                BCvals(numPtcls, 2) = 1.;
            end
            if((i == 0) && (j == 0))
                PtclToMonitor = numPtcls;
            end
        end
    end
end
g = 9.81; % m/s^2
loadPerPtcl = loadFactor*mPtcl*g;
BCvals = BCvals*loadPerPtcl; 
xCur = Xref;      %initial positions of the particles
velocity = zeros(numPtcls, 2);
massBeam = mPtcl*numPtcls;
I = beamThickness*(beamHeight^3)/12;
deltaAnalytical = loadFactor*5*massBeam*g*(beamSpan^3)/(384*E*I)
fundAngularFreqBeam = sqrt((48*E*I/(beamSpan^3))/(massBeam/2));
fundPeriod = 2*pi/fundAngularFreqBeam;
dampCoeff = -2*mPtcl*dampRatio*fundAngularFreqBeam;
bondList = setupBondList(Xref, latBmatrix, numPtcls, L);
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
                internForce(iPtcl, :) = internForce(iPtcl, :) + dc*(c*stretch);
            end
        end
    end
    % Update the particle positions x:
    for iPtcl = 1 : numPtcls
        for iDir = 1 : 2
            if(BCcodes(iPtcl, iDir) == 0)
                externalForce = BCvals(iPtcl, iDir);
                dampForce = dampCoeff*velocity(iPtcl, iDir);
                forcePtcl = internForce(iPtcl, iDir) + externalForce ...
                    + dampForce;
                accel = forcePtcl/mPtcl;
                velocity(iPtcl, iDir) = velocity(iPtcl, iDir) + dT*accel;
                xCur(iPtcl, iDir)= xCur(iPtcl, iDir) ...
                    + dT*velocity(iPtcl, iDir);
            else
                xCur(iPtcl, iDir)= Xref(iPtcl, iDir) + BCvals(iPtcl, iDir);
                velocity(iPtcl, iDir) = 0.;
            end
        end
    end
    if(mod(iStep, writeInterval) == 0)
        plot(xCur(:, 1), xCur(:, 2),'Marker', ...
            'o','MarkerFaceColor', 'c','MarkerEdgeColor', 'none', ...
            'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 3);
        hold on
        plot(Xref(:, 1), Xref(:, 2),'Marker', ...
            'o','MarkerFaceColor', 'k','MarkerEdgeColor', 'none', ...
            'Marker', 'o', 'LineStyle', 'none', ...
            'LineWidth', 1,'MarkerSize', 1);
        plot(xCur(PtclToMonitor, 1), xCur(PtclToMonitor, 2),'Marker', ...
            'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'none', ...
            'Marker', 'o', 'LineStyle', 'none', ...
            'LineWidth', 1,'MarkerSize', 6);
        plot(xCur(leftSupport, 1), xCur(leftSupport, 2),'Marker', ...
            'o','MarkerFaceColor', 'b','MarkerEdgeColor', 'none', ...
            'Marker', 'o', 'LineStyle', 'none', ...
            'LineWidth', 1,'MarkerSize', 6);
        plot(xCur(rightSupport, 1), xCur(rightSupport, 2),'Marker', ...
            'o','MarkerFaceColor', 'b','MarkerEdgeColor', 'none', ...
            'Marker', 'o', 'LineStyle', 'none', ...
            'LineWidth', 1,'MarkerSize', 6);
        axis equal;
        axis([(min(Xref(:,1))-2*L) (max(Xref(:,1))+2*L) ...
             (min(Xref(:,2)) - beamSpan) (max(Xref(:,2))+3*L)]);
        hold off;
        title('unmagnified deformed shape');
        drawnow;
    end
    history(iStep, 1:2) = xCur(PtclToMonitor, :) - Xref(PtclToMonitor, :);
    if (iStep == 1)
        initial_xpos = xCur(PtclToMonitor, 1)
        initial_ypos = xCur(PtclToMonitor, 2)
    end
end
xpos = xCur(PtclToMonitor, 1)
ypos = xCur(PtclToMonitor, 2)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save as output %%%%%%
set(gcf,'paperunits','inches');
set(gcf,'papersize',[6,4]);
set(gcf,'paperposition',[0,0,6,4]);
fname = strcat('4_def');
fpath = 'C:\Users\Student\Documents\CE_598\Assignment #7\figs\';
saveas(gcf,strcat(fpath,fname), 'pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save as output %%%%%%
set(gcf,'paperunits','inches');
set(gcf,'papersize',[6,4]);
set(gcf,'paperposition',[0,0,6,4]);
fname = strcat('4_hist');
fpath = 'C:\Users\Student\Documents\CE_598\Assignment #7\figs\';
saveas(gcf,strcat(fpath,fname), 'pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%

deltaAnalytical = loadFactor*massBeam*g*5*(beamSpan^3)/(384*E*I)
deltaSPLM = history(numTimeSteps, 2)
percentDifference = 100*(deltaSPLM-deltaAnalytical)/deltaAnalytical

return
end

%----------------------------------------------------------------------------------------
function latB = setupLattice(L, LatRotAngle)
%----------------------------------------------------------------------------------------
%
% This function sets up the rotated lattice matrix latB for a hexagonal lattice
%
latRotZ = zeros(2);
LatRotAngle = LatRotAngle*pi/180.;  %convert to radians
latB(1,1) = L*1.;
latB(1,2) = L/2.;
latB(2,1) = 0.;
latB(2,2) = L*sqrt(3.)/2.;
latRotZ(1,1) = cos(LatRotAngle);
latRotZ(1,2) = -sin(LatRotAngle);
latRotZ(2,1) = sin(LatRotAngle);
latRotZ(2,2) = cos(LatRotAngle);
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
