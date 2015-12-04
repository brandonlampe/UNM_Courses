function cableExample
% peridynamic cable example in 2D space
L = 1; % lattice spacing
M = 1; % mass of each particle
K = 100; % stiffness (force per unit stretch) of bond between particles.
dampCoeff = 0.05; % external damping coefficient
dT = sqrt(M/K)/3.; % time step
numTimeSteps = 1000; % number of time steps
numPtcls = 12; % number of lattice particles
X = zeros(numPtcls, 2); % reference positions of particles
X(:, 1) = L*(0 : numPtcls - 1); %initial horizontal positions of particles
X(:, 2) = 0.; % initial vertical positions of particles
BCcodes = zeros(numPtcls, 2); % boundary condition codes: 
                              %  0 = free to move
                              %  1 = applied displacement
BCcodes(1, 1:2) = 1; %fix the left-most particle
BCvals = zeros(numPtcls, 2); % applied force or displacement values.
BCvals(2:numPtcls, 2) = -1.; %apply a vertical force free particles
x = X; %initial positions of the particles
velocity = zeros(numPtcls, 2);
figure('name', 'Particle Positions');
for iStep = 1 : numTimeSteps
    % Find internal forces acting upon each particle:
    internForce = zeros(numPtcls, 2);
    for iBond = 1 : numPtcls - 1
        dx = x(iBond + 1, :) - x(iBond, :);
        curBondLength = sqrt(dx(1)^2 + dx(2)^2);
        dc = dx/curBondLength;
        stretch = (curBondLength - L)/L;
        internForce(iBond, :) = internForce(iBond, :) + dc*(K*stretch);
        internForce(iBond+1, :) = internForce(iBond+1, :) - dc*(K*stretch);
    end
    % Update the particle positions x:
    for iPtcl = 1 : numPtcls
        for iDir = 1 : 2
            if(BCcodes(iPtcl, iDir) == 0)
                externalForce = BCvals(iPtcl, iDir);
                dampForce = -dampCoeff*velocity(iPtcl, iDir);
                forcePtcl = internForce(iPtcl, iDir) + externalForce ...
                                                            + dampForce;
                accel = (forcePtcl)/M;
                velocity(iPtcl, iDir) = velocity(iPtcl, iDir) + dT*accel;
                x(iPtcl, iDir)= x(iPtcl, iDir) + dT*velocity(iPtcl, iDir);
            else
                x(iPtcl, iDir)= X(iPtcl, iDir) + BCvals(iPtcl, iDir);
                velocity(iPtcl, iDir) = 0.;
            end
        end
    end
    plot(x(:,1), x(:,2), 'Marker', 'o');
    axis(2*L*numPtcls*[-1, 1, -.5, .1]);
    drawnow;
    xHist(iStep) = x(numPtcls, 2); 
end
timeSteps = 0 : dT : (numTimeSteps-1)*dT;
figure('name', 'Time History');
plot(timeSteps, xHist);
title('vertical position of rightmost particle versus time');
xlabel('time');
ylabel('vertical position of righmost particle');
return
end

