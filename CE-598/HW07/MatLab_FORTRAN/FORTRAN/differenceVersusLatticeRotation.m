function differenceVersusLatticeRotation
% plot the difference between analytical and SPLM vertical deflection as
% as a function of lattice rotation angle
%
maxLatticeRotAngle = 30.;
numAngles = 7;
numTimeSteps = 4000;
for i = 0 : numAngles - 1

% fill in one or two lines here    
    
    percentDifference(i+1) = beamExample ...
                     (LatRotAngle, numberOfMovieFrames, numTimeSteps);
        
% fill in one or two lines here                 
                 
end

figure('name', 'Percent Difference in Vertical Deflection Vs. Lattice Rotation');
hold on
plot(rotAngle, percentDifference);
title('Percent Difference in Vertical Deflection Vs. Lattice Rotation');
xlabel('Lattice Rotation, Degrees');
ylabel('Percent Difference in Vertical Deflection');
return
end

