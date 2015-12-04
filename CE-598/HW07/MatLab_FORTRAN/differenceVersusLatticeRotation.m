function differenceVersusLatticeRotation
% plot the difference between analytical and SPLM vertical deflection as
% as a function of lattice rotation angle
%
maxLatticeRotAngle = 30.;
numAngles = 7;
numTimeSteps = 4000;
numberOfMovieFrames = 20;

percentDifference = zeros(numAngles,1);
RotAngle = zeros(numAngles,1);
for i = 0 : numAngles - 1

% fill in one or two lines here    
    RotAngle(i+1) = maxLatticeRotAngle / (numAngles-1) * i;
    LatRotAngle = RotAngle(i+1)
    percentDifference(i+1) = beamExample ...
                     (LatRotAngle, numberOfMovieFrames, numTimeSteps);
        
% fill in one or two lines here                 
                 
end

figure('name', 'Percent Difference in Vertical Deflection Vs. Lattice Rotation');
hold on
plot(RotAngle, percentDifference);
title('Percent Difference in Vertical Deflection Vs. Lattice Rotation');
xlabel('Lattice Rotation, Degrees');
ylabel('Percent Difference in Vertical Deflection');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save as output %%%%%%
set(gcf,'paperunits','inches');
set(gcf,'papersize',[6,4]);
set(gcf,'paperposition',[0,0,6,4]);
fname = strcat('3_Rot-PerDiff');
fpath = 'C:\Users\Student\Documents\CE_598\Assignment #7\figs\';
saveas(gcf,strcat(fpath,fname), 'pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%

return
end

