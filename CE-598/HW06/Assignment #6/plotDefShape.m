function h = plotDefShape(plot_timeStep, magFactor, numbersFlag, ...
         dispReferenceFlag)
if(nargin == 0)
    plot_timeStep = 1000;
    magFactor = 1.;  % magnification factor
    numbersFlag = 'On';
    dispReferenceFlag = 'On';
    postProcfile = fopen('postProcParams.pdq','r');
    postParams = fscanf(postProcfile,'%g');
    fclose(postProcfile);
    endTimeStep = postParams(1);
    restartWriteInterval = postParams(2);
    scrsz = get(0,'ScreenSize');
    figure('Position', [.1*scrsz(3) .1*scrsz(4) .8*scrsz(3) .8*scrsz(4)]);
    set(gcf, 'color', 'white');
    axis equal
end
postProcfile = fopen('postProcParams.pdq','r');
if(postProcfile == -1)
    ['no postProcParams.pdq file found']
    return
end
postParams = fscanf(postProcfile,'%g');
fclose(postProcfile);
endTimeStep = postParams(1);
restartWriteInterval = postParams(2);
h = [];
dispDeformedFlag = 'On';
%
% Open restart file
%
filename = sprintf('%s%d','Restart.', plot_timeStep,'.pdq');
Resfile = fopen(filename,'r');
currentTime = fscanf(Resfile,'%g',1);
numAlterDOF = 7;
curState = fscanf(Resfile,'%e',[numAlterDOF inf]);
fclose(Resfile);
numPtcls = size(curState, 2); % number of particles
if(numPtcls == 0)
    return
end
curPos = curState(4:5, 1:numPtcls);
refPos = curState(6:7, 1:numPtcls);
displ = curPos - refPos;
magPos(1:2,:) = refPos(1:2,:)+ magFactor*displ(1:2,:);
hold on;
for iPtcl = 1: numPtcls
    iGlID = curState(1, iPtcl);
    if( strncmpi(dispReferenceFlag, 'On', 2 ) )
        h1 = plot(refPos(1,iPtcl), refPos(2,iPtcl) ,'Marker', ...
            'o','MarkerFaceColor', 'k','MarkerEdgeColor', 'none', ...
            'LineStyle', 'none', 'LineWidth', 3,'MarkerSize', 6);
        h = [h h1];
    end
    if( strncmpi(numbersFlag, 'On', 2 ) )
        text(refPos(1,iPtcl), refPos(2,iPtcl), num2str(iGlID));
    end
    if(  strncmpi(dispDeformedFlag, 'On', 2 ) )
        
        h1 = plot(magPos(1,iPtcl), magPos(2,iPtcl),'Marker', ...
            'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'none', ...
            'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 6);
        h = [h h1];
    end
end
if(nargout == 0)
    h = [];
end
return
