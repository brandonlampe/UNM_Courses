function h = plotDefShape(plot_timeStep, magFactor, numbersFlag, ...
         dispReferenceFlag)
if(nargin == 0)
    plot_timeStep = 2000;
    magFactor = 1.;  % magnification factor
    numbersFlag = 'Off';
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
dispDeformedFlag = 'On';
postProcfile = fopen('postProcParams.pdq','r');
if(postProcfile == -1)
    ['no postProcParams.pdq file found']
    return
end
postParams = fscanf(postProcfile,'%g');
fclose(postProcfile);
endTimeStep = postParams(1);
restartWriteInterval = postParams(2);
ptclToMonitor = postParams(7);
h = [];
%
% Open restart file
%
filename = sprintf('%s%d','restart.', plot_timeStep,'.pdq');
Resfile = fopen(filename,'r');
currentTime = fscanf(Resfile,'%g',1);
numAlterDOF = 7;
curState = fscanf(Resfile,'%e',[numAlterDOF inf]);
fclose(Resfile);
numPtcls = size(curState, 2); % number of particles
if(numPtcls == 0)
    return
end
hold on;
for iPtcl = 1: numPtcls
    iGlID = curState(1, iPtcl);
    bcCodes = curState(2:3, iPtcl);
    curPos = curState(4:5, iPtcl);
    refPos = curState(6:7, iPtcl);
    displ = curPos - refPos;
    magPos = refPos + magFactor*displ;
    if( strncmpi(dispReferenceFlag, 'On', 2 ) )
        h1 = plot(refPos(1), refPos(2) ,'Marker', ...
            'o','MarkerFaceColor', 'k','MarkerEdgeColor', 'none', ...
            'LineStyle', 'none', 'LineWidth', 3,'MarkerSize', 6);
        h = [h h1];
    end
    if( strncmpi(numbersFlag, 'On', 2 ) )
        text(refPos(1), refPos(2), num2str(iGlID));
    end
    if(  strncmpi(dispDeformedFlag, 'On', 2 ) )
        if(iGlID == ptclToMonitor)
            h1 = plot(magPos(1), magPos(2),'Marker', ...
                'o','MarkerFaceColor', 'c','MarkerEdgeColor', 'none', ...
                'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 8);
        elseif((bcCodes(1) == 0) && (bcCodes(2) == 0))
            h1 = plot(magPos(1), magPos(2),'Marker', ...
                'o','MarkerFaceColor', 'g','MarkerEdgeColor', 'none', ...
                'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 6);
        else
            h1 = plot(magPos(1), magPos(2),'Marker', ...
                'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'none', ...
                'LineStyle', 'none', 'LineWidth', 1,'MarkerSize', 8);
        end
        h = [h h1];
    end
end
if(nargout == 0)
    h = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save as output %%%%%%
set(gcf,'paperunits','inches');
set(gcf,'papersize',[6,4]);
set(gcf,'paperposition',[0,0,6,4]);
fname = '5_def';
fpath = 'C:\Users\Student\Documents\CE_598\Assignment #7\figs\';
saveas(gcf,strcat(fpath,fname), 'pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%
return
