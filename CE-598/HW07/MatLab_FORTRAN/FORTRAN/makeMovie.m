function makeMovie

plotDim = 2;
axisMarginFactor = 0.05;
trajectoryFlag = 'Off';
numbersFlag = 'Off';
linksInternalFlag = 'Off';
linksExternalFlag = 'Off';
dispReferenceFlag = 'Off'
magFactor = 1.0;
postProcfile = fopen('postProcParams.pdq','r');
if(postProcfile == -1)
    ['no postProcParams.pdq file found']
    return
end
postParams = fscanf(postProcfile,'%g');
fclose(postProcfile);
endTimeStep = postParams(1);
restartWriteInterval = postParams(2);
xmin = postParams(3) - axisMarginFactor;
xmax = postParams(4) + axisMarginFactor;
ymin = postParams(5) - axisMarginFactor;
ymax = postParams(6) + axisMarginFactor;
axisWindow = [xmin xmax ymin ymax];
plotTimeSteps = [0 : restartWriteInterval : endTimeStep];
plotTimeSteps(1) = 1;
numFrames = size(plotTimeSteps, 2);
writerObj = VideoWriter('movie.avi');
writerObj.FrameRate = 8;
writerObj.Quality = 10;
open(writerObj);
scrsz = get(0,'ScreenSize');
figure('Position', [.05*scrsz(3) .05*scrsz(4) .8*scrsz(3) .8*scrsz(4)]);
set(gcf, 'color', 'white');
axis(axisWindow);
axis equal
axis manual
% axis tight
for k = 1 : numFrames
    frame = k
    hold on
    h = plotDefShape(plotTimeSteps(k), magFactor, numbersFlag, ...
        dispReferenceFlag);
    frame = getframe;
    if(strncmpi(trajectoryFlag, 'Off', 2 ) && (k ~= numFrames))
      delete(h);
    end
    writeVideo(writerObj,frame);
end
close(writerObj);
fclose all
return
end