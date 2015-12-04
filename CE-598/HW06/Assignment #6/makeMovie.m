function makeMovie

plotDim = 2;
axisMarginFactor = 0.1;
trajectoryFlag = 'Off';
numbersFlag = 'Off';
linksInternalFlag = 'Off';
linksExternalFlag = 'Off';
dispReferenceFlag = 'On'
magFactor = 1.;
postProcfile = fopen('postProcParams.pdq','r');
if(postProcfile == -1)
    ['no postProcParams.pdq file found']
    return
end
postParams = fscanf(postProcfile,'%g');
fclose(postProcfile);
endTimeStep = postParams(1);
restartWriteInterval = postParams(2);
axisWindow = [-5 15 -15 2];
plotTimeSteps = [0 : restartWriteInterval : endTimeStep];
plotTimeSteps(1) = 1;
numFrames = size(plotTimeSteps, 2);
writerObj = VideoWriter('movie.avi');
writerObj.FrameRate = 8;
writerObj.Quality = 10;
open(writerObj);
scrsz = get(0,'ScreenSize');
figure('Position', [.1*scrsz(3) .1*scrsz(4) .8*scrsz(3) .8*scrsz(4)]);
set(gcf, 'color', 'white');
axis(axisWindow);
axis equal
axis tight
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