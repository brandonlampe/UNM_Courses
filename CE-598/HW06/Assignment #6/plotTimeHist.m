function plotTimeHist

    displacementsFlag = 'On';
    velocityFlag = 'Off';
    plastStretchFlag = 'Off';
    damageFlag = 'Off';
    IDFlag = 'Off'
    forcesIntFlag = 'Off'
    forceVsDispFlag = 'Off'

posCurPtr = 1;   %pointer to first position coordinate
posRefPtr = 3;   %pointer to first position coordinate
TimeHistfile = fopen('timehist.pdq','r');
dT = fscanf(TimeHistfile,'%g', 1);
particleID = fscanf(TimeHistfile,'%g', 1);
data = (fscanf(TimeHistfile, '%g', [4, inf]))';
numSteps = size(data, 1);
steps = [1 : numSteps];
times = steps*dT;
x = data(:, posCurPtr : posCurPtr + 1);
x_init = data(:, posRefPtr : posRefPtr + 1);
deformation = x - x_init;

fclose('all');
AxisFont = 8;
TextFont = 8;
LineWidth = .8;

if(  strncmpi(displacementsFlag, 'On', 2 ) )
    figname = ['plots of deformations of node #' num2str(particleID), ...
        [', initially located at ', num2str(x_init(1, 1)), '  ', num2str(x_init(1, 2))], ' , versus time'];
    scrsz = get(0,'ScreenSize');
    figure('name', figname, 'Position', [.25*scrsz(3) .25*scrsz(4) scrsz(3)/2 scrsz(4)/2]);
    
    subplot(2,1,1);
    set(gca,'FontSize',AxisFont);
    set(gcf,'Color','white');
    plot(times, deformation(:, 1),'LineWidth',LineWidth);
    ylabel('X - Disp','FontWeight','b','FontSize',TextFont);
    xlabel('Time, s','FontWeight','b','FontSize',TextFont);
    
    subplot(2,1,2);
    set(gca,'FontSize',AxisFont);
    plot(times, deformation(:, 2),'LineWidth',LineWidth);
    ylabel('Y - Disp','FontWeight','b','FontSize',TextFont);
    xlabel('Time, s','FontWeight','b','FontSize',TextFont);


end

return
