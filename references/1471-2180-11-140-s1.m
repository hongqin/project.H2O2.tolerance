% scGrowthCurveSynchronization --- Program to perform growth curve
% synchronization from data obtained from diluting Pseudomonas aeruginosa
% PA14 WT cells.

% start by cleaning the matlab workspace of any variables and closing
% everything
clear all;
close all;

% create a new InfiniteM1000TimeSeries java object
myData = TimeSeriesData('101117-GcsWt32h.csv');

% set variables for the wavelengths. Here we have od600 and GFP
od600 = 1;
gfp   = 2;

% set the name  of samples - name, column in plate
myData = myData.setSampleName('2.5E-3', 2);
myData = myData.setSampleName('1.3E-3', 3);
myData = myData.setSampleName('6.3E-4', 4);
myData = myData.setSampleName('3.1E-4', 5);
myData = myData.setSampleName('1.6E-4', 6);
myData = myData.setSampleName('7.8E-5', 7);
myData = myData.setSampleName('3.9E-5', 8);
myData = myData.setSampleName('2.0E-5', 9);
myData = myData.setSampleName('blank', 10);

% Blank correction - subtract the median of the blank at each point
myData = myData.performBlankCorrectionPerLine('blank');

%import the sulfuric acid anthrone assay data
anthrone = importdata('Rhamnose.csv', ',');
% NOTE: the anthrone data was originally below zero for the last three
% samples, which were hardly going up in GFP RFUs, so the lowest value
% has been set to 0 and all the other values have been adjusted
% accordingly.

% NOTE 2: The first row is the rhamnose itself, the second row the positive
% error and the third row is the negative error.

% Calculate the medians, upper errors and lower errors for the GFP
% fluorescence, as well as upper, lower and median values for the rhamnose
% measurements.
for i = [1:8]
    y = [1:8]+(i)*8;
    medianODMatrix(:,i) = median(myData.wavelenghtData(1).data(:,y),2);
    maxODMatrix(:,i) = max(myData.wavelenghtData(1).data(:,y),[],2);
    minODMatrix(:,i) = min(myData.wavelenghtData(1).data(:,y),[],2);
    medianMatrix(:,i) = median(myData.wavelenghtData(2).data(:,y),2);
    maxMatrix(:,i) = max(myData.wavelenghtData(2).data(:,y),[],2);
    minMatrix(:,i) = min(myData.wavelenghtData(2).data(:,y),[],2);
    rhamnoseMatrix(:,i) = anthrone(1,i);
    rhamnoseMaxMatrix(:,i) = rhamnoseMatrix(:,i) + (anthrone(2,i));
    rhamnoseMaxError(:,i) = rhamnoseMaxMatrix(:,i) - rhamnoseMatrix(:,i);
    rhamnoseMinMatrix(:,i) = rhamnoseMatrix(:,i) - (anthrone(3,i));
    rhamnoseMinError(:,i) = rhamnoseMatrix(:,i) - rhamnoseMinMatrix(:,i);
end;

% Determine the lag phases between the different curves. The first number
% myDaya.optimizeTauArray defines the wavelength number on which the
% alignment will be based (1 = od600, 2 = gfp). With the determined
% timeShift the alignment is performed later on. The first graph, of the
% culture inoculated at 0.0025 od600, is considered to be 'time zero'.
tauArray = myData.optimizeTauArray(1, [1:8]);
timeShift = [0 tauArray];
timeMinusDelay = myData.wavelenghtData(2).times(:,end) - timeShift;

% The calculations that are necessary to make the graph of tau as a
% function of ln (X2/X1). The lnValues will be defined here, since they are
% common between experiments.
lnValues = [0.0025 0.0025/2 0.0025/4 0.0025/8 0.0025/16 0.0025/32 0.0025/64 0.0025/128];

for i= [1:8];
    lnValues(i) = [0.0025/lnValues(i)];
    lnValues(i) = [log(lnValues(i))];
end;

lnValues = lnValues';

tau = [timeShift];
tau = tau';
[b,bint,r,rint,stats] = regress (lnValues, [ones(size(tau(:,1))) tau(:,1)]);

% Make a subplot of 3-by-2 with OD600, OD600 synchronized, GFP, GFP
% synchronized, rhamnolipid secretion and tau vs ln.

% First figure: OD600 before alignment. In this graph the medians of the 8
% replicates will be plotted as a thick line and the entire range of the
% experiments will be plotted as thin lines in the same color.
figure(1);
set(gcf,...
    'Position', [600 731 1400 1600]);
subplot(3, 2, 1);
myData.plotManySamples(od600, [1:8]);
set(gca,...
    'YScale', 'log',...
    'YLim', [0.01 1],...
    'XLim', [0 35],...
    'FontSize', 16,...
    'Box', 'on',...
    'Color', 'none');
set(legend,...
    'Location', 'NorthWest',...
    'fontsize', 14,...
    'Color', 'none');
ylabel('OD_{600}',...
    'fontsize', 18);
xlabel('time (h)',...
    'fontsize', 18);
title('OD_{600} before alignment',...
    'fontsize', 22);

% Second figure: OD600 after alignment. The timeshift that is used to plot
% the various samples stems from the optimizeTauArray function. Plotted,
% again, are the medians in thick lines and the ranges of the experiments
% in thin lines.
subplot(3,2,2);
cmap = jet(8);

for i = 1:8;
    timeOnPlot = myData.wavelenghtData(2).times - timeShift(1, i);
    h = plot(timeOnPlot, medianODMatrix(:, i));
    hold on;
    set(h,...
        'Color', cmap(i, :),...
        'LineWidth', 2);
    h = plot(timeOnPlot, maxODMatrix(:, i));
    set(h,...
        'Color', cmap(i, :),...
        'LineWidth', 1);
    h = plot(timeOnPlot, minODMatrix(:, i));
    set(h,...
        'Color', cmap(i, :),...
        'LineWidth', 1);
end;
set(gca,...
    'XLim', [10 35],...
    'YScale', 'log',...
    'YLim', [0.01 1],...
    'FontSize', 16,...
    'Box', 'on',...
    'Color', 'none');
ylabel('OD_{600}',...
    'fontsize', 18);
xlabel('time (h)',...
    'fontsize', 18);
title('OD_{600} after alignment',...
    'fontsize', 22);
hold off;

% Third figure: GFP before alignment. As before, the median values of the 8
% replicates are plotted in thick lines and the range of the measurements
% in thin lines.
subplot(3, 2, 3);
myData.plotManySamples(gfp, [1:8]);
set(gca,...
    'YScale', 'lin',...
    'YLim', [0 30000],...
    'XLim', [0 35],...
    'FontSize', 16,...
    'Box', 'on',...
    'Color', 'none');
set(legend,...
    'Location', 'NorthWest',...
    'fontsize', 14,...
    'Color', 'none');
ylabel('GFP fluorescence (RFU)',...
    'fontsize', 18);
xlabel('time (h)',...
    'fontsize', 18);
title('GFP before alignment',...
    'fontsize', 22);

% Fourth figure: GFP after alignment. This alignment is done with the exact
% same time shifts as were used for the od600 alignment (timeOnPlot). The
% thick lines once again represent the medians and the thin lines the
% ranges of the measurements.
subplot(3, 2, 4);
cmap = jet(8);

for i = 1:8;
    timeOnPlot = myData.wavelenghtData(2).times - timeShift(1, i);
    h = plot(timeOnPlot, medianMatrix(:, i));
    hold on;
    set(h, 'Color', cmap(i, :), 'LineWidth', 2);
    h = plot(timeOnPlot, maxMatrix(:, i));
    set(h, 'Color', cmap(i, :), 'LineWidth', 1);
    h = plot(timeOnPlot, minMatrix(:, i));
    set(h, 'Color', cmap(i, :), 'LineWidth', 1);
end;
set(gca,...
    'XLim', [10 35],...
    'YScale', 'lin',...
    'YLim', [0 30000],...
    'FontSize', 16,...
    'Box', 'on',...
    'Color', 'none');
ylabel('GFP fluorescence (RFU)',...
    'fontsize', 18);
xlabel('time (h)',...
    'fontsize', 18);
title('GFP after alignment',...
    'fontsize', 22);
hold off;


% Fifth subplot: Rhamnolipid secretion. The entire range of the
% measurements is plotted as the error bars and the median values represent
% the squares.
subplot(3, 2, 5);
h = errorbar(timeMinusDelay,...
    rhamnoseMatrix,...
    rhamnoseMinError, ...
    rhamnoseMaxError,...
    's-',...
    'Color', 'black',...
    'LineWidth', 2);
set(gca,...
    'XLim', [10 35],...
    'YLim', [0 0.1],...
    'YTick', [0:0.025:0.1],...
    'FontSize', 16,...
    'Color', 'none');
ylabel('Rhamnose (g/L)',...
    'fontsize', 18);
xlabel('time (h)',...
    'fontsize', 18);
title('Rhamnolipid secretion',...
    'fontsize', 22);

% Sixth figure: tau as a function of ln (X2/X1). This graph should allow
% determination of the reproducibility of the lag phase, independent of the
% inoculum density, provided that the different points are on a straight
% line. The tau values, as determined by optimizeTauArray are plotted
% against the natural logarithm of the ratio of the two densities of the
% inocula. The points are then regressed by the regress function of the
% statistics toolbox in Matlab.
subplot(3, 2, 6);

h = plot(lnValues,...
    tau(:,1),...
    'sk',...
    'MarkerEdgeColor', 'none',...
    'MarkerFaceColor', [0 0 0],...
    'MarkerSize', 8);

hold on;

plot(lnValues, b(1)+(1/b(2))*lnValues,...
    '--',...
    'Color', [0 0 0]);

set(gca,...
    'XLim', [-0.5 5.5],...
    'YLim', [-2.5 22.5],...
    'YTick', [0:5:20],...
    'Color', 'none',...
    'FontSize', 16);
xlabel('ln (X_2/X_1)',...
    'fontsize', 18);
ylabel('\tau',...
    'fontsize', 18);
title('Reproducibility lag phase',...
    'fontsize', 22);

% This part is used to express the regression values and statistics on the
% last graph. A range is given for the error in determining the mu_max
% (maximum growth rate; slope = 1/mu_max), which is defined by the highest
% error compared to the average.
if stats(3)<0.0001;
    stats(3) = 0.0001;
end;
muError = [];
if b(2)-bint(2,1) > bint(2,2)-b(2);
    muError = [b(2)-bint(2,1)];
else;
    muError = [bint(2,2)-b(2)];
end;
str1(1) = {['R^2 = ', TimeSeriesData.num2strRound(stats(1), 3)]};
str1(2) = {['p < ', num2str(stats(3))]};
str1(3) = {['\mu_{max} = ', TimeSeriesData.num2strRound(b(2), 2),...
    '\pm', TimeSeriesData.num2strRound(muError(1), 2), ' h^{-1}']};
text(0.5, 19,...
    str1,...
    'verticalAlignment', 'top',...
    'FontSize', 16);
hold off;