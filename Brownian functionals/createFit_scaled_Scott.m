function [y,x]=createFit_scaled_Scott(xf,sx,sy,marker,color)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(XF)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with distributionFitter
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  0
%
%   See also FITDIST.

% This function was automatically generated on 09-Sep-2021 19:54:47

% Data from dataset "xf data":
%    Y = xf

cc1=[0.49 0.18 0.56]; %purple
cc2=[0.3 0.75 0.93]; %light blue
cc3=[0.93 0.69 0.13]; %gold
% Force all inputs to be column vectors
xf = xf(:);

% Prepare figure
figure(1)
hold on;
LegHandles = []; LegText = {};
% The following lines set log scales on both x and y axis.
% Comment where necessary
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';


% --- Plot data originally in dataset "xf data"
[CdfF,CdfX] = ecdf(xf,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(xf,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = plot(BinCenter*sx,BinHeight*sy);
set(hLine,'LineStyle','none', 'Marker',marker,...
    'Markersize',7,'color',color);
xlabel('Data');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'xf data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% Adjust figure
box on;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical',...
    'FontSize', 15, 'Location', 'northeast');
set(hLegend,'Interpreter','none');

y=BinHeight*sy;
x=BinCenter*sx;