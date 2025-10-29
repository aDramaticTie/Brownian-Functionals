function [y,x]=createFit_FD(xf,marker,color)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(ZET_)
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

% This function was automatically generated on 02-Jul-2023 10:03:17

% Data from dataset "zet_ data":
%    Y = zet_

% Force all inputs to be column vectors
xf = xf(:);

% Prepare figure
figure(100)
hold on;
LegHandles = []; LegText = {};
% The following lines set log scales on both x and y axis.
% Comment where necessary
ax = gca;
%ax.YScale = 'log';
ax.XScale = 'log';


% --- Plot data originally in dataset "xf data"
[CdfF,CdfX] = ecdf(xf,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(xf,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = plot(BinCenter,BinHeight);
set(hLine,'LineStyle','none', 'Marker',marker,...
    'Markersize',7,'color',color);
xlabel('Data');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'zet_ data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% Adjust figure
box on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical',...
    'FontSize', 15, 'Location', 'northeast');
set(hLegend,'Interpreter','none');

y=BinHeight;
x=BinCenter;
