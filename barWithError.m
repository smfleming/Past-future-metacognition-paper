function barWithError(values,errors,width)
% function barWithError(values,errors,width)
%
% Similar to errobar except produces bar chart with error bars
% Col is line color (defaults to black)
% Face is face color (defaults to white; can be matrix)
%
% Uses Matlab's plot function
%
% SF 2012

numgroups = size(values, 1); % number of groups
numbars = size(values, 2); % number of bars in a group

% Plot bars
handles.bar = bar(values,width,'edgecolor','k', 'linewidth', 2);
if numbars == 1
    colormap white
else
    colormap gray
end
% Plot errors
x = handles.bar.XData;   % gives coordinates of centres
if ~isempty(errors)
    hold on
    for i = 1:numbars
        errorbar(x, values(:,i), errors(:,i),'k', 'linestyle', 'none', 'linewidth', 2);
    end
end