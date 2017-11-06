%PLOT_STATES   Plot States
%   PLOT_STATES(t1,x1,t2,x2,....,tn,xn,states, options) creates a subplot 
%   for each state in the same figure, with n different lines per subplot.
%
% Options:
%
% - title       figure title                string
% - ylabel      ylabel per subplot          cell of strings
% - legend      legend for each subplot     cell of strings
% - fontsize    fontsize                    integer
% - linewidth   linewidth                   integer
% - bgcolor     background color            array of rgb values
% - events      mark events with lines      array of time instants
%
% Example:
%
% plot_states(t1,x1,t2,x2,[1 5])
% shows the first column of x1 and x2 in the first subplot and the fifth
% column of x1 and x2 in the second subplot.
%
% Erik Vlasblom
% Last modified: 12-06-2015

% ============== PLOT STATES ===================

function PlotStates(varargin)

% ----- READ OPTIONS
% defaults
fontsize = 15;
linewidth = 2;
bgcolor = [1 1 1];
clrs = {'k','g','b','r','c','m','y','w'};


% get options
opt_provided = false;
if isa(varargin{end},'struct')
    options_user = varargin{end};
    opt_provided = true;
else
    options_user = struct();
end

% get states
if ~opt_provided
    st = varargin{end};
else
    st = varargin{end-1};
end
if (size(st,1)~=1 && size(st,2)~=1)
    error('The states should be provided as a vector')
end

% get subplots
plts = length(st);
if opt_provided
    endplt = nargin - 2;
else
    endplt = nargin - 1;
end

% set options
plot_options = struct('title','',...
    'ylabel',[],...
    'legend',[],...
    'fontsize',fontsize,...
    'linewidth',linewidth,...
    'bgcolor',bgcolor,...
    'events',[]);
setbyuser = fieldnames(options_user);
for i = 1:length(setbyuser)
    plot_options.(setbyuser{i}) = options_user.(setbyuser{i});
end
events_present = false;
if ~isempty(plot_options.events)
    events_present = true;
end




% ----- PLOT
scrsz = get(0,'ScreenSize');
try ti = plot_options.title{1};
catch; ti = plot_options.title;
end
StatefigHandle = figure('Name',ti,'OuterPosition',[1 scrsz(4)*(1/8) scrsz(3)/2 scrsz(4)*(13/16)],...
    'Color',bgcolor);

% For each state
for i = 1:plts 
    ax(i) = subplot(plts,1,i);
    
    % plot lines
    for j = 1:2:endplt
        plot_next(varargin{j},varargin{j+1},st(i),clrs{ceil(j/2)});
        hold on
    end
    
    % plot events
    if events_present
        events = plot_options.events;
        for ii = 1:length(events)
            line([events(ii),events(ii)],ylim,'color','k','linewidth',1);
        end
    end
    hold off
    
    % title 
    if i == 1
        pt = title(plot_options.title);
        set(pt,'fontsize',fontsize);
    end
    
    % ylabel
    if ~isempty(plot_options.ylabel)
        ylabel(plot_options.ylabel(i))
    else
        ylabel(['State ',num2str(st(i))])
    end
    
    % xlabel
    if i == plts
        xlabel('Time (s)')
    end
    
    % legend
    if isempty(plot_options.legend)
        for ii = 1:endplt/2;
            plot_options.legend{ii} = ['State ',num2str(st(ii))];
        end
    end
    if events_present
        plot_options.legend{end+1} = 'Events';
    end
    if ~isempty(plot_options.legend) && i == 1
        legend(plot_options.legend)
    end    
end
linkaxes(ax,'x');

% ---------- SUBFUNCTIONS ------------------------------------
    function plot_next(tin,xin,currst,c)
        if (length(tin) > 1) && (length(xin) > 1)
            plot(tin,xin(:,currst),c,'linewidth',plot_options.linewidth)
            set(gca,'linewidth',plot_options.linewidth)
            set(gca,'fontsize',plot_options.fontsize)
        end
    end

end