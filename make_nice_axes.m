% Function to help initialise axes
function ax = make_nice_axes(ax,xlab,ylab,titletext,varargin)
    cla;
    hold(ax, 'on');
    set(ax, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
    xlabel(ax, xlab, 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(ax, ylab, 'FontSize', 14, 'Interpreter', 'latex');
    title(ax, titletext, 'FontSize', 14, 'Interpreter', 'latex');

    if nargin > 4
        xmax = varargin{1};
        xlim(ax,[0 xmax]);
    end
end