function printfig(hf, showMarker)

if nargin == 0
    hf = gcf;
end

if nargin < 1
    showMarker = 0;
end

set(hf, 'units','normalized','outerposition',[0 0 1 1])
% set(hf, 'Position', [100, 100, 1300, 800]);

ax = get(hf, 'Children');
pl = get(ax, 'Children');

set(ax, 'MinorGridLineStyle', 'none', 'XGrid', 'off', ...
    'YGrid', 'off', 'ZGrid', 'off', 'Box', 'off', 'FontSize', 46, ...
    'FontWeight', 'normal', 'FontName', 'Serif')

axis tight
styles = {'-', '--', '-.'};
mark  = {'o', 's', 'd', '^', '>', '<', '*', '+', 'x', 'p', 'h'};
c = colormap(lines(length(pl)));
for i=1:length(pl)
    cur_style = styles{mod(i-1,3)+1};
    cur_mark = mark{mod(i-1, length(mark))+1};
    set(pl(length(pl)-i+1), 'Color', c(i,:), 'LineStyle', cur_style)
    if showMarker
        set(pl(length(pl)-i+1), 'Marker', cur_mark, 'MarkerFaceColor', c(i,:), ...
            'MarkerEdgeColor', c(i,:), 'MarkerSize', 12, 'LineWidth', 3)
    end
end

box on
end

