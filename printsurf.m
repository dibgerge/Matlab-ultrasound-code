function printsurf(hf)

if nargin == 0
    hf = gcf;
end

set(hf, 'Position', [2200, 100, 1500, 900]);

ax = get(hf, 'Children');

set(ax, 'MinorGridLineStyle', 'none', 'XGrid', 'off', ...
    'YGrid', 'off', 'ZGrid', 'off', 'Box', 'off', 'FontSize', 46, ...
    'FontWeight', 'normal', 'Position', [0.125, 0.2, 0.75, 0.8])

axis tight

end

