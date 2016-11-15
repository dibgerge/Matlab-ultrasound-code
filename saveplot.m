function saveplot(fname, fileformat)


if strcmp(fileformat, 'tikz')
    matlab2tikz([fname '.tikz'], 'width', '\figurewidth', 'height', '\figureheight', ...
        'externalData', false, 'floatFormat', '%.5g');
else
    export_fig(fname, ['-' fileformat], '-transparent');
end

end

