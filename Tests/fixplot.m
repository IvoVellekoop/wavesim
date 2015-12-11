function fixplot(xlab, ylab, label)
    set(gca(), 'FontSize', 16);
    for curve=get(gca(), 'children')
        set(curve, 'LineWidth', 1);
    end
    annotation('textbox', [0.6, 0.8, 0.3, 0.1], 'String', label, 'LineStyle', 'none',...
        'FontSize', 20, 'HorizontalAlignment', 'right');
    
    %%reposition axes to have space for xlabel
%    pos = get(gca(), 'Position');
%    pos(2) = pos(2)+0.1;
%    pos(4) = pos(4)-0.1;
%    set(gca(), 'Position', pos);
    h = xlabel(xlab, 'FontSize', 30);
%    pos = get(h, 'Position');
%    pos(2) = pos(2)+0.05;
%    set(h, 'Position', pos);
    ylabel(ylab, 'FontSize', 30);
    
end