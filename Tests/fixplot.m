function fixplot(xlab, ylab, size, label)
    %% Set paper size correctly
    fig = gcf();
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 size(1), size(2)];
    fig.PaperPositionMode = 'manual';
    fig.Position = 100*fig.PaperPosition+50;
    
    %% Set font size correctly
    set(gca(), 'FontSize', 16);
    
    %% Set all curves to minimum width of 1
    children = get(gca(), 'children');
    for c=1:length(children)
        curve=children(c);
        if isfield(curve,'LineWidth')
            width = get(curve, 'LineWidth');
            if (width < 1.0)
                set(curve, 'LineWidth', 1);
            end
        end
    end
    annotation('textbox', [0.6, 0.8, 0.3, 0.1], 'String', label, 'LineStyle', 'none',...
        'FontSize', 20, 'HorizontalAlignment', 'right');
    
    %%reposition axes to have space for xlabel
%    pos = get(gca(), 'Position');
%    pos(2) = pos(2)+0.1;
%    pos(4) = pos(4)-0.1;
%    set(gca(), 'Position', pos);
    xlabel(xlab, 'FontSize', 30,'FontName','Times New Roman');
%    pos = get(h, 'Position');
%    pos(2) = pos(2)+0.05;
%    set(h, 'Position', pos);
    ylabel(ylab, 'FontSize', 30,'FontName','Times New Roman');
    
end