function fixplot()
    set(gca(), 'FontSize', 16);
    for curve=get(gca(), 'children')
        set(curve, 'LineWidth', 1);
    end
end