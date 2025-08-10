function myCmap = asymColorMapWhiteZero(clims, N_colors_standard)
Cmap_standard = diverging_map((1:N_colors_standard)./N_colors_standard,...
    [0.230, 0.299, 0.754],[0.706, 0.016, 0.150]);

clims(1) = min(clims(1),0);
clims(2) = max(clims(2),0);

if clims(1) < 0 && clims(2) > 0
    posratio = clims(2)/(clims(2) - clims(1));
    
    if posratio == 0.5
        N_pos = N_colors_standard/2;
        N_neg = N_colors_standard/2;
    elseif posratio > 0.5
        N_pos = N_colors_standard/2;
        N_neg = (1-posratio)*N_colors_standard/(2*posratio);
    else
        N_pos = posratio*N_colors_standard/(2*(1-posratio));
        N_neg = N_colors_standard/2;
    end
    
    Cmap_neg = diverging_map((1:N_neg)./N_neg,...
        [0.230, 0.299, 0.754], Cmap_standard(N_colors_standard/2,:));
    Cmap_pos = diverging_map((1:N_pos)./N_pos,...
        Cmap_standard(N_colors_standard/2,:), [0.706, 0.016, 0.150]);
    myCmap = cat(1, Cmap_neg, Cmap_pos);
elseif clims(1) == 0
    N_pos = N_colors_standard/2;
    myCmap = diverging_map((1:N_pos)./N_pos,...
        Cmap_standard(N_colors_standard/2,:), [0.706, 0.016, 0.150]);
elseif clims(2) == 0
    N_neg = N_colors_standard/2;
    myCmap = diverging_map((1:N_neg)./N_neg,...
        [0.230, 0.299, 0.754], Cmap_standard(N_colors_standard/2,:));
else
    error('Colormap adjustment issue')
end
end