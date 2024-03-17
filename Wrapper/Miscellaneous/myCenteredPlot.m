function [hAx,h1,h2] = myCenteredPlot(TIME_1, SIG_1, TIME_2, SIG_2, C_1, C_2, LX, LY1, LY2, cmap)

% % -------------------------------------------------------------------------
% yyaxis left
% h1 = plot(TIME_1, SIG_1, ...
%     'Color', cmap(1,:));
% 
% xlabel(LX)
% 
% ylim([C_1, max(SIG_1) + abs(max(SIG_1) - C_1)*0.1 + eps])
% ylabel(LY1) % left y-axis
% 
% hAx(1) = gca();
% set(gca, 'YColor', cmap(1,:), 'Ydir', 'reverse');
% 
% % -------------------------------------------------------------------------
% yyaxis right
% h2 = plot(TIME_2, SIG_2, ...
%     'Color', cmap(2,:));
% 
% % ylim([C_2, max(SIG_2) + abs(max(SIG_2) - C_2)*0.1 + eps])
% ylim([C_2, max(SIG_2) + abs(max(SIG_2) - C_2)])
% ylabel(LY2) % right y-axis
% 
% hAx(2) = gca();
% set(gca, 'YColor', cmap(2,:), 'YTick', []);



% The version above using yyaxis is not yet supported by matlab2tikz ...

Q2 = C_2 - abs(C_2*0.05);

SIG_2(SIG_2 < Q2) = Q2;

if exist('cmap', 'var')
    set(gcf,'defaultAxesColorOrder', cmap);
end

[hAx,h1,h2] = plotyy(TIME_1, SIG_1, ... 
                    TIME_2, SIG_2 ); 

xlabel(LX)

ylabel(hAx(1),LY1) % left y-axis 
ylabel(hAx(2),LY2) % right y-axis 


ylim(hAx(1), [C_1, max(SIG_1) + abs(max(SIG_1) - C_1)*0.1 + eps])
ylim(hAx(2), [C_2, max(SIG_2) + abs(max(SIG_2) - C_2)]) 

hAx(1).YDir = 'reverse';

end