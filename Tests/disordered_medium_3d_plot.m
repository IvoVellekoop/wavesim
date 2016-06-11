%%% script used for plotting the results of the random medium
%%% simulations.

close all;
addpath('..');
fp = '../../wavesimpaper/figures/'; % filepath for figures
%% Load random medium simulation data
load('disordered_medium_3d.mat');

%% Field solutions (wavesim & PSTD)
% wavesim solution
fig_size = [8 7];

x = (-63:64)/PPW;
y = (-63:64)/PPW;
z = (-63:64)/PPW;

figure(1);
slice(x,y,z, log(abs(E_wavesim)), [0 8 16], 17, 17);

set(findobj(gca,'Type','Surface'),'EdgeColor','none')
h = colorbar;
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
zlabel('z (\lambda)', 'FontSize', 30,'FontName','Times New Roman');
set(get(h,'Title'),'String','log|\psi|','FontSize',18, 'FontName', 'Times New Roman');
colormap(jet);
axis square;
xlim([0 16])

print([fp 'disordered_3d_solution.eps'], '-depsc2'); % print figure to eps file