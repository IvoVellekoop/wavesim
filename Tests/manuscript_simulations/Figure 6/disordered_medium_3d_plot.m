%%% script used for plotting the results of the random medium
%%% simulations.
%Note: uncomment print commands to save figures in eps format
close all;
addpath('..');
fp = '../../figures/'; % filepath for figures
%% Load random medium simulation data
load('disordered_medium_3d.mat');

%% Field solutions (wavesim & PSTD)
% wavesim solution
fig_size = [8 7];

x = ( -N(1)/2+1:N(1)/2 )/PPW;
y = ( -N(2)/2+1:N(2)/2 )/PPW;
z = ( -N(3)/2+1:N(3)/2 )/PPW;

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

%print([fp 'disordered_3d_solution.eps'], '-depsc2'); % print figure to eps file