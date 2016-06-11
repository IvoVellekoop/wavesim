%%% script used for plotting the results of the random medium
%%% simulations.

close all;
addpath('..');
fp = '../../wavesimpaper/figures/'; % filepath for figures
%% Load random medium simulation data
load('random_medium_results.mat');

%% Field solutions (wavesim & PSTD)
% wavesim solution
x = (-N(2)/2 + 1 : N(2)/2) /PPW;
y = (-N(1)/2 + 1 : N(1)/2) /PPW;

figure(1);
imagesc(x,y,log(abs(E_wavesim)))
h = colorbar;

fig_size = [8 7];
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
set(get(h,'Title'),'String','log|\psi|','FontSize',18,'FontName','Times New Roman');
axis square;
colormap(jet);

print([fp 'random_solution.eps'], '-depsc2'); % print figure to eps file

% PSTD solution
% figure(2);
% imagesc(x,y,log(abs(E_PSTD)))
% h = colorbar;
% set(get(h,'Title'),'String','log|\psi|','FontSize',16);
% fixplot('x (\lambda)','y (\lambda)',fig_size,'');
% colormap(jet);
% axis square;