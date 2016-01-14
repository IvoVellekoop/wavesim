%%% Script to create figures of adipose dopc results

clear all, close all;
%% load dopc results (requires acces to BMPI server)
path = '\\ad.utwente.nl\TNW\BMPI\Users\Gerwin Osnabrugge\Wavesim\';
load([path, 'dopc_results.mat']);

%% plot properties
x = mopt.pixel_size * (1:N(2));     % horizontal axis
y = mopt.pixel_size * (1:N(2));     % vertical axis
fig_size = [8 7];
fp = '../../wavesimpaper/figures/'; % filepath for figures

%% create figure of refractive index map
figure(1)
imagesc(x,y,n)
h = colorbar;
set(get(h,'Title'),'String','refractive index','FontSize',16);
fixplot('x (\mum)','y (\mum)',fig_size,'');
axis square;
colormap(gray)

print([fp 'adipose_nmap.eps'], '-depsc2'); % print figure to eps file
%% create figure of field solutions (point source and dopc)
figure(2)
imagesc(x,y,log(abs(E1)))
h = colorbar;
set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
fixplot('x (\mum)','y (\\mum)',fig_size,'');
axis square;
colormap(jet)
print([fp 'adipose_PS.eps'], '-depsc2'); % print figure to eps file

figure(3)
imagesc(x,y,log(abs(E2)))
h = colorbar;
set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
fixplot('x (\mum)','y (\\mum)',fig_size,'');
axis square;
colormap(jet)
print([fp 'adipose_dopc.eps'], '-depsc2'); % print figure to eps file