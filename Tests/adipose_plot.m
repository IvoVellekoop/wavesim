%%% Script to create figures of adipose dopc results

clear all, close all;
%% load dopc results (requires acces to BMPI server)
path = '\\ad.utwente.nl\TNW\BMPI\Users\Gerwin Osnabrugge\Data\';
load([path, 'dopc_results.mat']);

%% plot properties
x = mopt.pixel_size * (1:N(2));     % horizontal axis
y = mopt.pixel_size * (1:N(1));     % vertical axis
fig_size = [8 7];
fp = '../../wavesimpaper/figures/'; % filepath for figures

%% create figure of refractive index map
figure(1); clf;
imagesc(x,y,n)
h = colorbar;
fixplot('x (\mum)','y (\mum)',fig_size,'');
set(get(h,'Title'),'String',sprintf('Refractive\nindex'),'FontSize',18,'FontName','Times New Roman');
axis square;
colormap(gray)

hold on;
dim_box = [0.475,0.4,0,0];
annotation('textbox',dim_box, 'String', 'Point Source' ,'FontSize',16,'FontName','Times New Roman' , ...
    'FontWeight','Bold', 'BackgroundColor', [1 1 1], 'FaceAlpha', 0.65, 'FitBoxToText', 'on',...
    'HorizontalAlignment','center');
plot(y(end)/2,x(end*3/4),'*r','MarkerSize', 22,'LineWidth',2);
plot([x(1),x(end)],[y(10),y(10)],'r','LineWidth',3);
text(x(end)/2 - 370,-45,'Phase Conjugating Mirror', 'FontSize',18,'FontName','Times New Roman' , 'FontWeight', 'Bold');

% Draw arrow
ha = annotation('arrow',[0.465,0.505], [0.96,0.927],'LineWidth',2);

print([fp 'adipose_nmap.eps'], '-depsc2'); % print figure to eps file
%% create figure of field solutions (point source and dopc)
figure(2)
imagesc(x,y,log(abs(E1)))
h = colorbar;
fixplot('x (\mum)','y (\mum)',fig_size,'');
set(get(h,'Title'),'String','log|\psi|','FontSize',18,'FontName','Times New Roman');
axis square;
colormap(jet)
print([fp 'adipose_PS.eps'], '-depsc2'); % print figure to eps file

figure(3)
imagesc(x,y,log(abs(E2)))
h = colorbar;
fixplot('x (\mum)','y (\mum)',fig_size,'');
set(get(h,'Title'),'String','log|\psi|','FontSize',18,'FontName','Times New Roman');
axis square;
colormap(jet)
print([fp 'adipose_dopc.eps'], '-depsc2'); % print figure to eps file