clear all;
close all;

alpha=1.2; %reaction rates
beta=1.2;
gamma=1.0;
size=300; %width of square simulation grid
steps=300; %number of timesteps
[a,b,c]=BZSim(alpha,beta,gamma,size,steps); %run simulation
BZTrack(alpha,beta,gamma,a,b,c,size); %export parameteric a,b,c(t) at a point
BZShow(alpha,beta,gamma,a,b,c,size,steps-1); %export a frame
BZAnimate(alpha,beta,gamma,a,b,c,size,steps); %export animation
%Note: Animation export requires third-party function "gif":
%https://www.mathworks.com/matlabcentral/fileexchange/63239-gif

function [a,b,c]=BZSim(alpha,beta,gamma,size,steps)
% Runs a forward model of the three mutually-symmetric reactions
% a_t+1 = a + a * [alpha*b - gamma*c],
% b_t+1 = b + b * [beta *c - alpha*a],
% c_t+1 = c + c * [gamma*a - beta *b],
% while simultaneously allowing diffusion over 3 grid spaces per timestep.

%concentrations
a=zeros(size,size,steps);
b=zeros(size,size,steps);
c=zeros(size,size,steps);
apad=zeros(size+2,size+2,steps);
bpad=zeros(size+2,size+2,steps);
cpad=zeros(size+2,size+2,steps);
%initial condition: random uniform concentrations in range [0,1]
a(:,:,1)=rand(size,size);
b(:,:,1)=rand(size,size);
c(:,:,1)=rand(size,size);
%diffusion kernel
kernel=ones(3,3)/9;

for t=1:steps-1
    %constrains concentrations to [0,1] so they don't blow up
    a(:,:,t)=min(max(a(:,:,t),0),1);
    b(:,:,t)=min(max(b(:,:,t),0),1);
    c(:,:,t)=min(max(c(:,:,t),0),1);
    %sets toroidal boundary condition by adding padding
    apad(:,:,t)=padarray(a(:,:,t),[1 1],'circular','both');
    bpad(:,:,t)=padarray(b(:,:,t),[1 1],'circular','both');
    cpad(:,:,t)=padarray(c(:,:,t),[1 1],'circular','both');
    %averages concentrations over neighboring cells (diffusion)
    apad(:,:,t)=conv2(apad(:,:,t),kernel,'same');
    bpad(:,:,t)=conv2(bpad(:,:,t),kernel,'same');
    cpad(:,:,t)=conv2(cpad(:,:,t),kernel,'same');
    %removes padding
    a(:,:,t)=apad(2:size+1,2:size+1,t);
    b(:,:,t)=bpad(2:size+1,2:size+1,t);
    c(:,:,t)=cpad(2:size+1,2:size+1,t);
    %calculates new concentrations
    a(:,:,t+1)=a(:,:,t)+a(:,:,t).*(alpha*b(:,:,t)-gamma*c(:,:,t));
    b(:,:,t+1)=b(:,:,t)+b(:,:,t).*(beta*c(:,:,t)-alpha*a(:,:,t));
    c(:,:,t+1)=c(:,:,t)+c(:,:,t).*(gamma*a(:,:,t)-beta*b(:,:,t));
end
end


function BZAnimate(alpha,beta,gamma,a,b,c,size,steps)
% Animates and exports the BZ reaction simulation with realistic colors.

%produce colormap
r=linspace(0.85,0.75,256);
g=linspace(0.95,0.33,256);
b=linspace(1.00,0.25,256);
cmap=[r' g' b'];
%plot first frame
f=figure();
h=surf(log(a(:,:,1))/log(10),'FaceColor','interp','EdgeColor','none');
title(['log(a) at t = 1; \alpha' num2str(alpha,'%2.1f') ' \beta' num2str(beta,'%2.1f') ' \gamma' num2str(gamma,'%2.1f')]);
caxis([-3 0]);
xlim([0 size]);
ylim([0 size]);
xticklabels('');
yticklabels('');
colormap(cmap);
colorbar;
grid off;
box on;
view(0,90);
ax=gca;
ax.DataAspectRatio=[1,1,1];
gif(['BZ-' num2str(10*alpha,'%2.0f') '-' num2str(10*beta,'%2.0f') '-' num2str(10*gamma,'%2.0f') '.gif'],'DelayTime',0.042,'LoopCount',Inf,'frame',gcf)
%plot remaining frames
for i=2:1:steps
    h=surf(log(a(:,:,i))/log(10),'FaceColor','interp','EdgeColor','none');
    title(['log(a) at t = ' num2str(i) '; \alpha' num2str(alpha,'%2.1f') ' \beta' num2str(beta,'%2.1f') ' \gamma' num2str(gamma,'%2.1f')]);
    caxis([-3 0]);
    xlim([0 size]);
    ylim([0 size]);
    xticklabels('');
    yticklabels('');
    colormap(cmap);
    colorbar;
    grid off;
    box on;
    view(0,90);
    ax=gca;
    ax.DataAspectRatio=[1,1,1];
    gif
end
end


function BZShow(alpha,beta,gamma,a,b,c,size,step)
% Plots and exports a single frame at t=step.

%produce colormap
r=linspace(0.85,0.75,256);
g=linspace(0.95,0.33,256);
b=linspace(1.00,0.25,256);
cmap=[r' g' b'];
%plot frame
f=figure();
h=surf(log(a(:,:,step))/log(10),'FaceColor','interp','EdgeColor','none');
title(['log(a) at t = ' num2str(step) '; \alpha' num2str(alpha,'%2.1f') ' \beta' num2str(beta,'%2.1f') ' \gamma' num2str(gamma,'%2.1f')]);
caxis([-3 0]);
xlim([0 size]);
ylim([0 size]);
xticklabels('');
yticklabels('');
colormap(cmap);
colorbar;
grid off;
box on;
view(0,90);
ax=gca;
ax.DataAspectRatio=[1,1,1];
saveas(f,['BZ-' num2str(10*alpha,'%2.0f') '-' num2str(10*beta,'%2.0f') '-' num2str(10*gamma,'%2.0f') '-frame.jpg']);
end


function BZTrack(alpha,beta,gamma,a,b,c,size)
% Parametrically plots and exports the path of one point of the simulation
% in a,b,c(t) phase space and exports the plot. 

%plot figure
f=figure();
f.Color='w';
hold on;
plot3(squeeze(a(size/2,size/2,:)),squeeze(b(size/2,size/2,:)),squeeze(c(size/2,size/2,:)),'-r');
title(['\alpha' num2str(alpha,'%2.1f') ' \beta' num2str(beta,'%2.1f') ' \gamma' num2str(gamma,'%2.1f')]);
text(1.1,0,0,'[a]','FontSize',12,'HorizontalAlignment','center');
text(0,1.1,0,'[b]','FontSize',12,'HorizontalAlignment','center');
text(0,0,1.1,'[c]','FontSize',12,'HorizontalAlignment','center');
ax=gca;
ax.DataAspectRatio=[1,1,1];
view(-167,10);
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);
xticks([0 0.5 1]);
yticks([0 0.5 1]);
zticks([0 0.5 1]);
ax.XRuler.FirstCrossoverValue=0;
ax.YRuler.FirstCrossoverValue=0;
ax.ZRuler.FirstCrossoverValue=0;
ax.ZRuler.SecondCrossoverValue=0;
ax.XRuler.SecondCrossoverValue=0;
ax.YRuler.SecondCrossoverValue=0;
ax.XAxis.LineWidth=1.0;
ax.YAxis.LineWidth=1.0;
ax.ZAxis.LineWidth=1.0;
ax.FontSize=12;
ax.TickLength=[0.001,0.01];
ax.TickDir='both';
ax.Box='off';
saveas(f,['BZ-' num2str(10*alpha,'%2.0f') '-' num2str(10*beta,'%2.0f') '-' num2str(10*gamma,'%2.0f') '.jpg']);
end
