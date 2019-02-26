clear all;
close all;

%timesteps
nsteps=100;
%size of grid box
w=200+2;
h=200+2;
%concentrations
a=zeros(h,w,nsteps);
b=zeros(h,w,nsteps);
c=zeros(h,w,nsteps);
%initial perturbation magnitude
sd=.1;
%initial condition: concentrations are 0.5, 
a(:,:,1)=0.5*(ones(h,w)+sd*randn(h,w));
b(:,:,1)=0.5*(ones(h,w)+sd*randn(h,w));
c(:,:,1)=0.5*(ones(h,w)+sd*randn(h,w));
%reaction rates
alpha=1.0;
beta=1.0;
gamma=1.0;
%range of allowable concentrations
min=0;
max=1;

for t=1:nsteps-1
    %constrains concentrations so they don't blow up
    for i=1:h 
        for j=1:w
            if a(i,j,t)<min
                a(i,j,t)=min;
            else if a(i,j,t)>max
                    a(i,j,t)=max;
                end
            end
            if b(i,j,t)<min
                b(i,j,t)=min;
            else if b(i,j,t)>max
                    b(i,j,t)=max;
                end
            end
            if c(i,j,t)<min
                c(i,j,t)=min;
            else if c(i,j,t)>max
                    c(i,j,t)=max;
                end
            end
        end
    end
    %averages concentrations over neighboring cells (diffusion)
    for i=2:h-1
        for j=2:w-1
            a(i,j,t)=mean(a(i-1:i+1,j-1:j+1,t),'all');
            b(i,j,t)=mean(b(i-1:i+1,j-1:j+1,t),'all');
            c(i,j,t)=mean(c(i-1:i+1,j-1:j+1,t),'all');
        end
    end
    %calculates new concentrations
    a(:,:,t+1)=a(:,:,t)+a(:,:,t).*(alpha*b(:,:,t)-gamma*c(:,:,t));
    b(:,:,t+1)=b(:,:,t)+b(:,:,t).*(beta*c(:,:,t)-alpha*a(:,:,t));
    c(:,:,t+1)=c(:,:,t)+c(:,:,t).*(gamma*a(:,:,t)-beta*b(:,:,t));
end

f=figure(1);
f.Position=[100 100 575 500];
h=heatmap(a(:,:,1),'GridVisible','off','ColorLimits',[0 1],'Colormap',copper);
gif('BZ.gif','DelayTime',0.042,'LoopCount',Inf,'frame',gcf)
for i=2:1:nsteps
    h=heatmap(a(:,:,i),'GridVisible','off','ColorLimits',[0 1],'Colormap',copper);
    gif
end