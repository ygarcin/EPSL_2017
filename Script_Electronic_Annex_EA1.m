%% Examples of erosion step models
clear all

%% Parameters Cosmogenic dating
% density, g/cm3
param.ro = 2.7;
% attenuation length, g/cm2
param.atten = 160;
% decay constant, ln2/half-life
param.decay = log(2)/1380000;
% surface production rate, atm/(g*yr)
param.Pzero = 100;
param.mu = param.ro/param.atten;

%% Parameters GEOLOGY
T0=16; %in kyr bp (time of the begin of the last period)
param.erosion_init = 0.0005;
param.erosion = [0.05,0.0005];
% time in yr (for each erosion rate specifed in 'erosion')
param.time = [5000,40000]; %time = [width of the step function, duration after it];
param.ct=10;
param.dz=1;

%% MODELING
[time_axis_all,surfaceEall]=Be_err_Lal(param);
% remove first two elements from time axis all (some zeros were required)
ind = [1 2];
time_axis_all(ind) = [];
time_vector=-((time_axis_all/1000))+param.time(1)/1000;

%% DISPLAY
%building the model of True errosion rate (for display)
TIME=[-2:0.1:T0+param.time(1)/1000+2];
True_err=param.erosion(2).*ones(1,length(TIME));
True_err(TIME>T0)=param.erosion(1);
True_err(TIME>T0+param.time(1)/1000)=param.erosion_init;

figure1=figure;
cindex=0.5;
axes1=axes('parent',figure1);
set(axes1,'ylim',[0 0.055],'nextplot','add')
set(axes1,'xlim',[-2 25],'nextplot','add')
p1=plot(TIME,True_err,'color',cindex.*[1 1 1],'linewidth',1,'parent',axes1);
p2=plot(time_vector+T0,surfaceEall,'color',cindex.*[1 1 1],'linestyle','--','parent',axes1);

%DATA display
sample_ages =15- [-0.06 8.75 11.74 11.76];
sample_Es = [0.0013 0.0024 0.0086 0.0035];
sample_Es_sigma = [0.0001 0.0002 0.0006 0.0003];
time_sample=15-(sample_ages);
p3=errorbar(time_sample,sample_Es,2*sample_Es_sigma,'.k','parent',axes1);
% set(gca,'XDir','reverse')
legend([p3,p2,p1],'Meas. ^{10}Be erosion rate','Mod. ^{10}Be erosion rate','True erosion rate','Location','northwest')
ylabel('Erosion rate (cm yr^{-1})')
xlabel('Time (kyr)')
title('Examples of erosion step models')
drawnow