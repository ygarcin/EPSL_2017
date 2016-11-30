%% This script loads the table ('Table_response.res') computed with Build_Ensemble.m
% It compares the modeled 10Be erosion rates to the measured 10Be erosion rates
% and generates Figure 7

clear all %close all;
RESULT1=load('Table_response.res');

%% DATA display
sample_ages = 15-[-0.06 8.75 11.74 11.76];
sample_Es = [0.0013 0.0024 0.0086 0.0035];% 10Be-measured erosion rates: samples BA-4, BA-3, BA-2, BA-1 
sample_Es_sigma = [0.0001 0.0002 0.0006 0.0003];
time_sample=-(sample_ages-15); % vector time  display
figure7=figure('units','normalized','position',[0.1 0.07 0.33 0.79]);
axes1=axes('parent',figure7,'position',[0.12 0.70 0.8 0.22])
axes2=axes('parent',figure7,'position',[0.12 0.38 0.8 0.22])
axes3=axes('parent',figure7,'position',[0.12 0.07 0.8 0.22])
set(axes1,'ylim',[0 0.03])

set(get(axes1,'ylabel'),'string','Erosion rate (cm yr^{-1})')
set(get(axes1,'title'),'string','a) ^{10}Be-derived erosion rate: data and models')
set(get(axes1,'xlabel'),'string','Time (kyr)')
set(get(axes2,'ylabel'),'string','Erosion rate (cm yr^{-1})')
set(get(axes2,'xlabel'),'string','Time (kyr)')

set(axes1,'nextplot','add')
set(axes2,'nextplot','add')
set(get(axes3,'ylabel'),'string','Erosion rate (cm yr^{-1})')
set(get(axes3,'xlabel'),'string','Time (kyr)')
set(axes3,'nextplot','replace')
% set(gca,'XDir','reverse')
drawnow

%% CONSTRAINT 1 : begin of periode 2
index2=find(RESULT1(:,4)+RESULT1(:,5)<15000); % constraints for the beginning of periode 2: 
                                              % 15000 (Figure 7 A-C) or 12000 (Figure 7 D-F)
RESULT=RESULT1(index2,:);

%% CONSTRAINT 2 : data fit
DATAM=RESULT(:,6:9);
DATAD=0*DATAM;%initialisation
UNCER=0*DATAM;%initialisation
for i1=1:length(sample_Es)
    DATAD(:,i1)=sample_Es(i1).*ones(size(DATAM,1),1);
    UNCER(:,i1)=3*sample_Es_sigma(i1).*ones(size(DATAM,1),1); 
end
DIFF=DATAM-DATAD;
CRIT=100*abs(DIFF)./UNCER;
th=100;
index=find((CRIT(:,1)<th)&(CRIT(:,2)<th)); % retains only the data that fits the data points BA-3 and BA-4 within UNCER

%% DATA ANALYSIS
MODELS=RESULT(index,:);
cmap=jet;
TIME=[-5:0.1:30]; %vector time of the display
Err_ratevec=[0.002:0.002:0.16]; %vector erosion rate for the display
% matrix proba
PROB_density=zeros(length(Err_ratevec),length(TIME)); %matrix proba domain (TIME;Err_rate)
MAT_err=PROB_density; % matrix erosion rate in the domain  (TIME;Err_rate)

% Coordinates matrixes in the domain (TIME;Err_rate)
ERR_rate_mat=zeros(length(Err_ratevec),length(TIME));
TIME_mat=zeros(length(Err_ratevec),length(TIME));
for i1=1:length(TIME)
    ERR_rate_mat(:,i1)=Err_ratevec';
end
for i1=1:length(Err_ratevec)
    TIME_mat(i1,:)=TIME;
end
jet2=flipud(gray);
jet2(1,:)=[1 1 1];

dinc=1;
for i1=1:dinc:size(MODELS,1) %loop over all the models which fulfill constraints 1 and 2
    i1
    %% Imposed erosion rates
    % regenerates the associated imposed erosion rates with a fine mesh in the
    % domain (time, err_rate)
    i2=size(MODELS,1)-i1+1;
    param.erosion_init = MODELS(i2,1);
    param.erosion = [MODELS(i2,2),MODELS(i2,3)];
    % time in yr (for each erosion rate specifed in 'erosion')
    param.time = [MODELS(i2,4),21000];
    True_err=param.erosion(2).*ones(1,length(TIME));
    T0=MODELS(i2,5)/1000;
    True_err(TIME>T0)=param.erosion(1);
    True_err(TIME>T0+param.time(1)/1000)=param.erosion_init;
    
    %% Graphics A,B and D,E
    linew=0.5;    
    periode=(param.time(1));
    % define the colormap for the curves in graphics A,B and D,E
    max_periode=6500;
    min_periode=100;
    cscale=ceil(size(cmap,1)*(periode-min_periode)/(max_periode-min_periode));
    if cscale>size(cmap,1)
        cscale=size(cmap,1)
    end
    if mod(i1,2)==0 % here, only one curve per two computed is plotted (for computer memory save)
    hk=plot(TIME,True_err,'color',cmap(cscale,:),'linewidth',linew,'parent',axes2);
    set(hk,'visible','on')
    end
    NT=MODELS(i2,10);
    TIME_disp=MODELS(i2,10+1:10+NT);
    datam_disp=MODELS(i2,10+1+NT:end);
    if mod(i1,2)==0 % here, only one curve per two computed is plotted (for computer memory save) 
    hk2=plot(TIME_disp,datam_disp,'color',cmap(cscale,:),'linewidth',0.5,'parent',axes1);
    set(axes1,'ylim',[0 0.01],'xlim',[-0.5 16])
    set(axes2,'ylim',[0 0.16],'xlim',[-0.5 16])
    end
    
    %% Graphic C and F
    MAT_err=0.*MAT_err;% Building Matrix err
    for i6=1:length(TIME)
        MAT_err(Err_ratevec<=True_err(i6),i6)=1;
    end
    PROB_density=PROB_density+MAT_err; 
    if i1==size(MODELS,1) % display the final matrix only at the end of the loop 
    hk3=imagesc(TIME,Err_ratevec,PROB_density./size(MODELS,1),'parent',axes3);%proba=sum(Matt_rer{i1})/Nb_models
    set(axes3,'ylim',[0 0.16],'xlim',[-0.5 16],'ydir','normal')
    set(axes3,'clim',[0 0.4])
    shading(axes3,'flat')
    box(axes3,'on')
    box(axes2,'on')
    colormap(axes3,jet2);
    colorbar(axes3)
    drawnow
    end

end
set(axes3,'nextplot','Add')
set(get(axes3,'xlabel'),'string','Time (kyr)')
set(get(axes3,'ylabel'),'string','Erosion rate (cm yr^{-1})')
colormap(axes2,cmap)
colormap(axes1,cmap)
hb1=errorbar(time_sample,sample_Es,2*sample_Es_sigma,'.k','linewidth',1,'parent',axes1)
legend([hb1,hk],'Meas. ^{10}Be erosion rate','Mod. ^{10}Be erosion rate','Location','northwest')

%display options
ticksvalues=[min_periode,500:1000:max_periode];
%ticklabel=[];
for i1=1:length(ticksvalues)
ticklabel{1,i1}=num2str(ticksvalues(i1));
end
tickpos=(ticksvalues-min_periode)/(max_periode-min_periode);
colorbar(axes2,'Ticks',tickpos,...
         'TickLabels',ticklabel)
colorbar(axes1,'Ticks',tickpos,...
         'TickLabels',ticklabel)
     set(get(axes3,'title'),'string','c) Probability')
     set(get(axes2,'title'),'string','b) Imposed erosion rate models')
     box(axes1,'on');
