%% This script generate a file ('Table_response.res') containing the 10Be-modeled erosion rates
% as response of imposed erosion rates

%%
clear all
sample_ages = 15-[-0.06 8.75 11.74 11.76];
sample_Es = [0.0013 0.0024 0.0086 0.0035];% 10Be-measured erosion rates: samples BA-4, BA-3, BA-2, BA-1 
sample_Es_sigma = [0.0001 0.0002 0.0006 0.0003];
time_sample=-(sample_ages-15);

fid=fopen('Table_response.res','w');% name of the output file
%% five parameters are tested 
EROSION_INIT=[0.0005];% Erosion rate (cm/yr) before the step (periode 1)
EROSION2=[0.0005];% Erosion rate (cm/yr) after the step (periode 3)
EROSION1=[0.01:0.002:0.16]; %Erosion rate during the step (periode 2)
PERIODE=[100:100:8300]; % Duration of the step (of periode 2)
T0=[6.7:0.1:15]; % Time (in kyr) before present
nb_test=length(EROSION_INIT)*length(EROSION1)*length(EROSION2)*length(PERIODE)*length(T0);
%% Fixed parameters
 % density, g/cm3
param.ro = 2.7;
% attenuation length, g/cm2
param.atten = 160;
% decay constant, ln2/half-life
param.decay = log(2)/1380000;
% surface production rate, atm/(g*yr)
param.Pzero = 100;
param.mu = param.ro/param.atten;
% imposed.
%% LOOP 
%  Five loops that test all the combinations of the five degrees of freedom 
%  defined by the vector  
%  Vector time for the output
TIME_display=[-1:0.4:6.4,6.7:0.1:15,15.5,16]; % vector optimizing the vector time
inc=0;
for b1=1:length(EROSION_INIT)
   for b3=1:length(EROSION2)
      for b2=1:length(EROSION1)
         for b4=1:length(PERIODE)      
            param.erosion_init = EROSION_INIT(b1);
            param.erosion = [EROSION1(b2),EROSION2(b3)];
            % time in yr (for each erosion rate specifed in 'erosion')
            param.time = [PERIODE(b4),21000]; 
            param.ct=5;
            param.dz=1;
            [time_axis_all,surfaceEall]=Be_err_Lal(param);
            % remove first two elements from time axis all (some zeros were required)
            ind = [1 2];
            time_axis_all(ind) = [];
            time_vector=-((time_axis_all/1000)) +param.time(1)/1000;
            % T0=(T1-PERIODE(b4))/1000;
            for b5=1:length(T0) 
                datam=interp1(time_vector+T0(b5),surfaceEall,time_sample);  
                datam2=interp1(time_vector+T0(b5),surfaceEall,TIME_display);
                if (T0(b5)*1000+param.time(1))<15000 % only the models with T0<15000 are retained
                    fprintf(fid,'%2.5f %2.5f %2.5f %5.0f %5.0f %2.6f %2.6f %2.6f %2.6f ',...
                    param.erosion_init,param.erosion(1),param.erosion(2),param.time(1),T0(b5)*1000,datam(1),datam(2),datam(3),datam(4));
                    fprintf(fid,'% 3.0f ',length(TIME_display));
                    for i1=1:length(TIME_display)-1
                        fprintf(fid,'% 5.1f',TIME_display(i1));
                    end
                    i1=i1+1;
                    fprintf(fid,'% 5.1f\t',TIME_display(i1));
                    for i1=1:length(TIME_display)-1
                        fprintf(fid,'% 2.6f ',datam2(i1));
                    end
                    i1=i1+1;
                    fprintf(fid,'% 2.6f \n',datam2(i1));
                end
            end
            inc=(inc+length(T0));
            waitbar_percent=100*inc/nb_test
        end
      end
   end
end
fclose(fid)
