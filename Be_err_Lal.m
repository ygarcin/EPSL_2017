%%%%%%%%%%%%%%%%%%%%%%%%%%%%Be_err_Lal.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model cosmogenic nuclide concentrations and erosion rates at the surface 
% as a function of changing erosion rates. This script allows you to run 
% any number of time steps with differing erosion rates.
% Equations from Lal (1991).  Neutron spallation only (no muogenic
% production). Taylor Schildgen 22.08.2015.
% May 2016 : implemented as a function. Julien Guillemoteau
%
% INPUT: structure param:
%			r0 : density g/cm^3
%  			atten : attenuation length, g/cm^2
%			decay : decay constant, ln2/half-life
%			Pzero : surface production rate, atm/(g*year)
%			mu :  in cm^(-1)
%			erosion_init : erosion rates in cm/yr; first initializes the model at a steady-state
%			erosion : erosion rate vector
%			time : duration vector
%			ct: factor controlling the time step (implemented as variable here but kept equal to 5)
%			dz: depth resolution (again kept equal to 1 cm)
% OUPUT: 
%			time_axis_all : vector time
%			surfaceEall : 10Be derived erosion rate
%			surfaceNall : nuclide concentration at the surface 
function [time_axis_all,surfaceEall,surfaceNall]=Be_err_Lal(param)
	decay = param.decay;
	Pzero = param.Pzero;
	mu = param.mu;
	erosion_init = param.erosion_init;
	erosion = param.erosion;
	time = param.time; 
	% depths in 1-cm increments
	% first calculate depth required; extra 15 elements added just in case.
	depth_needed = 500;
	ct=param.ct;
	for m = 1:1:length(erosion) %% change increment to 2 for speed
	depth_needed = depth_needed +(erosion(1,m)*time(1,m)*ct);
	end
	dz=param.dz;
	depths = [0:dz:depth_needed];
	% initial nuclide concentrations at depth (saturation at initial rate).
	N = (Pzero*exp(-mu.*depths/ct))/(mu*erosion_init + decay);
	% alternative: start with concentration of 0 at all depths
	% N = zeros(1,length(depths));
	% initialize variables that will change with each time step
	totaltime = 0;
	surfaceNall = [];
	surfaceEall = [];
	time_axis_all = [0,0];
	erosion_imposed_all = [];
	% This loop iterates over each erosion-rate step.
	for n = 1:1:length(erosion)
		% time step is the time required to erode 1 cm of material; can later 
		% change to mm if needed?
		%tstep = round((1/erosion(1,n))/5);
		tstep = (dz/erosion(1,n))/ct;
		t = time(1,n);
		% Surface concentration of nuclides and surface erosion rates
		Ns = round(t/tstep);
		surfaceN = zeros(1,Ns);
		surfaceE = zeros(1,Ns);
		  % For each time step, calculate change in N, surfaceN, surfaceE, and adjust
		  % depths (make erosion happens)
		  for i = 0:tstep:t
			   % calculate change in N for each depth for a given time step 
			   % (decay and new production)
			   dN = -N*decay*tstep + Pzero*exp(-mu*depths/ct)*tstep;
			   N = N + dN;
			   Nnew = ones(1,length(N)-1);
				for j = 1:length(Nnew)
					 Nnew(1,j) = N(1,j+1);% fill new N values into the Nnew array
					 depthsnew = [0:1:length(Nnew)-1];% reduce depth array by 1
				end 
			   incr = round(i/tstep);
			   surfaceN(1,incr+1) = Nnew(1,1);
			   % convert surface concentrations (surfaceN) to erosion rates
			   surfaceE(1,incr+1) = (Pzero-surfaceN(1,incr+1)*decay)/(surfaceN(1,incr+1)*mu);
			   depths = depthsnew;
			   N = Nnew;
		  end
		% for this step... 
		time_axis = [1:tstep:length(surfaceN)*tstep];
		erosion_imposed_added = ones(1,length(time_axis))*erosion(1,n);
		% add to previous steps... 
		surfaceNall = cat(2,surfaceNall,surfaceN);
		surfaceEall = cat(2,surfaceEall,surfaceE);
		erosion_imposed_all = cat(2,erosion_imposed_all,erosion_imposed_added);
		totaltime = totaltime + length(surfaceN)*tstep;
		time_axis_mod = time_axis + max(time_axis_all);
		time_axis_all = cat(2,time_axis_all,time_axis_mod);
	end
end