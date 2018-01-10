function Connolly15_AMPKexcito_PMCA08ATP(figurename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If you use this code, please cite:
%%% "Connolly N.M.C, D´Orsi B., Monsefi N., Huber H.J. & Prehn J.H.M. (2016) 
%%% Computational analysis of AMPK-mediated neuroprotection suggests acute 
%%% excitotoxic bioenergetics and glucose dynamics are regulated by a 
%%% minimal set of critical reactions. PLoS ONE 11(2):e0148326. 
%%% doi: 10.1371/journal.pone.0148326 Date: Feb 3 2016 (PMID: 26840769)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call eg.: Connolly15_AMPKexcito_PMCA08ATP('Fig1') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fig1: Fits parameters to single-cell data (.txt files) = Fig 1BCDE
%%%     Randomly varies independent parameters within defined range, runs
%%%     curve fitting function. Must have .txt files in same folder as .m
%%%     function
%%% Fig2: Cell-to-cell heterogeneity = Fig. 2ABCDE
%%%     Randomly vary parameters
%%%     Generate scatter and boxplots comparing model variation
%%%     single-cell experimental variation  
%%%     NB. Random variation of parameters - data points may differ
%%%     from figures!
%%% Fig3: Sensitivity Analysis = Fig 3ABC
%%%     Vary parameters one by one +/- X% 
%%%     Output min/max values and recovery duration of Ca, ATP, AMPK, Gluc  
%%% Fig 4AB: Necrosis = Fig 4AB
%%%     Vary calcium duration and plot ATP profiles. Highlights necrotic 
%%%     ATP collapse 
%%% Fig 4C: ATP output on 60 minute calcium influx = Fig 4C
%%% Fig4DEF: Necrosis = Fig4CD
%%%     Vary calcium magnitude and plot calcium profiles. Highlights 
%%%     occurrence of necrotic calcium profile and ATP collapse
%%% Fig 4H: Population ~ Fig 4H
%%%     Simulates population, determines % of viable & necrotic cells
%%%     following mild and severe excitotoxic stress. Outputs numbers (no
%%%     figure) NB. Random variation of parameters - output may differ
%%%     from figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global tt t_final options stepsize ca_mag ca_duration ca_init CaFold 

t_final = 6000;             % Duration of simulation (sec)
stepsize = 1;               % Simulation Step Size [s]
tt = 0:stepsize:t_final;    % Allocate same # of data points to each result
% ODE solver options
options = odeset('AbsTol',1.00E-019,'RelTol',1E-8,'MaxStep',50); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcium influx parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca_mag = 34.5;
ca_init = 540;             % Time of Calcium influx initiation (s)
ca_duration = 600;         % Duration (s) of calcium influx
% CaC Fold Change (InitConc -> plateau during influx) 
CaFold = 3.7;               
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define initial conditions and kinetic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[InitConc kon koff] = defineParam;      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ode solver in isolation and plot output (if desired)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,kon,koff,...
%    ca_mag,ca_duration,ca_init); 
%figure, plot(t/60,y(:,9)/y(1,9))
%xlabel('Time (min)'),ylabel('ATP concentration (fold change over baseline)')
%ylim([0.5 1.2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(figurename,'Fig1')
    % Number of parameter variations (set to 2 to plot figures with
    % optimised parameters)
    setCond = 2;    
    curveFitting(InitConc, kon, koff, setCond)
elseif strcmp(figurename,'Fig2')
    setCond = 100;  % # simulations to be run
    c2cheterogeneity(InitConc, kon, koff, setCond)
elseif strcmp(figurename,'Fig3')
    setCond = [1 2 3]; % 1 = vary initConc, 2 = vary k_on, 3 = vary k_off
    Sensitivity_Analysis(InitConc, kon, koff, setCond)
elseif strcmp(figurename,'Fig4AB')
    setCond = 2; % 1 = vary ca_mag (Fig 4AB), 2 = vary ca_duration (Fig 4CD)
    necrosis(InitConc, kon, koff, setCond)
elseif strcmp(figurename,'Fig4C')
    ca_duration = 3600;     % 60 min calcium influx
    % Call ODE solver
    [t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,kon,koff,...
        ca_mag,ca_duration,ca_init); 
    % Plot ATP output
    figure, plot(t/60,y(:,9)/y(1,9))
    xlabel('Time (min)')
    ylabel('ATP concentration (fold change over baseline)')
    ylim([0.5 1.2])
elseif strcmp(figurename,'Fig4DEF')
    setCond = 1; % 1 = vary ca_mag (Fig 4AB), 2 = vary ca_duration (Fig 4CD)
    necrosis(InitConc, kon, koff, setCond)
elseif strcmp(figurename, 'Fig4H')
    setCond = 1000;   % number of simulations ('cells') in each experiment
    temp_print = 0; %Set to 1 if you want to print indivicual cell fate
    runPopulation(InitConc, kon, koff, setCond, temp_print)
else
    disp('Incorrect function call. Try eg.: Connolly15_AMPKexcito_PMCA08ATP(''Fig1'')')
end

end

function curveFitting(InitConc, k_on, k_off, setCond)

global tt t_final options ca_mag ca_duration ca_init CaFold
minerror = 0;                 % Initialising values
minerror_individ = 0;
sim_number = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain data from single-cell experimental data
% .txt files must be saved in same folder as .m file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileid = fopen('Survival_data.txt');
C = textscan(fileid,'%f %f %f %f %f');
fclose(fileid);
            
data_time = C{1};
data = zeros(size(data_time,1),size(InitConc,2));  % Matrix to store exptal data 
data(:,1) = C{2};                     %Ca_c
data(:,11) = C{3};                    %pAMPKFret
data(:,8) = C{4};                     %Glucose
data(:,9) = C{5};                     %ATP
            
[data_Glucmin data_ATPmin data_AMPKmin data_CaCmin data_Glucduration...
    data_ATPduration data_AMPKduration data_CaCduration data_AMPKplot... 
    data_Glucplot data_ATPplot data_Fluo4plot] = getExptData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulations multiple times, varying parameters, and fitting curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for x = 1:setCond        
         
   fprintf('Simulation %i: Solving...',sim_number)
               
   % Calling ode solver (first run is with original parameters)
   [t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,k_on,k_off,ca_mag,ca_duration,ca_init); 
        
   % Convert all model-generated curves to set number of data points (set by tt)
   for i=1:size(InitConc,2)
       values(:,i)=spline(t,y(:,i),tt);
   end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%Calculate fold change for each model species
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1:size(InitConc,2)
     for m = 1:size(tt,2)
        % Baseline defined as mean of first 5 values
        foldchange(m,i) = (values(m,i)/values(mean(1:5),i));    
     end
   end
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%Fit curves
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [species_error,totalerror] = curvefit(foldchange,data,data_time,InitConc);
                   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%Calculate error
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   errortype = 1;      %Set to 1 if minimising TOTAL error, set to 2 if minimising one individual error (then set which individual error that is!)
   individerror = 8;    %1 for Ca, 8 for Gluc, 9 for ATP, 11 for pAMPKF
                
   switch errortype
       case 1
       %If TOTAL error is smallest yet, store values
          if minerror==0 || totalerror < minerror 
             minerror = totalerror;
             ca_mag_final = ca_mag;
             foldchange_final = foldchange;
          end
          fprintf(' Total Error=%1.2i, min=%1.2i\n',totalerror,minerror);
                        
          case 2
             %If INDIVIDUAL error is smallest yet, store values
             if minerror_individ==0 || species_error(individerror) < minerror_individ
                 minerror_individ = species_error(individerror);
                 minerror = totalerror;
                 ca_mag_final = ca_mag;
                 foldchange_final = foldchange;
             end
             fprintf(' Individual Error(%i) = %1.2i, min = %1.2i\n',individerror,species_error(individerror),minerror_individ);
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%Vary parameters
   %%%Can select which parameters to vary, and in which range
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Allowing simulation to run multiple times, using sim_number as index    
   sim_number = sim_number+1;       
   neg=1;
       
   while neg==1 
        neg=0;
       
       %k_off(2) = 23e-6 + (116e-6-23e-6)*rand(1);
       %k_off(3) = 11.6e-6 + (150e-6 -11.6e-6)*rand(1);
       %k_on(3) = 20e-3 + (140e-3 - 20e-3)*rand(1);      
       %k_on(10) = 0.1e-9+(10e-6-0.1e-9)*rand(1);
       %k_on(12) = 1e-6+(1e+3-1e-6)*rand(1);
       %k_off(11) = 0.1e-9+(7e-6-0.1e-9)*rand(1);
       %k_on(13) = 578e-6 + (69e-3-578e-6)*rand(1);
       %k_on(14) = 1 + (0.1 - 1)*rand(1);
       %k_on(15) = 1e-10 + (1e-5 - 1e-10)*rand(1);
                    
        % Vary concentrations between 1nM - 2uM (signalling proteins)
       %InitConc(1) = 50+(400-50)*rand(1);                  %Ca_c
       %InitConc(2) = 100+(500-100)*rand(1);                %Ca_m
       %InitConc(3) = 10000+(60000-10000)*rand(1);          %AMP 
       %InitConc(4) = 1+(2000-1)*rand(1);                   %AMPK
       %InitConc(5) = 1+(2000-1)*rand(1);                   %pAMPK
       %InitConc(6) = 1+(2000-1)*rand(1);                   %GLUT3
       %InitConc(7) = 1+(2000-1)*rand(1);                   %GLUT3m 
       %InitConc(8) = 3000000+(11000000-3000000)*rand(1);   %Glucose
       %InitConc(9) = 2000000+(3000000-2000000)*rand(1);    %ATP
       %InitConc(10) = 100000+(700000-100000)*rand(1);      %ADP
       %InitConc(11) = InitConc(5);     % AMPKAR should always = pAMPK
     
        %--------Update dependent kinetic parameters----------
        k_on = updateSS(InitConc,k_on,k_off);
        
        % Ensure Ca curve (input) fits Fluo4 regardless of parameter changes
        ca_mag = (k_on(15)*(InitConc(1)*InitConc(9)*CaFold*0.75))-(k_on(15)*InitConc(1)*InitConc(9));      
                               
        %Ensure no negative k values get through
        for m = 1:size(k_on,2)          
           if k_on(m)<0, neg = 1; end
           if k_off(m)<0, neg = 1; end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Plot output
%%%Model output overlaid on median and IQ regions of experimental data for
%%%CaC, ATP, pAMPK and Glucose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If initial values are 0, foldchanges cannot be calculated
for i=1:size(InitConc,2)
  if values(1,i)==0
     fprintf('\nError: Value(1,%i)==0, fold change caclulations are meaningless!!\n',i);
  end
end

figure
subplot(2,2,1), hold on
ylabel('CaC FoldChange'), xlabel('Simulation Time [s]'),
% Experimental data - IQ regions
area(data_Fluo4plot(:,1)/60,data_Fluo4plot(:,4),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[1 1 1]);
area(data_Fluo4plot(:,1)/60,data_Fluo4plot(:,3),'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
% Experimental data - median
plot(data_Fluo4plot(:,1)/60,data_Fluo4plot(:,2),'k:')
% Model Output
plot(tt/60,foldchange_final(:,1),'k','LineWidth',2)
xlim([0 t_final/60])
ylim([0 4.5])

subplot(2,2,3), hold on
ylabel('pAMPK FoldChange'), xlabel('Simulation Time [s]')
area(data_AMPKplot(:,1)/60,data_AMPKplot(:,4),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[1 1 1]);
area(data_AMPKplot(:,1)/60,data_AMPKplot(:,3),'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
plot(data_AMPKplot(:,1)/60,data_AMPKplot(:,2),'k:')
plot(tt/60,foldchange_final(:,11),'k','LineWidth',2)
xlim([0 t_final/60])
ylim([0 1.6])

subplot(2,2,4), hold on
ylabel('Glucose FoldChange'), xlabel('Simulation Time [s]')
area(data_Glucplot(:,1)/60,data_Glucplot(:,4),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[1 1 1])
area(data_Glucplot(:,1)/60,data_Glucplot(:,3),'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
plot(data_Glucplot(:,1)/60,data_Glucplot(:,2),'k:')
plot(tt/60,foldchange_final(:,8),'k','LineWidth',2)
xlim([0 t_final/60])
ylim([0 1.8])

subplot(2,2,2), hold on
ylabel('ATP FoldChange'), xlabel('Simulation Time [s]'),
area(data_ATPplot(:,1)/60,data_ATPplot(:,4),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[1 1 1])
area(data_ATPplot(:,1)/60,data_ATPplot(:,3),'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
plot(data_ATPplot(:,1)/60,data_ATPplot(:,2),'k:')
plot(tt/60,foldchange_final(:,9),'k','LineWidth',2)
xlim([0 t_final/60])
ylim([0 1.2])
fprintf('Ca-mag (optimal fit): %i\n',ca_mag_final)
 
end

function c2cheterogeneity(InitConc, k_on, k_off, setCond)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Independent parameters are varied randomly within a set range
%%% Generates cell-to-cell heterogeneity figures
%%% Plots mean and IQs of model variations over FWHM and min/max values 
%%% of single-cell traces used in Fig 2A
%%% Generates scatter and boxplots comparing model variation with single-cell
%%% variation = Fig. 2BCDE and Fig 3E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global tt t_final options ca_mag ca_duration ca_init

% Store initial values as variables will be changed within function
InitConc_orig = InitConc;          
k_on_orig = k_on;
k_off_orig = k_off;
ca_mag_orig = ca_mag;
ca_duration_orig = ca_duration;
ca_init_orig = ca_init;

% Check for negative k-values
neg = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain data from single-cell traces, min/max,durations, medians and
% IQ regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data_Glucmin data_ATPmin data_AMPKmin data_CaCmin data_Glucduration...
    data_ATPduration data_AMPKduration data_CaCduration data_AMPKplot... 
    data_Glucplot data_ATPplot data_Fluo4plot] = getExptData;

%pre-allocate matrices for speed
minGluc = zeros(setCond,1);
minATP = zeros(setCond,1);
minAMPKF = zeros(setCond,1);
minCaC = zeros(setCond,1);
durationGluc = zeros(setCond,1);
durationATP = zeros(setCond,1);
durationAMPKF = zeros(setCond,1);
durationCaC = zeros(setCond,1);

for i = 1:setCond                 % Number of times simulation will be run
        
    InitConcAll(i,:) = InitConc;  % To test distribution of varied parameters
    k_onAll(i,:) = k_on;
    k_offAll(i,:) = k_off;
        
    if neg == 1
        fprintf('Cell %i: Parameter set returns negative kinetic constants.\n',i);    % Don't perform calculations with negative parameter values!
        neg = 0;
         
    else
        fprintf('Cell %i: Solving...  ',i)
              
        % Call ode solver (function runs once with excel values)
        [t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,k_on,k_off,ca_mag,ca_duration,ca_init);  
                
            
        % Convert all curves to set number of data points (set by tt)
        for j = 1:size(InitConc,2)
            values(:,j) = spline(t,y(:,j),tt);
        end

        % Calculate fold change for each species (model output)
        for j = 1:size(InitConc,2)
           for m = 1:size(tt,2)
               foldchange(m,j) = (values(m,j)/values(1,j));      %Baseline is first value
           end
        end
             
        % Calculate min/max and duration of ATP, AMPKFret and Glucose
        % (and absolute Ca value) curves (all called min*) from model
        % predictions
        [minCaC(i,1), index_minCaC] = max(foldchange(:,1));
        [minGluc(i,1), index_minGluc] = max(foldchange(:,8));
        [minATP(i,1), index_minATP] = min(foldchange(:,9));
        [minAMPKF(i,1), index_minAMPKF] = max(foldchange(:,11));
            
        % Find index of first value before min/max > 0.98 (response
        % onset) and first value >1.02 after (response end)
        index_CaConset = find(foldchange(1:index_minCaC,1) < 1.02,1,'last');
        index_Gluconset = find(foldchange(1:index_minGluc,8) < 1.02,1,'last');
        index_ATPonset = find(foldchange(1:index_minATP,9) > 0.98,1,'last');
        index_AMPKFonset = find(foldchange(1:index_minAMPKF,11) < 1.02,1,'last');
            
        index_CaCend = find(foldchange(index_minCaC:end,1) < 1.02,1,'first')...
             + index_minCaC - 1;         % Need to get the ACTUAL index of end
        index_Glucend = find(foldchange(index_minGluc:end,8) < 1.02,1,'first')...
             + index_minGluc - 1;
        index_ATPend = find(foldchange(index_minATP:end,9) > 0.98,1,'first')...
             + index_minATP - 1;
        index_AMPKFend = find(foldchange(index_minAMPKF:end,11) < 1.02,1,'first')...
             + index_minAMPKF - 1;
                                 
        % If responses do not return to baseline, set values to NaN so
        % they won't be plotted 
        if isempty(index_CaConset) || isempty(index_Gluconset) || ...
           isempty(index_AMPKFonset) || isempty(index_ATPonset) ||...
           isempty(index_Glucend) || isempty(index_ATPend) ||...
           isempty(index_AMPKFend) || isempty(index_CaCend) || minATP(i,1)<0.35
           durationCaC(i,1) = NaN;
           durationGluc(i,1) = NaN;
           durationATP(i,1) = NaN;
           durationAMPKF(i,1) = NaN;
           minCaC(i,1) = NaN;
           minATP(i,1) = NaN;
           minGluc(i,1) = NaN;
           minAMPKF(i,1) = NaN;
           fprintf('Simulation #%i: Not all species returned to baseline (cell considered necrotic).\n',i);
       else
           % Recovery duration is TIME (tt) of end - TIME (tt) of ca_init_orig - ca_duration. 
           % ca_duration variability matches variability in times of MK801 addition in expts
           % Using ca_init_orig (same for all simulations) to
           % calclulate recovery time models variability in the
           % time of drug diffusion (i.e. drug addition time is the same, but onset of injury might not have been)                
           % This allows comparison with the recovery of single-cells as there is
           % variation in treatment time for experiments
           % This also models variations in drug difussion and glut
           % receptors etc
           fprintf('\n');
           durationCaC(i,1) = tt(index_CaCend) - tt(round(ca_init_orig)) - ca_duration;
           durationGluc(i,1) = tt(index_Glucend) - tt(round(ca_init_orig)) - ca_duration;
           durationATP(i,1) = tt(index_ATPend) - tt(round(ca_init_orig)) - ca_duration;
           durationAMPKF(i,1) = tt(index_AMPKFend) - tt(round(ca_init_orig)) - ca_duration;
                
           %Convert into minutes to align with trace data
           durationCaC(i,1) = durationCaC(i,1)/60;
           durationGluc(i,1) = durationGluc(i,1)/60;
           durationATP(i,1) = durationATP(i,1)/60;
           durationAMPKF(i,1) = durationAMPKF(i,1)/60;
       end
         
      %Save values to calculate median of all results
      foldchange_Ca(:,i) = foldchange(:,1);
      foldchange_ATP(:,i) = foldchange(:,9);    %ATP
      foldchange_AMPKF(:,i) = foldchange(:,11);   %pAMPKFret
      foldchange_Gluc(:,i) = foldchange(:,8);    %Gluc
            
      foldchange_Cam(:,i) = foldchange(:,2);
      foldchange_AMP(:,i) = foldchange(:,3);    %AMP
      foldchange_AMPK(:,i) = foldchange(:,4);   
      foldchange_pAMPK(:,i) = foldchange(:,5);
      foldchange_GLUT3(:,i) = foldchange(:,6);
      foldchange_GLUT3m(:,i) = foldchange(:,7);
      foldchange_ADP(:,i) = foldchange(:,10);
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Vary initial concentrations and kinetics
      % Parameter variations fluctuate around original concentrations
      % and kinetics. Calcium magnitude, duration and time of onset
      % are also varied. Dependent kinetic parameters are updated
      % according to steady-state constraints
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ca_vary = 0.1;    % % of calcium variations (see population fn)
      [InitConc,k_on,k_off,ca_mag,ca_duration,ca_init]...
          = population(InitConc_orig,k_on_orig,k_off_orig,ca_mag_orig,ca_duration_orig,ca_init_orig,i,ca_vary); 
                
    end
end                                     %End of parameter variations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Create Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign x-axis values for plotting scatter plot & boxplots
boxC = 1;
boxA = 3;
boxK = 5;
boxG = 7;

for x = 1:size(minATP,1)
   minATP(x,2) = boxA;
   minATP(x,3) = boxA + 0.3 + rand/10;
   minGluc(x,2) = boxG;
   minGluc(x,3) = boxG + 0.3 + rand/10;
   minAMPKF(x,2) = boxK;
   minAMPKF(x,3) = boxK + 0.3 + rand/10; 
   minCaC(x,2) = boxC;
   minCaC(x,3) = boxC + 0.3 + rand/10;
            
   durationATP(x,2) = boxA;
   durationATP(x,3) = boxA + 0.3 + rand/10;
   durationGluc(x,2) = boxG;
   durationGluc(x,3) = boxG + 0.3 + rand/10;
   durationAMPKF(x,2) = boxK;
   durationAMPKF(x,3) = boxK + 0.3 + rand/10;
   durationCaC(x,2) = boxC;
   durationCaC(x,3) = boxC + 0.3 + rand/10;
end
    
% In exptal data, replace min/max values = 0 with NaN so they are not plotted
for x = 1:size(data_Glucmin,1)
   if data_Glucmin(x,1) == 0
       data_Glucmin(x,:) = NaN;
       data_Glucduration(x,:) = NaN;
   else
       data_Glucmin(x,2) = 8;
       data_Glucmin(x,3) = 8 + 0.3 + rand/10;
       data_Glucduration(x,2) = 8;
       data_Glucduration(x,3) = 8 + 0.3 + rand/10;
   end
end
for x = 1:size(data_ATPmin,1)
   if data_ATPmin(x,1) == 0
        data_ATPmin(x,:) = NaN;
        data_ATPduration(x,:) = NaN;
   else
        data_ATPmin(x,2) = 4;
        data_ATPmin(x,3) = 4 + 0.3 + rand/10;
        data_ATPduration(x,2) = 4;
        data_ATPduration(x,3) = 4 + 0.3 + rand/10;
   end
end
for x = 1:size(data_AMPKmin,1)
   if data_AMPKmin(x,1) == 0
        data_AMPKmin(x,:) = NaN;
        data_AMPKduration(x,:) = NaN;
   else
        data_AMPKmin(x,2) = 6;
        data_AMPKmin(x,3) = 6 + 0.3 + rand/10;
        data_AMPKduration(x,2) = 6;
        data_AMPKduration(x,3) = 6 + 0.3 + rand/10;
   end
end
for x = 1:size(data_CaCmin,1)
   if data_CaCmin(x,1) == 0
        data_CaCmin(x,:) = NaN;
        data_CaCduration(x,:) = NaN;
   else
        data_CaCmin(x,2) = 2;
        data_CaCmin(x,3) = 2 + 0.3 + rand/10;
        data_CaCduration(x,2) = 2;
        data_CaCduration(x,3) = 2 + 0.3 + rand/10;
   end
end
        
% Combine model and traces to be plotted on same figure
minCaC_boxplot = [minCaC; data_CaCmin];
minATP_boxplot = [minATP; data_ATPmin];
minGluc_boxplot = [minGluc; data_Glucmin];
minAMPK_boxplot = [minAMPKF; data_AMPKmin];      
        
durationCaC_boxplot = [durationCaC; data_CaCduration];
durationATP_boxplot = [durationATP; data_ATPduration];       
durationGluc_boxplot = [durationGluc; data_Glucduration];
durationAMPK_boxplot = [durationAMPKF; data_AMPKduration];
        
min_boxplot = [minCaC_boxplot; minATP_boxplot; minAMPK_boxplot; minGluc_boxplot];
duration_boxplot = [durationCaC_boxplot; durationATP_boxplot; durationGluc_boxplot; durationAMPK_boxplot];
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plots of durations and min/max values for each simulation
% Overlay boxplots of values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
boxplot(min_boxplot(:,1),min_boxplot(:,2))
title('Ca_c (model/expt), ATP (model/expt), AMPK (model/expt), Glucose (model/expt)')
hold all, scatter(min_boxplot(:,3),min_boxplot(:,1),'.','k')
ylabel('Min. / Max. Values')
% Colour first, second, second last and last parameter sets
tempC = size(minCaC_boxplot,1);
tempA = tempC + size(minATP_boxplot,1);
tempK = tempA + size(minAMPK_boxplot,1);
tempx = [1 tempC+1 tempA+1 tempK+1];
tempx2 = [setCond tempC+setCond tempA+setCond tempK+setCond];
scatter(min_boxplot(tempx,3),min_boxplot(tempx,1),'.','g')  % Original parameter set
scatter(min_boxplot(tempx+1,3),min_boxplot(tempx+1,1),'.','y')
scatter(min_boxplot(tempx2,3),min_boxplot(tempx2,1),'.','r')
scatter(min_boxplot(tempx2-1,3),min_boxplot(tempx2-1,1),'.','c')
        
figure 
boxplot(duration_boxplot(:,1),duration_boxplot(:,2))
title('Ca_c (model/expt), ATP (model/expt), AMPK (model/expt), Glucose (model/expt)')
hold all, scatter(duration_boxplot(:,3),duration_boxplot(:,1),'.','k')
ylabel('Recovery Duration (min)')
scatter(duration_boxplot(tempx,3),duration_boxplot(tempx,1),'.','g')
scatter(duration_boxplot(tempx+1,3),duration_boxplot(tempx+1,1),'.','y')
scatter(duration_boxplot(tempx2,3),duration_boxplot(tempx2,1),'.','r')
scatter(duration_boxplot(tempx2-1,3),duration_boxplot(tempx2-1,1),'.','c')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate medians and IQs of model simulations
median_CaC = median(foldchange_Ca,2);
median_ATP = median(foldchange_ATP,2);
median_AMPKF = median(foldchange_AMPKF,2);
median_Gluc = median(foldchange_Gluc,2);
        
quart1_CaC = prctile(foldchange_Ca,25,2);
quart1_ATP = prctile(foldchange_ATP,25,2);
quart1_AMPKF = prctile(foldchange_AMPKF,25,2);
quart1_Gluc = prctile(foldchange_Gluc,25,2);
        
quart3_CaC = prctile(foldchange_Ca,75,2);
quart3_ATP = prctile(foldchange_ATP,75,2);
quart3_AMPKF = prctile(foldchange_AMPKF,75,2);
quart3_Gluc = prctile(foldchange_Gluc,75,2);
        
% Plot mean and IQs of predictions on top of mean and IQs of
% experiments
figure, subplot(2,2,1), hold on
plot(data_Fluo4plot(:,1),data_Fluo4plot(:,4),'b:');
plot(data_Fluo4plot(:,1),data_Fluo4plot(:,3),'b:');
plot(data_Fluo4plot(:,1),data_Fluo4plot(:,2),'k')
plot(tt,quart3_CaC,'r:'); 
plot(tt,quart1_CaC,'r:');
plot(tt,median_CaC,'k');
xlim([0 4200]);
ylim([0 5]);

subplot(2,2,2), hold on
plot(data_ATPplot(:,1),data_ATPplot(:,4),'b:')
plot(data_ATPplot(:,1),data_ATPplot(:,3),'b:')
plot(data_ATPplot(:,1),data_ATPplot(:,2),'k')
plot(tt,quart3_ATP,'r:'); 
plot(tt,quart1_ATP,'r:');
plot(tt,median_ATP,'k');
xlim([0 4200]);
ylim([0 1.2]);
        
subplot(2,2,3), hold on
plot(data_AMPKplot(:,1),data_AMPKplot(:,4),'b:');
plot(data_AMPKplot(:,1),data_AMPKplot(:,3),'b:');
plot(data_AMPKplot(:,1),data_AMPKplot(:,2),'k')
plot(tt,quart3_AMPKF,'r:'); 
plot(tt,quart1_AMPKF,'r:');
plot(tt,median_AMPKF,'k');
xlim([0 4200]);
ylim([0 1.6]);

subplot(2,2,4), hold on
plot(data_Glucplot(:,1),data_Glucplot(:,4),'b:')
plot(data_Glucplot(:,1),data_Glucplot(:,3),'b:')
plot(data_Glucplot(:,1),data_Glucplot(:,2),'k')
plot(tt,quart3_Gluc,'r:'); 
plot(tt,quart1_Gluc,'r:');
plot(tt,median_Gluc,'k');
xlim([0 4200]);
ylim([0 1.8]);
        
end

function Sensitivity_Analysis(InitConc, k_on, k_off, setCond)

global tt t_final options ca_mag ca_duration ca_init

paramChangeValue = [0.5 0.75 1 1.5 2];  % Parameter changes (multiples)
InitConc_orig = InitConc;          % Storing original values (values will change)
k_on_orig = k_on;
k_off_orig = k_off;

neg = 0;                     % Check for negative k-values
           
%pre-allocate matrices for speed
temp = zeros(size(paramChangeValue,2),1);
minGluc = temp;
minATP = temp;
minAMPKF = temp;
durationGluc = temp;
durationATP = temp;
durationAMPKF = temp;

for a = setCond    
   switch a
     case 1
        paramChangeID = 1:size(InitConc,2);
        % Create figure for plotting within for loops
        h(1) = figure('Name','ATP Minimum Value (Vary InitConc)');                  
        h(2) = figure('Name','pAMPK Maximum Value (Vary InitConc)');
        h(3) = figure('Name','Glucose Maximum Value (Vary InitConc)');
        h(4) = figure('Name','ATP Recovery (Vary InitConc)');
        h(5) = figure('Name','pAMPK Recovery (Vary InitConc)');
        h(6) = figure('Name','Glucose Recovery (Vary InitConc)');
        nameGraph = {'Ca_c'; 'Ca_m'; 'AMP'; 'AMPK'; 'pAMPK'; 'GLUT3';... 
                    'GLUT3m'; 'Glucose'; 'ATP'; 'ADP'; 'AMPKAR'};
                
     case 2
        paramChangeID = [3 10 12 13 15]; 
        h(1) = figure('Name','ATP Minimum Value (Vary k-on)');                  
        h(2) = figure('Name','pAMPK Maximum Value (Vary k-on)');
        h(3) = figure('Name','Glucose Maximum Value (Vary k-on)');
        h(4) = figure('Name','ATP Recovery (Vary k-on)');
        h(5) = figure('Name','pAMPK Recovery (Vary k-on)');
        h(6) = figure('Name','Glucose Recovery (Vary k-on)');
        nameGraph = {'k-on1'; 'k-on2'; 'k-on3'; 'k-on4'; 'k-on5';... 
              'k-on6'; 'k-on7'; 'k-on8'; 'k-on9'; 'k-on10'; 'k-on11';...
              'k-on12'; 'k-on13'; 'k-on14'; 'k-on15'; 'k-on16'};
                
     case 3
        paramChangeID = [2 3 11];
        h(1) = figure('Name','ATP Minimum Value (Vary k-off)');                  
        h(2) = figure('Name','pAMPK Maximum Value (Vary k-off)');
        h(3) = figure('Name','Glucose Maximum Value (Vary k-off)');
        h(4) = figure('Name','ATP Recovery (Vary k-off)');
        h(5) = figure('Name','pAMPK Recovery (Vary k-off)');
        h(6) = figure('Name','Glucose Recovery (Vary k-off)');
        nameGraph = {'k-off1'; 'k-off2'; 'k-off3'; 'k-off4'; 'k-off5';... 
             'k-off6'; 'k-off7'; 'k-off8'; 'k-off9'; 'k-off10';...
             'k-off11'; 'k-off12'; 'k-off13'; 'k-off14'; 'k-off15'};
   end
        
   for x = paramChangeID       
     for i = 1:size(paramChangeValue,2)  % For each changed parameter value...   
      
        % Ensure you're only changing one parameter value
        InitConc = InitConc_orig;   
        k_on = k_on_orig;
        k_off = k_off_orig;
        
        %%%%%%%% Vary initial concentrations and kinetics
        switch a
            case 1
                InitConc(x) = InitConc_orig(x) * paramChangeValue(i);
                InitConc(11) = InitConc(5); %AMPKAR should always = pAMPK
                
            case 2, k_on(x) = k_on_orig(x) * paramChangeValue(i);
            case 3, k_off(x) = k_off_orig(x) * paramChangeValue(i);
        end
        
        %--------Update dependent kinetic parameters----------
        k_on = updateSS(InitConc,k_on,k_off);
                          
        %Ensuring no negative k values get through
        for m = 1:size(k_on,2)          
           if k_on(m)<0, neg=1; end
           if k_off(m)<0, neg=1; end
        end
        
        if neg == 1
            % Don't perform calculations with negative parameter values!
            fprintf('Simulation %i: Parameter set returns negative kinetic constants.\n',i);    
            neg = 0;
         
        else
            switch a
             case 1, fprintf('Simulation %i.%i: InitConc(%i) = %i  ',x,i,x,InitConc(x))
             case 2, fprintf('Simulation %i.%i: k-on(%i) = %i  ',x,i,x,k_on(x))
             case 3, fprintf('Simulation %i.%i: k-off(%i) = %i  ',x,i,x,k_off(x))
            end
                          
            % Call ode solver 
            
            [t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,k_on,k_off,ca_mag,ca_duration,ca_init);  
                        
            % Convert all curves to set number of data points (set by tt)
            for j = 1:size(InitConc,2)
                values(:,j) = spline(t,y(:,j),tt);
            end

            % Calculate fold change for each species (i.e. normalise model output)
            for j = 1:size(InitConc,2)
                for m = 1:size(tt,2)
                    foldchange(m,j) = (values(m,j)/values(1,j));      %Baseline is first value
                end
            end
             
            % Calculate min/max and duration of ATP, AMPKFret and Glucose
            % (and absolute Ca value) curves (all called min*)
            [minCaC(i,1), index_minCaC] = max(foldchange(:,1));
            [minGluc(i,1), index_minGluc] = max(foldchange(:,8));
            [minATP(i,1), index_minATP] = min(foldchange(:,9));
            [minAMPKF(i,1), index_minAMPKF] = max(foldchange(:,11));
                        
            % Find index of first value before min/max > 0.95 (response
            % onset) and first value after (response end)
            index_Gluconset = find(foldchange(1:index_minGluc,8) < 1.02,1,'last');
            index_ATPonset = find(foldchange(1:index_minATP,9) > 0.98,1,'last');
            index_AMPKFonset = find(foldchange(1:index_minAMPKF,11) < 1.02,1,'last');
            
            index_Glucend = find(foldchange(index_minGluc:end,8) < 1.02,1,'first');
            index_ATPend = find(foldchange(index_minATP:end,9) > 0.98,1,'first');
            index_AMPKFend = find(foldchange(index_minAMPKF:end,11) < 1.02,1,'first');
            
            % If responses do not return to baseline, set values to NaN so
            % they won't be plotted 
            if isempty(index_Gluconset) || isempty(index_AMPKFonset) || isempty(index_ATPonset) || isempty(index_Glucend) || isempty(index_ATPend) || isempty(index_AMPKFend)
                durationGluc(i,1) = NaN;
                durationATP(i,1) = NaN;
                durationAMPKF(i,1) = NaN;
                minATP(i,1) = NaN;
                minGluc(i,1) = NaN;
                minAMPKF(i,1) = NaN;
                fprintf('Simulation #%i: Not all species returned to baseline.\n',i);
            else
                % Recovery duration is TIME (tt) of end - time of onset -
                % 600s (= treatment duration, to allow direct comparison with
                % recovery duration of single-cells as there is a lot of
                % variation in treatment time for experiments
                fprintf('\n');
                durationGluc(i,1) = tt(index_minGluc + index_Glucend - 1) - tt(index_Gluconset) - 600;
                durationATP(i,1) = tt(index_minATP + index_ATPend - 1) - tt(index_ATPonset) - 600;
                durationAMPKF(i,1) = tt(index_minAMPKF + index_AMPKFend - 1) - tt(index_AMPKFonset) - 600;
                
                %Convert into minutes
                durationGluc(i,1) = durationGluc(i,1)/60;
                durationATP(i,1) = durationATP(i,1)/60;
                durationAMPKF(i,1) = durationAMPKF(i,1)/60;
            end
            
        end
             
     end                        % End of parameter variations
     
     RankMinATP = num2str(max(minATP) - min(minATP));
     RankMinAMPKF = num2str(max(minAMPKF) - min(minAMPKF));
     RankMinGluc = num2str(max(minGluc) - min(minGluc));
     RankDurATP = num2str(max(durationATP) - min(durationATP));
     RankDurAMPKF = num2str(max(durationAMPKF) - min(durationAMPKF));
     RankDurGluc = num2str(max(durationGluc) - min(durationGluc));
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%Generate Figures
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     switch a
         case 1 % Varying InitConc
             
            figure(h(1)), s(1) = subplot(1,12,x); 
            bar(minATP,1), title(strcat(nameGraph(x),' (',RankMinATP,')'));
            figure(h(2)), s(2) = subplot(3,4,x);
            bar(minAMPKF,1), title(strcat(nameGraph(x),' (',RankMinAMPKF,')'));
            figure(h(3)), s(3) = subplot(3,4,x);
            bar(minGluc,1), title(strcat(nameGraph(x),' (',RankMinGluc,')'));
            figure(h(4)), s(4) = subplot(3,4,x);
            bar(durationATP,1), title(strcat(nameGraph(x),' (',RankDurATP,')'));
            figure(h(5)), s(5) = subplot(3,4,x);
            bar(durationAMPKF,1), title(strcat(nameGraph(x),' (',RankDurAMPKF,')'));
            figure(h(6)), s(6) = subplot(3,4,x);
            bar(durationGluc,1), title(strcat(nameGraph(x),' (',RankDurGluc,')'));
            
            set(s(4),'YLim',[-1 7]);
            set(s(5),'YLim',[0 7.2]);
            set(s(6),'YLim',[0 50]);
            
         case 2 % Varying k_on
            plotPos = find(paramChangeID == x);
            figure(h(1)), s(1) = subplot(3,4,plotPos);
            bar(minATP,1), title(strcat(nameGraph(x),' (',RankMinATP,')'));
            figure(h(2)), s(2) = subplot(3,4,plotPos);
            bar(minAMPKF,1), title(strcat(nameGraph(x),' (',RankMinAMPKF,')'));
            figure(h(3)), s(3) = subplot(3,4,plotPos);
            bar(minGluc,1), title(strcat(nameGraph(x),' (',RankMinGluc,')'));
            figure(h(4)), s(4) = subplot(3,4,plotPos);
            bar(durationATP,1), title(strcat(nameGraph(x),' (',RankDurATP,')'));
            figure(h(5)), s(5) = subplot(3,4,plotPos);
            bar(durationAMPKF,1), title(strcat(nameGraph(x),' (',RankDurAMPKF,')'));
            figure(h(6)), s(6) = subplot(3,4,plotPos);
            bar(durationGluc,1), title(strcat(nameGraph(x),' (',RankDurGluc,')'));
            
            set(s(4),'YLim',[-1 7]);
            set(s(5),'YLim',[0 7.2]);
            set(s(6),'YLim',[0 50]);
    
         case 3 % Varying k_off
            plotPos = find(paramChangeID == x);
            figure(h(1)), s(1) = subplot(3,4,plotPos);
            bar(minATP,1), title(strcat(nameGraph(x),' (',RankMinATP,')'));
            figure(h(2)), s(2) = subplot(3,4,plotPos);
            bar(minAMPKF,1), title(strcat(nameGraph(x),' (',RankMinAMPKF,')'));
            figure(h(3)), s(3) = subplot(3,4,plotPos);
            bar(minGluc,1), title(strcat(nameGraph(x),' (',RankMinGluc,')'));
            figure(h(4)), s(4) = subplot(3,4,plotPos);
            bar(durationATP,1), title(strcat(nameGraph(x),' (',RankDurATP,')'));
            figure(h(5)), s(5) = subplot(3,4,plotPos);
            bar(durationAMPKF,1), title(strcat(nameGraph(x),' (',RankDurAMPKF,')'));
            figure(h(6)), s(6) = subplot(3,4,plotPos);
            bar(durationGluc,1), title(strcat(nameGraph(x),' (',RankDurGluc,')'));
            
            set(s(4),'YLim',[-1 7]);
            set(s(5),'YLim',[0 7.2]);
            set(s(6),'YLim',[0 50]);
     end
   end                         % End of all species / kinetics
end  

end

function necrosis(InitConc, k_on, k_off, varytemp)

global options ca_mag ca_duration ca_init

% Longer simulation time so as not to exclude results from simulations 
% with longer recovery times (investigating necrotic situation so there are
% species that do not return to baseline)
t_final = 20000;
stepsize = 1;               
tt = 0:stepsize:t_final;    

switch varytemp
  case 1, ca_mag_vary = [5:34 35:0.5:48.5 49:0.1:50 51:70];   
       num_sims = size(ca_mag_vary,2);     
       fprintf('Varying Magnitude of Calcium Influx: %i (au)\n\n',ca_mag)
  case 2, ca_duration_vary = [20:10:200 210:20:1200 1500:300:4200]; 
        % (Duration in seconds)
       num_sims = size(ca_duration_vary,2);
       fprintf('Varying Duration of Calcium Influx: %i s\n\n',ca_duration)
end

%pre-allocate matrices for speed
minCaC = zeros(num_sims,1);
minCaM = zeros(num_sims,1);
minAMP = zeros(num_sims,1);
minAMPK = zeros(num_sims,1);
minpAMPK = zeros(num_sims,1);
minGLUT3 = zeros(num_sims,1);
minGLUT3m = zeros(num_sims,1);
minGluc = zeros(num_sims,1);
minATP = zeros(num_sims,1);
minADP = zeros(num_sims,1);
minAMPKF = zeros(num_sims,1);
     
durationCaM = zeros(num_sims,1);
durationAMP = zeros(num_sims,1);
durationAMPK = zeros(num_sims,1);
durationpAMPK = zeros(num_sims,1);
durationGLUT3 = zeros(num_sims,1);
durationGLUT3m = zeros(num_sims,1);
durationGluc = zeros(num_sims,1);
durationATP = zeros(num_sims,1);
durationADP = zeros(num_sims,1);
durationAMPKF = zeros(num_sims,1);

for i = 1:num_sims       
% First simulation uses parameter values defined outside for loop
                                             
switch varytemp
    case 1, fprintf('Simulation %i: Solving... (ca-mag = %i) ',i,ca_mag)
    case 2, fprintf('Simulation %i: Solving... (ca-duration = %i) ',i,ca_duration)
end
              
% Calling ode solver
[t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,k_on,k_off,...
    ca_mag,ca_duration,ca_init); 
                
% Convert all curves to set number of data points (set by tt)
for j = 1:size(InitConc,2)
    values(:,j) = spline(t,y(:,j),tt);
end

% Calculate fold change for each species (model output)
for j = 1:size(InitConc,2)
   for m = 1:size(tt,2)
       foldchange(m,j) = (values(m,j)/values(1,j));%Baseline is first value
   end
end
             
% Calculate min/max and duration of all species
% (and absolute Ca value) curves (all called min*)
[minCaC(i,1), index_minCaC] = max(foldchange(:,1));
[minCaM(i,1), index_minCaM] = max(foldchange(:,2));
[minAMP(i,1), index_minAMP] = max(foldchange(:,3));
[minAMPK(i,1), index_minAMPK] = min(foldchange(:,4));
[minpAMPK(i,1), index_minpAMPK] = max(foldchange(:,5));
[minGLUT3(i,1), index_minGLUT3] = min(foldchange(:,6));
[minGLUT3m(i,1), index_minGLUT3m] = max(foldchange(:,7));
[minGluc(i,1), index_minGluc] = max(foldchange(:,8));
[minATP(i,1), index_minATP] = min(foldchange(:,9));
[minADP(i,1), index_minADP] = max(foldchange(:,10));
[minAMPKF(i,1), index_minAMPKF] = max(foldchange(:,11));
            
minCaC_abs(i,1) = max(values(:,1));      % Maximum/minimum absolute values
minATP_abs(i,1) = min(values(:,9));
minGluc_abs(i,1) = max(values(:,8));
minAMPKF_abs(i,1) = max(values(:,11));
            
% Find index of first value before min/max > 0.95 (response
% onset) and first value after (response end)
index_CaConset = find(foldchange(1:index_minCaC,1) < 1.02,1,'last');
index_CaMonset = find(foldchange(1:index_minCaM,2) < 1.02,1,'last');
index_AMPonset = find(foldchange(1:index_minAMP,3) < 1.02,1,'last');
index_AMPKonset = find(foldchange(1:index_minAMPK,4) > 0.98,1,'last');
index_pAMPKonset = find(foldchange(1:index_minpAMPK,5) < 1.02,1,'last');
index_GLUT3onset = find(foldchange(1:index_minGLUT3,6) > 0.98,1,'last');
index_GLUT3monset = find(foldchange(1:index_minGLUT3m,7) < 1.02,1,'last');
index_Gluconset = find(foldchange(1:index_minGluc,8) < 1.02,1,'last');
index_ATPonset = find(foldchange(1:index_minATP,9) > 0.98,1,'last');
index_ADPonset = find(foldchange(1:index_minADP,10) < 1.02,1,'last');
index_AMPKFonset = find(foldchange(1:index_minAMPKF,11) < 1.02,1,'last');
            
index_CaCend = find(foldchange(index_minCaC:end,1) < 1.02,1,'first')...
       + index_minCaC - 1;         % Need to get the ACTUAL index of end
index_CaMend = find(foldchange(index_minCaM:end,2) < 1.02,1,'first')...
       + index_minCaM - 1;
index_AMPend = find(foldchange(index_minAMP:end,3) < 1.02,1,'first')...
       + index_minAMP - 1;
index_AMPKend = find(foldchange(index_minAMPK:end,4) > 0.98,1,'first')...
       + index_minAMPK - 1;
index_pAMPKend = find(foldchange(index_minpAMPK:end,5) < 1.02,1,'first')...
       + index_minpAMPK - 1;
index_GLUT3end = find(foldchange(index_minGLUT3:end,6) > 0.9,1,'first')...     
       + index_minGLUT3 - 1;
index_GLUT3mend = find(foldchange(index_minGLUT3m:end,7) < 1.02,1,'first')...
       + index_minGLUT3m - 1;
index_Glucend = find(foldchange(index_minGluc:end,8) < 1.02,1,'first')...
       + index_minGluc - 1;
index_ATPend = find(foldchange(index_minATP:end,9) > 0.98,1,'first')...
       + index_minATP - 1;
index_ADPend = find(foldchange(index_minADP:end,10) < 1.02,1,'first')...
       + index_minADP - 1;
index_AMPKFend = find(foldchange(index_minAMPKF:end,11) < 1.02,1,'first')...
       + index_minAMPKF - 1;
            
% If responses do not return to baseline within simulation time, set  
% durations to NaN so they won't be plotted 
if isempty(index_CaConset) || isempty(index_CaCend) || ...
       isempty(index_CaMonset) || isempty(index_CaMend) || ...
       isempty(index_AMPonset) || isempty(index_AMPend) || ...
       isempty(index_AMPKonset) || isempty(index_AMPKend) || ...
       isempty(index_pAMPKonset) || isempty(index_pAMPKend) || ...
       isempty(index_GLUT3onset) || isempty(index_GLUT3end) || ...
       isempty(index_GLUT3monset) || isempty(index_GLUT3mend) || ...
       isempty(index_Gluconset) || isempty(index_Glucend) || ...
       isempty(index_ATPonset) || isempty(index_ATPend) || ...
       isempty(index_ADPonset) || isempty(index_ADPend) || ...
       isempty(index_AMPKFonset) || isempty(index_AMPKFend);
                    
   durationCaC(i,1) = NaN;
   durationCaM(i,1) = NaN;
   durationAMP(i,1) = NaN;
   durationAMPK(i,1) = NaN;
   durationpAMPK(i,1) = NaN;
   durationGLUT3(i,1) = NaN;
   durationGLUT3m(i,1) = NaN;
   durationGluc(i,1) = NaN;
   durationATP(i,1) = NaN;
   durationADP(i,1) = NaN;
   durationAMPKF(i,1) = NaN;
   fprintf('Simulation #%i: Not all species returned to baseline.\n',i);
else
   % Recovery duration is TIME (tt) of end - time of onset -
   % 600s (= treatment duration, to allow comparison with
   % recovery duration of single-cells as there is a lot of
   % variation in treatment time for experiments
   fprintf('\n');
   durationCaC(i,1) = tt(index_CaCend) - tt(round(ca_init)) - ca_duration;
   durationCaM(i,1) = tt(index_CaMend) - tt(round(ca_init)) - ca_duration;
   durationAMP(i,1) = tt(index_AMPend) - tt(round(ca_init)) - ca_duration;
   durationAMPK(i,1) = tt(index_AMPKend) - tt(round(ca_init)) - ca_duration;
   durationpAMPK(i,1) = tt(index_pAMPKend) - tt(round(ca_init)) - ca_duration;
   durationGLUT3(i,1) = tt(index_GLUT3end) - tt(round(ca_init)) - ca_duration;
   durationGLUT3m(i,1) = tt(index_GLUT3end) - tt(round(ca_init)) - ca_duration;
   durationGluc(i,1) = tt(index_Glucend) - tt(round(ca_init)) - ca_duration;
   durationATP(i,1) = tt(index_ATPend) - tt(round(ca_init)) - ca_duration;
   durationADP(i,1) = tt(index_ADPend) - tt(round(ca_init)) - ca_duration;
   durationAMPKF(i,1) = tt(index_AMPKFend) - tt(round(ca_init)) - ca_duration;
                
   %Convert into minutes to align with trace data
   durationCaC(i,1) = durationCaC(i,1)/60;
   durationCaM(i,1) = durationCaM(i,1)/60;
   durationAMP(i,1) = durationAMP(i,1)/60;
   durationAMPK(i,1) = durationAMPK(i,1)/60;
   durationpAMPK(i,1) = durationpAMPK(i,1)/60;
   durationGLUT3(i,1) = durationGLUT3(i,1)/60;
   durationGLUT3m(i,1) = durationGLUT3m(i,1)/60;
   durationGluc(i,1) = durationGluc(i,1)/60;
   durationATP(i,1) = durationATP(i,1)/60;
   durationADP(i,1) = durationADP(i,1)/60;
   durationAMPKF(i,1) = durationAMPKF(i,1)/60;
end
            
switch varytemp
       case 1, ca_mag = ca_mag_vary(i);                        
       case 2, ca_duration = ca_duration_vary(i);
end

y3d_ATP(:,i) = foldchange(:,9);
y3d_CaC(:,i) = foldchange(:,1);
time_3d(:,i) = tt;

end                                     % End of parameter variations
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Figures
%%% Plot max/min and recovery duration vs max Ca value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch varytemp
 case 1
     x = minCaC(:,1);
     xtitle = 'Ca_c Influx magnitude (Fold Change)';
     figure, plot(time_3d(1:4200,1:2:end)/60, y3d_CaC(1:4200,1:2:end));
     xlabel('Time (min)')
     ylabel('Cytosolic Calcium (fold change over baseline)')
     figure, plot(time_3d(1:4200,1:2:end)/60, y3d_ATP(1:4200,1:2:end));
     xlabel('Time (min)')
     ylabel('ATP Concentration (fold change over baseline)')
     
 case 2
     x = ca_duration_vary(1,:)/60;
     xtitle = 'Calcium influx duration (min)';
     figure, plot(time_3d(1:6000,1:2:end)/60, y3d_ATP(1:6000,1:2:end));
     xlabel('Time (min)')
     ylabel('ATP Concentration (fold change over baseline)')
end

figure, subplot(2,2,1)
[ax,h1,h2] = plotyy(x,durationATP(:,1),x,minATP(:,1),'scatter');               
set(h1,'Marker','.','MarkerEdgeColor','k')
set(h2,'Marker','.','MarkerEdgeColor','r')
set(get(ax(1),'Ylabel'),'String','ATP recovery duration (min)')
set(get(ax(2),'Ylabel'),'String','ATP min value (norm to baseline)')
xlabel(xtitle)
                
subplot(2,2,2)  
[ax,h1,h2] = plotyy(x,durationGluc(:,1),x,minGluc(:,1),'scatter');
set(h1,'Marker','.','MarkerEdgeColor','k')
set(h2,'Marker','.','MarkerEdgeColor','r')
set(get(ax(1),'Ylabel'),'String','Glucose recovery duration (min)')
set(get(ax(2),'Ylabel'),'String','Glucose max value (norm to baseline)')
xlabel(xtitle)
        
subplot(2,2,3)  
[ax,h1,h2] = plotyy(x,durationAMPKF(:,1),x,minAMPKF(:,1),'scatter');
set(h1,'Marker','.','MarkerEdgeColor','k')
set(h2,'Marker','.','MarkerEdgeColor','r')
set(get(ax(1),'Ylabel'),'String','AMPKAR recovery duration (min)')
set(get(ax(2),'Ylabel'),'String','AMPKAR max value (norm to baseline)')
xlabel(xtitle)
                
end

function runPopulation(InitConc, k_on, k_off, num_sims, temp_print)
global tt t_final options ca_init

vary_ca_mag = [34.5 60];
vary_ca_duration = [600 3600];
ca_vary_ = [0.1 0.1];

% Increase the randomness of the population variations
% rand('seed',input('Input seed for random number generator: '))

for j = 1:size(vary_ca_mag,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store original initial concentrations as InitConc variable will changd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InitConc_orig = InitConc;          
    k_on_orig = k_on;
    k_off_orig = k_off;
    ca_init_orig = ca_init;
    ca_mag_orig = vary_ca_mag(j);
    ca_duration_orig = vary_ca_duration(j);
    ca_vary = ca_vary_(j);     % % of calcium variation (see population fn)
    
    fprintf('\nca_mag = %i\n', vary_ca_mag(j))
    fprintf('ca_duration = %i\n', vary_ca_duration (j))
    
    for m = 1:3         % Repeat each 'experiment' 3 times: n = 3
        
        status = zeros(size(tt,2),num_sims);    % Pre-define status for speed
        
    for x = 1:num_sims  % # cells for each experiment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Vary initial concentrations and kinetics
        % Parameter variations fluctuate around original concentrations
        % and kinetics. Calcium magnitude, duration and time of onset
        % are also varied. Dependent kinetic parameters are updated
        % according to steady-state constraints
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if temp_print == 1
            fprintf('Cell %i...\n',x)
        end
        
        [InitConc,k_on,k_off,ca_mag,ca_duration,ca_init]...
          = population(InitConc_orig,k_on_orig,k_off_orig,ca_mag_orig,...
          ca_duration_orig,ca_init_orig,x,ca_vary); 
      
        %%%%%%%%%%%%%%%%%%%%
        %%% Run ODE solver
        %%%%%%%%%%%%%%%%%%%%
        [t,y] = ode15s(@ODEs_Run,[0,t_final],InitConc,options,k_on,k_off,...
          ca_mag,ca_duration,ca_init); 
      
        for i = 1:size(InitConc,2)
            values(:,i) = spline(t,y(:,i),tt);    %Spline so all solutions have same number of data points
        end
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine status of cell (viable/necrotic) for each timepoint
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nec_check = 0;                        
        nec_status_count = 0;
        
        for i = 1:size(tt,2)       % go through all time-points
            
            status(i,x) = 1;       % Default status is viable (Status = 1)
        
            if nec_check == 1
                status(i,x) = 3;   % Once necrotic, always necrotic (Status = 3)
            
                % If ATP levels drop below 0.35 * baseline cell is 
                % immediately necrotic 
            elseif values(i,9) < 0.35 * values(1,9)      
                status(i,x) = 3;                        
                nec_check = 1;
                if temp_print ==1
                    fprintf('is necrotic at time %i\n',i)
                end
            else
                % if ATP drops below 0.45*baseline, increment necrotic
                % status_count
                if values(i,9) >= 0.35*values(1,9) && values(i,9) < 0.45*values(1,9)      
                    nec_status_count = nec_status_count + 1;                            
                end

                % if ATP stays < 0.35 for > 500s, cell is necrotic
                if nec_status_count > 500                   
                    status(i,x) = 3;                        
                    nec_check = 1;
                    necr_time(x,1) = i;
                    if temp_print ==1
                        fprintf('is necrotic at time %i\n',i)
                    end
                end
            end
        
        end
        finalStatus(x) = status(end,x);
    end  
    
    viable = 0;     % Initialise counts
    necrotic = 0;
    
    for i = 1:size(finalStatus,2)        % For all cells
       if finalStatus(i) == 1
           viable = viable + 1;
       elseif finalStatus(i) == 3;
           necrotic = necrotic + 1;
       else
           fprintf('Error! at time %i, cell %i\n',tt,status)
       end
    end
    viab_percent = (viable/num_sims) * 100;
    necr_percent = (necrotic/num_sims) * 100;
    fprintf('Experiment %i (%i cells): %i percent viable, %i percent necrotic.\n',m,num_sims,viab_percent,necr_percent)
    end
end

end

function k_on = updateSS(InitConc,k_on,k_off)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Updates steady-state constraints when parameters are varied
%%% Returns updated kinetic parameters (k_on only, as determined by
%%% steady-state relationships)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CaC = InitConc(1);
CaM = InitConc(2);
AMP = InitConc(3);
AMPK = InitConc(4);
pAMPK = InitConc(5);
GLUT3 = InitConc(6);
GLUT3m = InitConc(7);
Gluc = InitConc(8);
ATP = InitConc(9);
ADP = InitConc(10);

k_on(1) = CaC*ATP*k_on(15);
k_on(2) = AMPK*k_off(2);
k_on(4) = (1.25*((ADP*Gluc*k_on(10))-(ATP*CaC*k_on(15))))/CaC;
k_on(5) = (1.25*((ADP*Gluc*k_on(10))-(ATP*CaC*k_on(15))))/((ATP^0.8)*CaM);
k_on(6) = (ATP*pAMPK*k_on(12))/(AMP*AMPK);
k_on(7) = ((-GLUT3*k_off(3)) + k_on(3))/(GLUT3*pAMPK);
k_on(8) = ((-GLUT3*k_off(3)) + k_on(3))/GLUT3m;   
k_on(9) = (ADP*Gluc*k_on(10) + Gluc*k_on(13))/(25*GLUT3m);
k_on(11) = (AMP*ATP*k_off(11))/ADP^2;
end

function [species_error,totalerror] = curvefit(foldchange,data,data_time,InitConc)
        
global stepsize
% # data points, # species
error = zeros(size(data_time,1),size(InitConc,2));           
species_error = zeros(1,size(InitConc,2));
        
errorTemp = 'abs';  % Set to abs for absolute value or chi for chi-squared
            
for i = 1:size(InitConc,2)                  %for each species 
   switch errorTemp
      case 'abs'
       for x = 1:size(data_time,1)          %go through each data point 
          if data(x,i)>0                    %if data exists
             %%calculate error (absolute value of difference between
             %%model and experiment)
             error(x,i) = abs(foldchange((data_time(x)*(1/stepsize)),i)-data(x,i));
          end
       end
                
     case 'chi'
       for x = 1:size(data_time,1)          %go through each data point 
          if data(x,i)>0                    %if data exists
             %%calculate error (chi-squared: square the difference between
             %%model and experiment). NB. Differences < 1 will
             %%be reduced by squaring them.....
             error(x,i) = (foldchange((data_time(x)*(1/stepsize)),i) - data(x,i))^2;
          end
       end
   end
        
end

% Define weighting function
% Set weighting function (glucose) to 0 for acute response
weightFn =  ones(size(data_time,1),size(InitConc,2));
weightFn(10:20,8) = 0;      
weightFn(:,1) = 0.5;          
%weightFn(:,9) = 0.5;
%weightFn(:,11) = 0.5;
    
% Multiply raw error function * weighting function
error = error.*weightFn;    
    
% Normalise error to experimental data (so 3 fitted to 4 will have same
% error as 0.75 fitted to 1). Gives NaN values for any elements in 
% data = 0 (can't divide by 0). Increases error for data < 1 (only issue
% for ATP).
error = error./data;        
    
for i = 1:size(InitConc,2)                     % for each species
   species_error(i) = nansum(error(:,i));       % sum the errors 
end
        
totalerror = sum(species_error);   % sum all errors for all the species
        
end

function [data_Glucmin data_ATPmin data_AMPKmin data_CaCmin... 
    data_Glucduration data_ATPduration data_AMPKduration... 
    data_CaCduration data_AMPKplot data_Glucplot data_ATPplot... 
    data_Fluo4plot] = getExptData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain data from single-cell traces - min/max, durations and plots
% .txt files must be saved in same folder as .m file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Min_values.txt contains min/max values of experimental single-cell
    %%% traces, and durations of recovery for each of Glucose, ATP and AMPK
    %%% Data obtained from "Single-Cell Data for curve fitting.xls"
    fileid = fopen('Min_values.txt');
    C = textscan(fileid,'%f %f %f %f %f %f %f %f');
    fclose(fileid);
            
    data_Glucmin(:,1) = C{1};                    % Glucose min
    data_ATPmin(:,1) = C{3};                     % ATP min
    data_AMPKmin(:,1) = C{5};                    % AMPK min
    data_CaCmin(:,1) = C{7};
    data_Glucduration(:,1) = C{2};
    data_ATPduration(:,1) = C{4};
    data_AMPKduration(:,1) = C{6};
    data_CaCduration(:,1) = C{8};
        
    %%% data_*.txt contain values of experimental single-cell traces that
    %%% can be plotted - this file caontains mean and upper and lower IQ
    %%% values for the trace averaged over all included experiments
    %%% Data obtained from "Single-Cell Data for curve fitting.xls"
    fileid = fopen('data_AMPKplot.txt');
    C = textscan(fileid,'%f %f %f %f');
    fclose(fileid);
    for m = 1:4
        %pAMPK data to plot TIME, median, IQ1, IQ2
        data_AMPKplot(:,m) = C{m};             
    end
                
    fileid = fopen('data_Glucplot.txt');
    C = textscan(fileid,'%f %f %f %f');
    fclose(fileid);
    for m = 1:4
        data_Glucplot(:,m) = C{m};             %Gluc data 
    end
            
    fileid = fopen('data_ATPplot.txt');
    C = textscan(fileid,'%f %f %f %f');
    fclose(fileid);
    for m = 1:4
       data_ATPplot(:,m) = C{m};             %ATP data 
    end
            
    fileid = fopen('data_Fluo4plot.txt');
    C = textscan(fileid,'%f %f %f %f');
    fclose(fileid);
    for m = 1:4
       data_Fluo4plot(:,m) = C{m};             %Fluo4 data
    end
end

function [InitConc_new,k_on,k_off,ca_mag,ca_duration,ca_init] = ...
    population(InitConc_orig,k_on_orig,k_off_orig,ca_mag_orig,...
    ca_duration_orig,ca_init_orig,x,ca_vary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function generates a set of initial conditions for a cell.
%%% New initial concentrations are based on original
%%% concentrations read into program and alters them within X% based
%%% on steady-state relationships defined from Mathematica analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check = 1; 
count = 0;    
    
k_on = k_on_orig;   % Default values are original values 
k_off = k_off_orig; 
    
% Vary calcium input (modelling cell-to-cell variability in glutamate
% receptors). 
% For cell-to-cell heterogenetiy, we are varying around cells that 
% SURVIVED only (i.e. ATP only recorded if cells recovered), so limit 
% variation here to +/-10% (varyCm).
% For population model, where necrosis occurs, vary this +/- 0.85. 
    
varyCm = ca_vary;
varyC = 0.1;    % Vary ca_duration by +/-10%
    
ca_mag = ca_mag_orig + (ca_mag_orig * (-varyCm + (varyCm + varyCm) * rand(1)));  
ca_duration = ca_duration_orig + (ca_duration_orig * (-varyC + (varyC + varyC) * rand(1)));   
    
%ca_init = ca_init_orig + (ca_init_orig * (-varyC + (varyC + varyC) * rand(1)));
ca_init = ca_init_orig; %(To align traces)
         
rand(1);
rand(1);

while check == 1
    % Track the number of iterations required to get non-negative values
    count = count + 1;      
    if x == 1               % x is simulation number
       Cac = InitConc_orig(1); % First simulation uses original values
       Cam = InitConc_orig(2);
       AMP = InitConc_orig(3);
       AMPK = InitConc_orig(4);
       pAMPK = InitConc_orig(5);
       Glut3 = InitConc_orig(6);
       Glut3m = InitConc_orig(7);
       Gluc = InitConc_orig(8);        
       ATP = InitConc_orig(9);
       ADP = InitConc_orig(10);
       pAMPKFret = InitConc_orig(11);
            
    else              % Remaining simulations use varied values
       
        % Vary concentrations within +/- 20% of original values
       varyI = 0.20;           
       % Vary kinetics within +/- 20% of original values
       varyK = 0.20;
            
       Cac = InitConc_orig(1) + (InitConc_orig(1) * (-varyI + (varyI + varyI) * rand(1)));    
       Cam = InitConc_orig(2) + (InitConc_orig(2) * (-varyI + (varyI + varyI) * rand(1)));    
       AMP = InitConc_orig(3) + (InitConc_orig(3) * (-varyI + (varyI + varyI) * rand(1)));    
       AMPK = InitConc_orig(4) + (InitConc_orig(4) * (-varyI + (varyI + varyI) * rand(1)));    
       pAMPK = InitConc_orig(5) + (InitConc_orig(5) * (-varyI + (varyI + varyI) * rand(1)));    
       Glut3 = InitConc_orig(6) + (InitConc_orig(6) * (-varyI + (varyI + varyI) * rand(1)));    
       Glut3m = InitConc_orig(7) + (InitConc_orig(7) * (-varyI + (varyI + varyI) * rand(1)));    
       Gluc = InitConc_orig(8) + (InitConc_orig(8) * (-varyI + (varyI + varyI) * rand(1)));    
       ATP = InitConc_orig(9) + (InitConc_orig(9) * (-varyI + (varyI + varyI) * rand(1)));
       ADP = InitConc_orig(10) + (InitConc_orig(10) * (-varyI + (varyI + varyI) * rand(1))); 
       pAMPKFret = pAMPK;
            
       for i = [3 10 12 13 14 15]
          k_on(i) = k_on_orig(i) + (k_on_orig(i) * (-varyK + (varyK + varyK) * rand(1)));       
       end
            
       for i = [2 3 11]
          k_off(i) = k_off_orig(i) + (k_off_orig(i) * (-varyK + (varyK + varyK) * rand(1)));
       end
                         
    end

    InitConc_new(1) = Cac;
    InitConc_new(2) = Cam;
    InitConc_new(3) = AMP;
    InitConc_new(4) = AMPK;
    InitConc_new(5) = pAMPK;
    InitConc_new(6) = Glut3;
    InitConc_new(7) = Glut3m;
    InitConc_new(8) = Gluc;
    InitConc_new(9) = ATP;
    InitConc_new(10) = ADP;
    InitConc_new(11) = pAMPKFret;         
                   
    %%%%%%%%%%%%%%%%%%
    %%% Maintain steady-state constraints for remaining parameters
    %%%%%%%%%%%%%%%%%%
    k_on = updateSS(InitConc_new,k_on,k_off);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    check = 0;
         
    %Ensuring no negative k values get through
    for m = 1:size(k_on,2)          
      if k_on(m) < 0, check = 1; end
      if k_off(m) < 0, check = 1; end
    end
        
end        
end

function [InitConc kon koff] = defineParam
InitConc = [180 320 24000 144 81 400 233 112500000 2700000 200000 81];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CaC = InitConc(1);
CaM = InitConc(2);
AMP = InitConc(3);
AMPK = InitConc(4);
pAMPK = InitConc(5);
GLUT3 = InitConc(6);
GLUT3m = InitConc(7);
Gluc = InitConc(8);
ATP = InitConc(9);
ADP = InitConc(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables - independent kinetic parameters
% GLUT3m and Glucose degration rates (listed as koff8 and koff13 in the 
% publication, for clarity) are named kon8 and kon13 here and throughout 
% the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_AMPK = 360;       % Half-lives (min)
t_GLUT3 = 900;
t_Gluc = 1.24/25;

koff(2) = log(2)/(t_AMPK*60);
koff(3) = log(2)/(t_GLUT3*60);
kon(3) = koff(3)*GLUT3+0.12;
kon(10) = 2.7e-010;
koff(11) = 4.5e-008;
kon(12) = 0.02;
kon(13) = log(2)/(t_Gluc*60);
kon(14) = 0.17;
kon(15) = 4.0e-008; % Referred to as kon1a in manuscript
koff(12) = 0;       % koff and kon vectors must have same dimensions
koff(13) = 0;
koff(14) = 0;
koff(15) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES - steady-state constraints (dependent kinetic parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kon(1) = CaC*ATP*kon(15);
kon(2) = AMPK*koff(2);
kon(4) = (1.25*((ADP*Gluc*kon(10))-(ATP*CaC*kon(15))))/CaC;
kon(5) = (1.25*((ADP*Gluc*kon(10))-(ATP*CaC*kon(15))))/((ATP^0.8)*CaM);
kon(6) = (ATP*pAMPK*kon(12))/(AMP*AMPK);
kon(7) = ((-GLUT3*koff(3)) + kon(3))/(GLUT3*pAMPK);
kon(8) = ((-GLUT3*koff(3)) + kon(3))/GLUT3m;    % Referred to as koff8 in ms
kon(9) = (ADP*Gluc*kon(10) + Gluc*kon(13))/(25*GLUT3m);
kon(11) = (AMP*ATP*koff(11))/ADP^2;
end

function [dydt] = ODEs_Run(t,y,k_on,k_off,ca_mag,ca_duration,ca_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CaC = y(1);
CaM = y(2);
AMP = y(3);
AMPK = y(4);
pAMPK = y(5);
GLUT3 = y(6);
GLUT3m = y(7);
Gluc = y(8);
ATP = y(9);
ADP = y(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTION KINETICS / FLUXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1 = k_on(1);
R1a = k_on(15)*CaC*ATP;
R2 = k_on(2);
R2a = k_off(2)*AMPK;
R3 = k_on(3);
R3a = k_off(3)*GLUT3;
R4 = k_on(4)*CaC;
R5 = k_on(5)*CaM*(ATP^0.8);
R6 = k_on(6)*AMP*AMPK;
R7 = k_on(7)*pAMPK*GLUT3;
R8 = k_on(8)*GLUT3m;
R9 = k_on(9)*GLUT3m;
R10 = k_on(10)*Gluc*ADP;
R11 = k_on(11)*ADP*ADP;
R11a = k_off(11)*ATP*AMP;
R12 = k_on(12)*ATP*pAMPK;
R13 = k_on(13)*Gluc;    % Referred to as koff13 in manuscript (for clarity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt = double(y);
dydt(1) = ca_mag * (t > ca_init & t < (ca_init+ca_duration)) +R1-R1a-R4+R5;
dydt(2) = +R4-R5;
dydt(3) = +R11-R11a;
dydt(4) = +R2-R2a-R6+R12;
dydt(5) = +R6-R12;
dydt(6) = +R3-R3a-R7;
dydt(7) = +R7-R8;
dydt(8) = +25*R9-R10-R13;
dydt(9) = -R1a-(0.8*R5)+R10+R11-R11a;
dydt(10) = R1a+(0.8*R5)-R10-2*R11+2*R11a;
dydt(11) = k_on(14) * (+R6 - R12);

end