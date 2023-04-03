%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DietrichScript: VM
% Date created: 23.04.22
% Date last mostified: 22.07.22
% Purpose: To test the implementation of the Dietrich drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
cd('C://Users/roisi/Desktop/mP Model Code/')
% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
VM_Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

rho_p = table2array(VM_Dataset(:, "ParticleDensity"));
rho_f = table2array(VM_Dataset(:, "FluidDensity"));
vis_dyn = table2array(VM_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(VM_Dataset(:, "KinematicViscosity"));

d_equi = table2array(VM_Dataset(:, "ParticleSize"));
size_a = table2array(VM_Dataset(:, "a"));
size_b = table2array(VM_Dataset(:, "b"));
size_c = table2array(VM_Dataset(:, "c"));
shape = table2array(VM_Dataset(:, "Shape"));

shape_flt = table2array(VM_Dataset(:, "Flatness"));
shape_eln = table2array(VM_Dataset(:, "elongation"));
shape_del = table2array(VM_Dataset(:, "Dellino"));
shape_sph = table2array(VM_Dataset(:, "Sphericity"));
shape_cir = table2array(VM_Dataset(:, "Circularity"));
Reynolds = table2array(VM_Dataset(:, "Re"));
Powers = table2array(VM_Dataset(:, "Powers"));

wvel_meas = table2array(VM_Dataset(:, "Wmeasured"));

% Set up and calculate additional variables:
SA_mP = zeros(140, 1);
SA_EqSph = zeros(140, 1);
Vol_mP = zeros(140, 1);
Mass_mP = zeros(140, 1);
CSF = zeros(140, 1);
rho_rel = zeros(140, 1);
ProjA_ESD = zeros(140, 1);

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dietrich
% <<<<<<<<<<<<<<<<<
% Dietrich's model cannot be used when CSF<0.2
% This is not an interative procedure, it just calculates the terminal velocity.

%% TIMED SECTION:
% Going to create a loop to measure the time taken 100 times to ensure there
% is no variation. Note that it is only the calculation steps that are
% measured.
T_Start = zeros(1, 1000);
T_End = zeros(1, 1000);
T_CPU = zeros(1, 1000);
T_Elapsed = zeros(1, 1000);

% Start time for 1000 loops
T_1000Times_I = cputime;
tStart = tic;
for n = 1:1000
    % Set up variable arrays
    d_dim = zeros(140, 1);
    wdim_Dietrich = zeros(140, 1);
    wt_Dietrich = zeros(140, 1);
    R1 = zeros(140,1);
    R2 = zeros(140,1);
    R3 = zeros(140,1);
    g=9.81;
    % Take CPUTIME:
    T_Start(n) = cputime;
    % Take clock time
    tStart_Loop = tic;
    for i=1:140
        d_dim(i) = ((rho_p(i) - rho_f(i))*g*(d_equi(i)^3.0))/(rho_f(i)*(vis_kin(i)^2.0));
	
	    R1(i) = -3.76715 + 1.92944*(log10(d_dim(i))) - 0.09815*((log10(d_dim(i)))^2.0) ...
		    -0.00575*((log10(d_dim(i)))^3.0) + 0.00056*((log10(d_dim(i)))^4.0);
	    
        if(CSF(i)>0.15)
            R2(i) = (log10(1-((1-CSF(i)))/0.85)) - ((1-CSF(i))^2.3)*tanh(log10(d_dim(i))-4.6)...
		    + 0.3*(0.5-CSF(i))*((1-CSF(i))^2.0)*(log10(d_dim(i))-4.6);
        else
            R2(i)=nan;
        end
	    R3(i) = (0.65-((CSF(i)/2.83)*(tanh(log10(d_dim(i))-4.6))))^(1+((3.5-Powers(i))/2.5));
	    
	    wdim_Dietrich(i) = R3(i) * (10^(R1(i)+R2(i)));
	    
	    wt_Dietrich(i) = ((wdim_Dietrich(i)*(rho_p(i) - rho_f(i))*g*vis_kin(i))/rho_f(i))^(1.0/3.0);
    end
    T_End(n) = cputime;
    T_CPU(n) = T_End(n) -T_Start(n); 
    tEnd_Loop(n) = toc(tStart_Loop);
end
tEnd = toc(tStart);
T_1000Times_N = cputime;
T_1000Times = T_1000Times_N - T_1000Times_I;
T_1Time = T_1000Times/100.0;

%% Store output in one array

Results_Dietrich = zeros(1, 4);
Results_Dietrich(1) = T_1000Times;
Results_Dietrich(2) = mean(T_CPU);
Results_Dietrich(3) = tEnd;
Results_Dietrich(4) = mean(tEnd_Loop);

Results_Dietrich = transpose(Results_Dietrich);

Table_Dietrich_Timing = array2table(Results_Dietrich, "RowNames", ...
    {'CPU: 1000 tests', 'Average CPU Per test', 'Time for 1000 tests:', 'Average time for 1 test'});

%% Going to repeat the test three timee
% The goal is to determine whether running the script as soon as MATLAB is
% started will have an impact on the result.

% How to do this:
% Test 1: Close Matlab completeley and restart. Uncomment the line below
% for the file ending in '1' and run the programme.
% Test 2: After Test 1, uncomment the line below for the file endingin '2'
% and run the programme.
% Test 3: After Test 2, uncomment the line below for the file endingin '2'
% and run the programme.

% Ensure that when running the test, no other scripts are open on Matlab.

% % Test 1
% writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing1.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing1.xls', 'WriteRowNames', true);

% %Test 2
% writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing2.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing2.xls', 'WriteRowNames', true);

%Test 3
writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing3.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich_Timing, './DragModelsTest/Output/20230215/Dietrich_Timing3.xls', 'WriteRowNames', true);

