%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DioguardiScript: VM
% Date created: 15.02.23
% Date last mostified: 15.02.23
% Purpose: To estimate the computational time need to run the code for
%           the implementation of the Bagheri drag model on a range of
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
Cd_meas = table2array(VM_Dataset(:, 'CdMeasured'));

% Set up and calculate additional variables:
SA_mP = zeros(140, 1);
SA_EqSph = zeros(140, 1);
Vol_mP = zeros(140, 1);
Mass_mP = zeros(140, 1);
CSF = zeros(140, 1);
rho_rel = zeros(140, 1);
ProjA_ESD = zeros(140, 1);
g=9.81;

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dioguardi' method 1
% <<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

%% TIMED SECTION:
% Going to create a loop to measure the time taken 1000 times to ensure there
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
    Re_Dio = zeros(140, 10000);
    ReFinal_Dio = zeros(140, 1);
    CdFinal_Dio = zeros(140, 1);
    Cd_Dio = zeros(140, 10000);
    Fd_Dio = zeros(140, 10000);
    Fg_Dio = zeros(140, 10000);
    Fb_Dio = zeros(140, 10000);
    Fnet_Dio = zeros(140, 10000);
    Dist1_Dio = zeros(140, 10000);
    DistTot_Dio=zeros(140, 1);
    Acc_Dio = zeros(140, 10000);	
    FinalTime_Dio = zeros(140, 1);
    FinalStep_Dio = zeros(140, 1);
    wtFinal_Dio = zeros(140, 1);

    % Set timestep
    timestep = 0.0002;
    
    % Set initial velocity and timestep
    wvel_Dio = zeros(140, 10000);
    wvel_Dio(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                                  % terminal velocity is independent of the
                                  % initial velocity.
    % Take CPUTIME:
    T_Start(n) = cputime;
    % Take clock time
    tStart_Loop = tic;
    % Begin calculation
    for i=1:140
        for t=1:10000
      
            Re_Dio(i,t) = abs((rho_p(i) * wvel_Dio(i, t) * d_equi(i))/ vis_dyn(i));
		    
		    Cd_Dio(i, t) = (24.0/Re_Dio(i,t))*(((1.0-shape_del(i))/(Re_Dio(i,t)+1.0))^0.25) ...
			         + (24.0/Re_Dio(i,t))*0.1806*(Re_Dio(i,t)^0.6459)*(shape_del(i)^(-1.0*(Re_Dio(i,t)^0.08))) ...
			         + 0.4251/(1.0+((6880.95/Re_Dio(i,t))*(shape_del(i)^5.05)));
	    
		    Fd_Dio(i,t) = 0.5*rho_f(i)*SA_mP(i)*(abs(wvel_Dio(i,t))*wvel_Dio(i,t))*Cd_Dio(i,t);
	    
		    Fg_Dio(i,t) = Vol_mP(i)*rho_p(i)*g;
	    
		    Fb_Dio(i,t) = Vol_mP(i)*rho_f(i)*g;
	    
		    Fnet_Dio(i,t) = Fg_Dio(i,t) - Fb_Dio(i,t) - Fd_Dio(i,t);
	    
		    wvel_Dio(i,t+1) = ((Fnet_Dio(i,t)/Mass_mP(i))*timestep)+wvel_Dio(i,t);
	    
            Dist1_Dio(i,t) = wvel_Dio(i,t) * timestep;
		    DistTot_Dio(i) = DistTot_Dio(i) + Dist1_Dio(i,t);
    
		    Acc_Dio(i,t) = (wvel_Dio(i,t+1) - wvel_Dio(i,t))/timestep;
		    
            if (Acc_Dio(i,t)< 0.001)
                FinalTime_Dio(i) = (t+1)*timestep;
                FinalStep_Dio(i) = t+1;
                wtFinal_Dio(i)=wvel_Dio(i, t+1);
                Dist1_Dio(i, t+1) = wvel_Dio(i, t+1) * timestep;
                DistTot_Dio(i) = DistTot_Dio(i) + Dist1_Dio(i, t+1);
                ReFinal_Dio(i) = abs((rho_p(i) * wvel_Dio(i, t+1) * d_equi(i))/ vis_dyn(i));
		        CdFinal_Dio(i) = (24.0/ReFinal_Dio(i))*(((1.0-shape_del(i))/(ReFinal_Dio(i)+1.0))^0.25) ...
			         + (24.0/ReFinal_Dio(i))*0.1806*(ReFinal_Dio(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_Dio(i)^0.08))) ...
			         + 0.4251/(1.0+((6880.95/ReFinal_Dio(i))*(shape_del(i)^5.05)));
                break
            end
        end
    end
     T_End(n) = cputime;
    T_CPU(n) = T_End(n) -T_Start(n); 
    tEnd_Loop(n) = toc(tStart_Loop);
end
tEnd = toc(tStart);
T_1000Times_N = cputime;
T_1000Times = T_1000Times_N - T_1000Times_I;
T_1Time = T_1000Times/100.0;

%% Store output:

Results_Dio_SA = zeros(1, 4);
Results_Dio_SA(1) = T_1000Times;
Results_Dio_SA(2) = mean(T_CPU);
Results_Dio_SA(3) = tEnd;
Results_Dio_SA(4) = mean(tEnd_Loop);

Results_Dio_SA = transpose(Results_Dio_SA);

%% Dioguardi Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area in the calculation of the drag force.

%% TIMED SECTION:
% Going to create a loop to measure the time taken 5 times to ensure there
% is no variation. Note that it is only the calculation steps that are
% measured.
T_Start = zeros(1, 1000);
T_End = zeros(1, 1000);
T_CPU = zeros(1, 1000);
T_Elapsed = zeros(1, 1000);

% Start time for 100 loops
T_1000Times_I = cputime;
tStart = tic;
for n = 1:1000
    Re_Dio = zeros(140, 10000);
    ReFinal_Dio = zeros(140, 1);
    Cd_Dio = zeros(140, 10000);
    CdFinal_Dio = zeros(140, 1);
    Fd_Dio = zeros(140, 10000);
    Fg_Dio = zeros(140, 10000);
    Fb_Dio = zeros(140, 10000);
    Fnet_Dio = zeros(140, 10000);
    Dist1_Dio = zeros(140, 10000);
    DistTot_Dio=zeros(140, 1);
    Acc_Dio = zeros(140, 10000);	
    FinalTime_Dio = zeros(140, 1);
    FinalStep_Dio = zeros(140, 1);
    wtFinal_Dio = zeros(140, 1);
    
    % Set timestep
    timestep = 0.0002;
    
    % Set initial velocity and timestep
    wvel_Dio = zeros(140, 10000);
    wvel_Dio(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                                  % terminal velocity is independent of the
                                  % initial velocity.
    % Take CPUTIME:
    T_Start(n) = cputime;
    % Take clock time
    tStart_Loop = tic;                             
    % Begin calculation
    for i=1:140
        for t=1:10000
            
            Re_Dio(i,t) = abs((rho_p(i) * wvel_Dio(i, t) * d_equi(i))/ vis_dyn(i));
		    
		    Cd_Dio(i, t) = (24.0/Re_Dio(i,t))*(((1.0-shape_del(i))/(Re_Dio(i,t)+1.0))^0.25) ...
			         + (24.0/Re_Dio(i,t))*0.1806*(Re_Dio(i,t)^0.6459)*(shape_del(i)^(-1.0*(Re_Dio(i,t)^0.08))) ...
			         + 0.4251/(1.0+((6880.95/Re_Dio(i,t))*(shape_del(i)^5.05)));
	    
		    Fd_Dio(i,t) = 0.5*rho_f(i)*ProjA_ESD(i)*(wvel_Dio(i,t)^2.0)*Cd_Dio(i,t);
	    
		    Fg_Dio(i,t) = Vol_mP(i)*rho_p(i)*g;
	    
		    Fb_Dio(i,t) = Vol_mP(i)*rho_f(i)*g;
	    
		    Fnet_Dio(i,t) = Fg_Dio(i,t) - Fb_Dio(i,t) - Fd_Dio(i,t);
	    
		    wvel_Dio(i,t+1) = ((Fnet_Dio(i,t)/Mass_mP(i))*timestep)+wvel_Dio(i,t);
	    
		    Dist1_Dio(i,t) = wvel_Dio(i,t) * timestep;
		    DistTot_Dio(i) = DistTot_Dio(i) + Dist1_Dio(i,t);
    
		    Acc_Dio(i,t) = (wvel_Dio(i,t+1) - wvel_Dio(i,t))/timestep;
		    
            if (Acc_Dio(i,t)< 0.001)
                FinalTime_Dio(i) = (t+1)*timestep;
                FinalStep_Dio(i) = t+1;
                wtFinal_Dio(i)=wvel_Dio(i, t+1);
                Dist1_Dio(i, t+1) = wvel_Dio(i, t+1) * timestep;
                DistTot_Dio(i) = DistTot_Dio(i) + Dist1_Dio(i, t+1);
                ReFinal_Dio(i) = abs((rho_p(i) * wvel_Dio(i, t+1) * d_equi(i))/ vis_dyn(i));
		        CdFinal_Dio(i) = (24.0/ReFinal_Dio(i))*(((1.0-shape_del(i))/(ReFinal_Dio(i)+1.0))^0.25) ...
			         + (24.0/ReFinal_Dio(i))*0.1806*(ReFinal_Dio(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_Dio(i)^0.08))) ...
			         + 0.4251/(1.0+((6880.95/ReFinal_Dio(i))*(shape_del(i)^5.05)));
                break
            end
        end
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

Results_Dio_Proj = zeros(1, 4);
Results_Dio_Proj(1) = T_1000Times;
Results_Dio_Proj(2) = mean(T_CPU);
Results_Dio_Proj(3) = tEnd;
Results_Dio_Proj(4) = mean(tEnd_Loop);

Results_Dio_Proj = transpose(Results_Dio_Proj);

%% Combine results arrays:
All_Results_Dio = [Results_Dio_SA Results_Dio_Proj];

%%
Table_Dio_Timing = array2table(All_Results_Dio, "RowNames", ...
    {'CPU: 1000 tests', 'Average CPU Per test', 'Time for 1000 tests:', 'Average time for 1 test'}, ...
    'VariableNames', {'SA', 'Proj'} );

%% Going to repeat the test three timee
% The goal is to determine whether running the script as soon as MATLAB is
% started will have an impact on the result.

% How to do this:
% Test 1: Close Matlab completeley and restart. Uncomment the line below
% for the file ending in '1' and run the programme.
% Test 2: After Test 1, uncomment the line below for the file endingin '2'
% and run the programme.
% Test 3: After Test 2, uncomment the line below for the file endingin '3'
% and run the programme.

% Ensure that when running the test, no other scripts are open on Matlab.

% % Test 1
% writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing1.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing1.xls', 'WriteRowNames', true);

% % Test 2
% writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing2.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing2.xls', 'WriteRowNames', true);

% % Test 3
writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing3.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dio_Timing, './DragModelsTest/Output/20230215/Dio_Timing3.xls', 'WriteRowNames', true);

