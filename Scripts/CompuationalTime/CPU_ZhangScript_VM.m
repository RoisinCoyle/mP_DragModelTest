%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: CPU_ZhangScript: VM
% Date created: 15.02.23
% Date last mostified: 15.02.23
% Purpose: To estimate the computational time need to run the code for
%           the implementation of the Bagheri drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
cd('C:\\Users\roisi\Desktop\mP Model Code\')
% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
VM_Dataset = readtable("C:\Users\roisi\Desktop\mP Model Code\SettlingVelocity calc\VanMelkebekeSIDataset.txt");

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
Cd_meas = table2array(VM_Dataset(:, "CdMeasured"));

% Set up and calculate additional variables:
SA_mP = zeros(140, 1);
SA_EqSph = zeros(140, 1);
Vol_mP = zeros(140, 1);
Mass_mP = zeros(140, 1);
CSF = zeros(140, 1);
rho_rel = zeros(140, 1);
ProjA_ESD = zeros(140, 1);
Zhang_EstVol = zeros(140, 1);
Zhang_deq = zeros(140, 1);
Zhang_Deqv = zeros(140, 1);
ZhangProjA = zeros(140, 1);
Zhang_EstMass = zeros(140, 1);
g=9.81;

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
    Zhang_EstVol(i) = size_a(i)*size_b(i)*size_c(i);
    Zhang_Deqv(i) = ((6*Zhang_EstVol(i))/pi())^(1/3);
    Zhang_deq(i) = sqrt((4*size_a(i)*size_b(i))/pi());
    ZhangProjA(i) = pi()*(Zhang_deq(i)^2)*0.25;
    Zhang_EstMass(i) = Zhang_EstVol(i)*rho_p(i);
end

%% Zhang and Choi's method 1
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    % Set timestep
    timestep = 0.0001;
    
    % Set up variable arrays
    Re_ZC = zeros(140, 10000);
    ReFinal_ZC = zeros(140, 1);
    shape_ASF = zeros(140, 1);
    wvel_ZC = zeros(140, 10000);
    Cd_ZC = zeros(140, 10000);
    CdFinal_ZC = zeros(140, 1);
    Fd_ZC = zeros(140, 10000);
    Fg_ZC = zeros(140, 10000);
    Fb_ZC = zeros(140, 10000);
    Fnet_ZC = zeros(140, 10000);
    Dist1_ZC = zeros(140, 10000);
    DistTot_ZC = zeros(140, 1);
    Acc_ZC = zeros(140, 10000);
    FinalTime_ZC = zeros(140, 1);
    FinalStep_ZC = zeros(140, 1);
    wtFinal_ZC = zeros(140, 1);
    
    % Set initial velocity and timestep
    wvel_ZC(:, 1) = 0.0001; % Note that earlier tests have shown that the
                            % terminal velocity is independent of the
                            % initial velocity.
    % Take CPUTIME:
    T_Start(n) = cputime;
    % Take clock time
    tStart_Loop = tic;
    % Begin calculation
    for i=1:140	
        
        shape_ASF(i) = (size_a(i) * size_c(i))/(size_b(i)^2);
        
        for t=1:10000
		    
		    Re_ZC(i, t) = abs((rho_p(i) * wvel_ZC(i, t) * Zhang_deq(i))/ vis_dyn(i));
		    
		    Cd_ZC(i,t) = (58.58*(shape_ASF(i)^0.1936))/(Re_ZC(i,t)^0.8273);
	    
		    Fd_ZC(i,t) = 0.5*rho_f(i)*SA_mP(i)*(abs(wvel_ZC(i,t))*wvel_ZC(i,t))*Cd_ZC(i,t);
	    
		    Fg_ZC(i,t) = Vol_mP(i)*rho_p(i)*g;
	    
		    Fb_ZC(i,t) = Vol_mP(i)*rho_f(i)*g;
	    
		    Fnet_ZC(i,t) = Fg_ZC(i,t) - Fb_ZC(i,t) - Fd_ZC(i,t);
	    
		    wvel_ZC(i,t+1) = ((Fnet_ZC(i,t)/Mass_mP(i))*timestep)+wvel_ZC(i,t);
    
            Dist1_ZC(i,t) = wvel_ZC(i,t) * timestep;
		    DistTot_ZC(i) = DistTot_ZC(i) + Dist1_ZC(i,t);
		    Acc_ZC(i,t) = (wvel_ZC(i, t+1) - wvel_ZC(i,t))/timestep;
		    
            if (Acc_ZC(i,t)< 0.001)
                FinalTime_ZC(i) = (t+1)*timestep;
                FinalStep_ZC(i) = t+1;
                wtFinal_ZC(i)=wvel_ZC(i, t+1);
                Dist1_ZC(i, t+1) = wvel_ZC(i, t+1) * timestep;
                DistTot_ZC(i) = DistTot_ZC(i) + Dist1_ZC(i, t+1);
                ReFinal_ZC(i) = abs((rho_p(i) * wvel_ZC(i, t+1) * d_equi(i))/ vis_dyn(i));
		        CdFinal_ZC(i) = (58.58*(shape_ASF(i)^0.1936))/(ReFinal_ZC(i)^0.8273);
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

Results_ZC_SA = zeros(1, 4);
Results_ZC_SA(1) = T_1000Times;
Results_ZC_SA(2) = mean(T_CPU);
Results_ZC_SA(3) = tEnd;
Results_ZC_SA(4) = mean(tEnd_Loop);

Results_ZC_SA = transpose(Results_ZC_SA);

%% Zhang Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area, using Newtons Drag formula

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
    % Set timestep
    timestep = 0.00015;
    
    % Set up variable arrays
    Re_ZC = zeros(140, 10000);
    shape_ASF = zeros(140, 1);
    wvel_ZC = zeros(140, 10000);
    Cd_ZC = zeros(140, 10000);
    Fd_ZC = zeros(140, 10000);
    Fg_ZC = zeros(140, 10000);
    Fb_ZC = zeros(140, 10000);
    Fnet_ZC = zeros(140, 10000);
    Dist1_ZC = zeros(140, 10000);
    DistTot_ZC = zeros(140, 1);
    Acc_ZC = zeros(140, 10000);
    FinalTime_ZC = zeros(140, 1);
    FinalStep_ZC = zeros(140, 1);
    wtFinal_ZC = zeros(140, 1);
    
    % Set initial velocity and timestep
    wvel_ZC(:, 1) = 0.0001; % Note that earlier tests have shown that the
                            % terminal velocity is independent of the
                            % initial velocity.
    % Take CPUTIME:
    T_Start(n) = cputime;
    % Take clock time
    tStart_Loop = tic;
    % Begin calculation
    for i=1:140	
        
        shape_ASF(i) = (size_a(i) * size_c(i))/(size_b(i)^2);
        
        for t=1:10000
		    
		    Re_ZC(i, t) = abs((rho_p(i) * wvel_ZC(i, t) * Zhang_deq(i))/ vis_dyn(i));
		    
		    Cd_ZC(i,t) = (58.58*(shape_ASF(i)^0.1936))/(Re_ZC(i,t)^0.8273);
	    
		    Fd_ZC(i,t) = 0.5*rho_f(i)*ZhangProjA(i)*(abs(wvel_ZC(i,t))*wvel_ZC(i,t))*Cd_ZC(i,t);
	    
		    Fg_ZC(i,t) = Vol_mP(i)*rho_p(i)*g;
	    
		    Fb_ZC(i,t) = Vol_mP(i)*rho_f(i)*g;
	    
		    Fnet_ZC(i,t) = Fg_ZC(i,t) - Fb_ZC(i,t) - Fd_ZC(i,t);
	    
		    wvel_ZC(i,t+1) = ((Fnet_ZC(i,t)/Mass_mP(i))*timestep)+wvel_ZC(i,t);
    
            Dist1_ZC(i,t) = wvel_ZC(i,t) * timestep;
		    DistTot_ZC(i) = DistTot_ZC(i) + Dist1_ZC(i,t);
		    Acc_ZC(i,t) = (wvel_ZC(i, t+1) - wvel_ZC(i,t))/timestep;
		    
            if (Acc_ZC(i,t)< 0.001)
                FinalTime_ZC(i) = (t+1)*timestep;
                FinalStep_ZC(i) = t+1;
                wtFinal_ZC(i)=wvel_ZC(i, t+1);
                Dist1_ZC(i, t+1) = wvel_ZC(i, t+1) * timestep;
                DistTot_ZC(i) = DistTot_ZC(i) + Dist1_ZC(i, t+1);
                ReFinal_ZC(i) = abs((rho_p(i) * wvel_ZC(i, t+1) * d_equi(i))/ vis_dyn(i));
		        CdFinal_ZC(i) = (58.58*(shape_ASF(i)^0.1936))/(ReFinal_ZC(i)^0.8273);
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

Results_ZC_Proj = zeros(1, 4);
Results_ZC_Proj(1) = T_1000Times;
Results_ZC_Proj(2) = mean(T_CPU);
Results_ZC_Proj(3) = tEnd;
Results_ZC_Proj(4) = mean(tEnd_Loop);

Results_ZC_Proj = transpose(Results_ZC_Proj);

%% Combine results arrays:
All_Results_ZC = [Results_ZC_SA Results_ZC_Proj];

%%
Table_ZC_Timing = array2table(All_Results_ZC, "RowNames", ...
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
% Test 3: After Test 2, uncomment the line below for the file endingin '2'
% and run the programme.

% Ensure that when running the test, no other scripts are open on Matlab.

% % Test 1
% writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing1.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing1.xls', 'WriteRowNames', true);

% % Test 2
% writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing2.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing2.xls', 'WriteRowNames', true);

% Test 3
writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing3.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_Timing, './DragModelsTest/Output/20230215/ZC_Timing3.xls', 'WriteRowNames', true);

