%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: BagheriScript: VM
% Date created: 23.04.22
% Date last mostified: 22.07.22
% Purpose: To test the implementation of the Bagheri drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
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
Cd_meas = table2array(VM_Dataset(:, "CdMeasured"));

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

%% Bagheri' method 1
% <<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

timestep = 0.0002;

Re_BB = zeros(140, 10000);
ReFinal_BB = zeros(140, 1);
CdFinal_BB = zeros(140, 1);
FormFactor_S = zeros(140, 1);
FormFactor_N = zeros(140, 1);
Correction_S = zeros(140, 1);
Correction_N = zeros(140, 1);
ratio_density = zeros(140,1);
alpha2= zeros(140,1);
beta2 = zeros(140,1);
wvel_BB = zeros(140, 10000);
Cd_BB = zeros(140, 10000);
Fd_BB = zeros(140, 10000);
Fg_BB = zeros(140, 10000);
Fb_BB = zeros(140, 10000);
Fnet_BB = zeros(140, 10000);
Dist1_BB = zeros(140, 10000);
DistTot_BB = zeros(140, 1);
Acc_BB = zeros(140, 10000);
FinalTime_BB = zeros(140, 1);
FinalStep_BB = zeros(140, 1);
wtFinal_BB = zeros(140, 1);

% Set initial velocity
wvel_BB(:, 1) = 0.0001;
for i=1:140	
    
    FormFactor_S(i) = shape_flt(i)*(shape_eln(i)^1.3)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
	FormFactor_N(i) = (shape_flt(i)^2.0)*shape_eln(i)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
		
	Correction_S(i) = 0.5*((FormFactor_S(i)^(1.0/3.0))+(FormFactor_S(i)^(-1.0/3.0)));
		
	ratio_density(i) = rho_p(i) / rho_f(i);
	alpha2(i) = 0.45 + (10.0/((exp(2.5*log10(ratio_density(i))))+30));		
    beta2(i) = 1.0 -  (37.0/((exp(3.0*log10(ratio_density(i))))+100));	
	Correction_N(i) = 10.0^(alpha2(i)*((-1.0*log10(FormFactor_N(i)))^beta2(i)));
    
    for t=1:10000
		
		Re_BB(i, t) = abs((rho_p(i) * wvel_BB(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_BB(i,t) = ((24.0*Correction_S(i))/Re_BB(i,t))*(1+ 0.125*(((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))^(2.0/3.0))) ...
			+ (0.46*Correction_N(i))/(1 + (5330/((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))));
	
		Fd_BB(i,t) = 0.5*rho_f(i)*SA_mP(i)*(abs(wvel_BB(i,t))*wvel_BB(i,t))*Cd_BB(i,t);
	
		Fg_BB(i,t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_BB(i,t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_BB(i,t) = Fg_BB(i,t) - Fb_BB(i,t) - Fd_BB(i,t);
	
		wvel_BB(i,t+1) = ((Fnet_BB(i,t)/Mass_mP(i))*timestep)+wvel_BB(i,t);

        Dist1_BB(i,t) = wvel_BB(i,t) * timestep;
		DistTot_BB(i) = DistTot_BB(i) + Dist1_BB(i,t);
		Acc_BB(i,t) = (wvel_BB(i, t+1) - wvel_BB(i,t))/timestep;
		
        if (Acc_BB(i,t)< 0.001)
            FinalTime_BB(i) = (t+1)*timestep;
            FinalStep_BB(i) = t+1;
            wtFinal_BB(i)=wvel_BB(i, t+1);
            ReFinal_BB(i) = abs((rho_p(i) * wvel_BB(i, t+1) * d_equi(i))/ vis_dyn(i));
            CdFinal_BB(i) = ((24.0*Correction_S(i))/ReFinal_BB(i))*(1+ 0.125*(((ReFinal_BB(i)*Correction_N(i))/(Correction_S(i)))*(2.0/3.0))) ...
			                 + (0.46*Correction_N(i))/(1 + (5330/((ReFinal_BB(i)*Correction_N(i))/(Correction_S(i)))));
			break
        end
    end
end

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_BB(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
    hold on
    plot(timesec(1:(FinalStep_BB(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_BB = zeros(140, 11);

for i=1:140
    Results_BB(i, 1) = d_equi(i);
    Results_BB(i, 2) = CSF(i);
    Results_BB(i, 3) = wtFinal_BB(i);
    Results_BB(i, 4) = wvel_meas(i);
    Results_BB(i, 5) = FinalTime_BB(i);
    Results_BB(i, 6) = DistTot_BB(i);
    Results_BB(i, 7) = timestep;
    Results_BB(i, 8) = Reynolds(i);
    Results_BB(i, 9) = ReFinal_BB(i);
    Results_BB(i, 10) = Cd_meas(i);
    Results_BB(i, 11) = CdFinal_BB(i);
end 

Table_BB_SA = array2table(Results_BB, "VariableNames", ...
    {'ESD', 'CSF', 'Wt_Calc','Wt_Meas', 'Time', ...
    'Distance', 'timestep', 'Re_Meas', 'Re_Calc'...
    'Cd_Meas', 'Cd_Calc'});


Table_BB_SA = [VM_Dataset.Shape Table_BB_SA];
Table_BB_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_BB_SA, './DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_SA, './DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_SA.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error(i) = ((residual(i)/wvel_meas(i))*100);
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i)^2);
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE_SA = AE_Sum/140;
Abs_AE_SA = Abs_AE_Sum/140;
RMSE_SA = sqrt(RMSE_Sum/140);

% B) Fragments
residual_F3= zeros(80, 1);
Percentage_Error_F3 = zeros(80, 1);
AE_Sum_F3 = 0.0;
Abs_AE_Sum_F3 = 0.0;
Percentage_Error_sq_F3= zeros(80, 1);
RMSE_Sum_F3 = 0.0;

for i=1:80
    residual_F3(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i) / wvel_meas(i))*100);
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Abs_AE_Sum_F3 = Abs_AE_Sum_F3 + abs(Percentage_Error_F3(i));
    Percentage_Error_sq_F3(i) = (Percentage_Error_F3(i)^2);
    RMSE_Sum_F3 = RMSE_Sum_F3 + Percentage_Error_sq_F3(i);
end

AE_SA_F3 = AE_Sum_F3/80;
Abs_AE_SA_F3 = Abs_AE_Sum_F3/80;
RMSE_SA_F3 = sqrt(RMSE_Sum_F3/80);

% C) Fibres 
residual_F2= zeros(20, 1);
Percentage_Error_F2 = zeros(20, 1);
AE_Sum_F2 = 0.0;
Abs_AE_Sum_F2 = 0.0;
Percentage_Error_sq_F2= zeros(20, 1);
RMSE_Sum_F2 = 0.0;

for i=81:100
    residual_F2(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i) / wvel_meas(i))*100);
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i));
    Percentage_Error_sq_F2(i) = (Percentage_Error_F2(i)^2);
    RMSE_Sum_F2 = RMSE_Sum_F2 + Percentage_Error_sq_F2(i);
end

AE_SA_F2 = AE_Sum_F2/20;
Abs_AE_SA_F2 = Abs_AE_Sum_F2/20;
RMSE_SA_F2 = sqrt(RMSE_Sum_F2/20);

% D) Films
residual_F1= zeros(40, 1);
Percentage_Error_F1 = zeros(40, 1);
AE_Sum_F1 = 0.0;
Abs_AE_Sum_F1 = 0.0;
Percentage_Error_sq_F1= zeros(40, 1);
RMSE_Sum_F1 = 0.0;

for i=101:140
    residual_F1(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F1(i) = ((residual_F1(i) / wvel_meas(i))*100);
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i);
    Abs_AE_Sum_F1 = Abs_AE_Sum_F1 + abs(Percentage_Error_F1(i));
    Percentage_Error_sq_F1(i) = (Percentage_Error_F1(i)^2);
    RMSE_Sum_F1 = RMSE_Sum_F1 + Percentage_Error_sq_F1(i);
end

AE_SA_F1 = AE_Sum_F1/40;
Abs_AE_SA_F1 = Abs_AE_Sum_F1/40;
RMSE_SA_F1 = sqrt(RMSE_Sum_F1/40);

Error_table_shape = ["All"; "Fragment"; "Fibre"; "Film"];
Error_table_AE = [AE_SA; AE_SA_F3; AE_SA_F2; AE_SA_F1];
Error_table_Abs_AE = [Abs_AE_SA; Abs_AE_SA_F3; Abs_AE_SA_F2; Abs_AE_SA_F1];
Error_table_RMSE = [RMSE_SA; RMSE_SA_F3; RMSE_SA_F2; RMSE_SA_F1];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_Abs_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_SA.xls', 'WriteRowNames', true);
%% Bagheri Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective area
% in the calculation of the drag force.

% Set timestep
timestep = 0.0002;

% Set up variable arrays
Re_BB = zeros(140, 10000);
ReFinal_BB = zeros(140, 1);
CdFinal_BB = zeros(140, 1);
FormFactor_S = zeros(140, 1);
FormFactor_N = zeros(140, 1);
Correction_S = zeros(140, 1);
Correction_N = zeros(140, 1);
ratio_density = zeros(140,1);
alpha2= zeros(140,1);
beta2 = zeros(140,1);
wvel_BB = zeros(140, 10000);
Cd_BB = zeros(140, 10000);
Fd_BB = zeros(140, 10000);
Fg_BB = zeros(140, 10000);
Fb_BB = zeros(140, 10000);
Fnet_BB = zeros(140, 10000);
Dist1_BB = zeros(140, 10000);
DistTot_BB = zeros(140, 1);
Acc_BB = zeros(140, 10000);
FinalTime_BB = zeros(140, 1);
FinalStep_BB = zeros(140, 1);
wtFinal_BB = zeros(140, 1);

% Set initial velocity
wvel_BB(:, 1) = 0.0001;
for i=1:140	
    
    FormFactor_S(i) = shape_flt(i)*(shape_eln(i)^1.3)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
	FormFactor_N(i) = (shape_flt(i)^2.0)*shape_eln(i)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
		
	Correction_S(i) = 0.5*((FormFactor_S(i)^(1.0/3.0))+(FormFactor_S(i)^(-1.0/3.0)));
		
	ratio_density(i) = rho_p(i) / rho_f(i);
	alpha2(i) = 0.45 + (10.0/((exp(2.5*log10(ratio_density(i))))+30));		
    beta2(i) = 1.0 -  (37.0/((exp(3.0*log10(ratio_density(i))))+100));	
	Correction_N(i) = 10.0^(alpha2(i)*((-1.0*log10(FormFactor_N(i)))^beta2(i)));
    
    for t=1:10000
		
		Re_BB(i, t) = abs((rho_p(i) * wvel_BB(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_BB(i,t) = ((24.0*Correction_S(i))/Re_BB(i,t))*(1+ 0.125*(((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))^(2.0/3.0))) ...
			+ (0.46*Correction_N(i))/(1 + (5330/((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))));
	
		Fd_BB(i,t) = 0.5*rho_f(i)*ProjA_ESD(i)*(abs(wvel_BB(i,t))*wvel_BB(i,t))*Cd_BB(i,t);
	
		Fg_BB(i,t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_BB(i,t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_BB(i,t) = Fg_BB(i,t) - Fb_BB(i,t) - Fd_BB(i,t);
	
		wvel_BB(i,t+1) = ((Fnet_BB(i,t)/Mass_mP(i))*timestep)+wvel_BB(i,t);

        Dist1_BB(i,t) = wvel_BB(i,t) * timestep;
		DistTot_BB(i) = DistTot_BB(i) + Dist1_BB(i,t);
		Acc_BB(i,t) = (wvel_BB(i, t+1) - wvel_BB(i,t))/timestep;
		
        if (Acc_BB(i,t)< 0.001)
            FinalTime_BB(i) = (t+1)*timestep;
            FinalStep_BB(i) = t+1;
            wtFinal_BB(i)=wvel_BB(i, t+1);
            ReFinal_BB(i) = abs((rho_p(i) * wvel_BB(i, t+1) * d_equi(i))/ vis_dyn(i));
            CdFinal_BB(i) = ((24.0*Correction_S(i))/ReFinal_BB(i))*(1+ 0.125*(((ReFinal_BB(i)*Correction_N(i))/(Correction_S(i)))*(2.0/3.0))) ...
			+ (0.46*Correction_N(i))/(1 + (5330/((ReFinal_BB(i)*Correction_N(i))/(Correction_S(i)))));
			break
        end
    end
end

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_BB(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
    hold on
    plot(timesec(1:(FinalStep_BB(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_BB = zeros(140, 11);

for i=1:140
    Results_BB(i, 1) = d_equi(i);
    Results_BB(i, 2) = CSF(i);
    Results_BB(i, 3) = wtFinal_BB(i);
    Results_BB(i, 4) = wvel_meas(i);
    Results_BB(i, 5) = FinalTime_BB(i);
    Results_BB(i, 6) = DistTot_BB(i);
    Results_BB(i, 7) = timestep;
    Results_BB(i, 8) = Reynolds(i);
    Results_BB(i, 9) = ReFinal_BB(i);
    Results_BB(i, 10) = Cd_meas(i);
    Results_BB(i, 11) = CdFinal_BB(i);
end 

Table_BB_Proj = array2table(Results_BB, "VariableNames", ...
    {'ESD', 'CSF', 'Wt_Calc','Wt_Meas', 'Time', ...
    'Distance', 'timestep', 'Re_Meas', 'Re_Calc', ...
    'Cd_Meas', 'Cd_Calc'});

Table_BB_Proj = [VM_Dataset.Shape Table_BB_Proj];
Table_BB_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_Proj.xls', 'WriteRowNames', true);

%% Distance assumption calculation and plot

DistConst_BB = zeros(140, 1);

for i=1:140
    DistConst_BB(i) = FinalTime_BB(i) * wtFinal_BB(i);
end

% Fit linear model through the intercept: SA
lm_BBDist = fitlm(DistTot_BB, DistConst_BB, 'y~-1+x1');
m_BBDist = lm_BBDist.Coefficients.Estimate(1);
fitY_BBDist = zeros(1000, 1);
% Generate data using linear model:
n1=[max(DistTot_BB), max(DistConst_BB)] ;
nMax = max(n1);
nVal=linspace(0.00001, nMax, 1000);
r_sq_Dist = lm_BBDist.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBDist(i) = m_BBDist * nVal(i);
end

subplot(1, 2, 2)
plot(DistTot_BB, DistConst_BB, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[1, 1, 0]')
ylabel('Distance travelled at constant velocity (m)')
xlabel('Distance travelled in attaining terminal velocity (m)')
title('Bagheri and Bonadonna (2016): Using Particle Projected Area.')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBDist, '--k')
plot(nVal, 0.7*nVal, ':k')
plot(nVal, 1.3*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBDist, r_sq_Dist), '', '', 'location', 'best');
set(gca,'YLim', [0.00005, nMax*1.3] )
set(gca,'XLim', [0.00005, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'Maximized')
exportgraphics(gcf, './DragModelsTest/Output/20220621/Distance/BB_DistanceProj.jpg', 'Resolution', 300)

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error(i) = (residual(i)/wvel_meas(i))*100;
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i))^2;
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE_Proj = AE_Sum/140;
Abs_AE_Proj = Abs_AE_Sum/140;
RMSE_Proj = sqrt(RMSE_Sum/140);

% B) Fragments
residual_F3= zeros(80, 1);
Percentage_Error_F3 = zeros(80, 1);
AE_Sum_F3 = 0.0;
Abs_AE_Sum_F3 = 0.0;
Percentage_Error_sq_F3= zeros(80, 1);
RMSE_Sum_F3 = 0.0;

for i=1:80
    residual_F3(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i) / wvel_meas(i))*100);
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Abs_AE_Sum_F3 = Abs_AE_Sum_F3 + abs(Percentage_Error_F3(i));
    Percentage_Error_sq_F3(i) = (Percentage_Error_F3(i)^2);
    RMSE_Sum_F3 = RMSE_Sum_F3 + Percentage_Error_sq_F3(i);
end

AE_Proj_F3 = AE_Sum_F3/80;
Abs_AE_Proj_F3 = Abs_AE_Sum_F3/80;
RMSE_Proj_F3 = sqrt(RMSE_Sum_F3/80);

% C) Fibres 
residual_F2= zeros(20, 1);
Percentage_Error_F2 = zeros(20, 1);
AE_Sum_F2 = 0.0;
Abs_AE_Sum_F2 = 0.0;
Percentage_Error_sq_F2= zeros(20, 1);
RMSE_Sum_F2 = 0.0;

for i=81:100
    residual_F2(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i) / wvel_meas(i))*100);
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i));
    Percentage_Error_sq_F2(i) = (Percentage_Error_F2(i)^2);
    RMSE_Sum_F2 = RMSE_Sum_F2 + Percentage_Error_sq_F2(i);
end

AE_Proj_F2 = AE_Sum_F2/20;
Abs_AE_Proj_F2 = Abs_AE_Sum_F2/20;
RMSE_Proj_F2 = sqrt(RMSE_Sum_F2/20);

% D) Films
residual_F1= zeros(40, 1);
Percentage_Error_F1 = zeros(40, 1);
AE_Sum_F1 = 0.0;
Abs_AE_Sum_F1 = 0.0;
Percentage_Error_sq_F1= zeros(40, 1);
RMSE_Sum_F1 = 0.0;

for i=101:140
    residual_F1(i) = (wtFinal_BB(i) - wvel_meas(i));
    Percentage_Error_F1(i) = ((residual_F1(i) / wvel_meas(i))*100);
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i);
    Abs_AE_Sum_F1 = Abs_AE_Sum_F1 + abs(Percentage_Error_F1(i));
    Percentage_Error_sq_F1(i) = (Percentage_Error_F1(i))^2;
    RMSE_Sum_F1 = RMSE_Sum_F1 + Percentage_Error_sq_F1(i);
end

AE_Proj_F1 = AE_Sum_F1/40;
Abs_AE_Proj_F1 = Abs_AE_Sum_F1/40;
RMSE_Proj_F1 = sqrt(RMSE_Sum_F1/40);

Error_table_shape = ["All"; "Fragment"; "Fibre"; "Film"];
Error_table_AE = [AE_Proj; AE_Proj_F3; AE_Proj_F2; AE_Proj_F1];
Error_table_Abs_AE = [Abs_AE_Proj; Abs_AE_Proj_F3; Abs_AE_Proj_F2; Abs_AE_Proj_F1];
Error_table_RMSE = [RMSE_Proj; RMSE_Proj_F3; RMSE_Proj_F2; RMSE_Proj_F1];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_Abs_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Bagheri output
% <<<<<<<<<<<<<<<<<<<
clear
Table_BB_SA= readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_SA.txt", "Delimiter", ",");
Table_BB_Proj= readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_Proj.txt", "Delimiter", ",");

%% A1) wt against ESD (ALL)
% ==========================

% Method 1: All
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

% Method 2: All
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_SA{1:80, "ESD"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "ESD"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "ESD"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_Proj{1:80, "ESD"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "ESD"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "ESD"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF (ALL)
% ==========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_BB_SA.('CSF'), Table_BB_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_SA.('CSF'), Table_BB_SA.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('CSF'), Table_BB_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_SA{1:80, "CSF"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "CSF"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "CSF"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_BB_Proj{1:80, "CSF"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "CSF"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "CSF"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C2) wt against wt measured using Matlab fitlm function
% ========================================================

% Fit linear model through the intercept: SA
lm_BBSA = fitlm(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'y~-1+x1');
m_BBSA = lm_BBSA.Coefficients.Estimate(1);
fitY_BBSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_BB_SA.Wt_Calc), max(Table_BB_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_BBSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBSA(i) = m_BBSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Particle Surface Area.')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBSA, '--k')
plot(nVal, 0.7*nVal, ':k')
plot(nVal, 1.3*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSA, r_sq_SA), '', '', 'location', 'best');
set(gca,'YLim', [0.0001, nMax*1.3] )
set(gca,'XLim', [0.0001, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_BBProj = fitlm(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_BB_Proj.Wt_Calc), max(Table_BB_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_Proj = lm_BBProj.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Projected Area of Volume Equivalent Sphere.')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBProj, '--k')
plot(nVal, 0.7*nVal, ':k')
plot(nVal, 1.3*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProj, r_sq_Proj), 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_MeasVsCalc_Fit.jpg', 'Resolution', 300);

%% C B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_BBSA = fitlm(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'y~-1+x1');
m_BBSA = lm_BBSA.Coefficients.Estimate(1);
fitY_BBSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_BB_SA.Wt_Calc), max(Table_BB_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_BBSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBSA(i) = m_BBSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_BB_SA{1:80, "Wt_Meas"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using particle Surface Area.')
hold on
plot(Table_BB_SA{81:100, "Wt_Meas"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "Wt_Meas"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBSA, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSA, r_sq_SA), 'location', 'best');
set(gca,'YLim', [0.0001, nMax*1.3] )
set(gca,'XLim', [0.0001, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_BBProj = fitlm(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_BB_Proj.Wt_Calc), max(Table_BB_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_Proj = lm_BBProj.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_BB_Proj{1:80, "Wt_Meas"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Projected Area of Volume Equivalent Sphere.')
hold on
plot(Table_BB_Proj{81:100, "Wt_Meas"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "Wt_Meas"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBProj, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProj, r_sq_Proj), 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 300);

%% C C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_BBSAF3 = fitlm(Table_BB_SA{1:80, "Wt_Meas"}, Table_BB_SA{1:80, "Wt_Calc"}, 'y~-1+x1');
m_BBSAF3 = lm_BBSAF3.Coefficients.Estimate(1);
fitY_BBSAF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_BB_SA{1:80, "Wt_Calc"}), max(Table_BB_SA{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_SAF3 = lm_BBSAF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBSAF3(i) = m_BBSAF3 * nVal_F3(i);
end

subplot(1, 2, 1)
plot(Table_BB_SA{1:80, "Wt_Meas"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Particle Surface Area.')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_BBSAF3, '--b')
plot(nVal_F3, 1.3*nVal_F3, ':k')
plot(nVal_F3, 0.7*nVal_F3, ':k')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSAF3, r_sq_SAF3), 'location', 'best');
set(gca,'YLim', [0.0001, nMax_F3*1.1] )
set(gca,'XLim', [0.0001, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_BBProjF3 = fitlm(Table_BB_Proj{1:80, "Wt_Meas"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'y~-1+x1');
m_BBProjF3 = lm_BBProjF3.Coefficients.Estimate(1);
fitY_BBProjF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_BB_Proj{1:80, "Wt_Calc"}), max(Table_BB_Proj{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_ProjF3 = lm_BBProjF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBProjF3(i) = m_BBProjF3 * nVal_F3(i);
end

subplot(1, 2, 2)
plot(Table_BB_Proj{1:80, "Wt_Meas"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Projected Area of Volume Equivalent Sphere.')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_BBProjF3, '--b')
plot(nVal_F3, 1.3*nVal_F3, ':k')
plot(nVal_F3, 0.7*nVal_F3, ':k')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProjF3, r_sq_ProjF3), 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_MeasVsCalc_FitF3.jpg', 'Resolution', 300);

%% C D) Plot fibres separately with fitted model

% Fit linear model through the intercept: SA
lm_BBSAF2 = fitlm(Table_BB_SA{81:100, "Wt_Meas"}, Table_BB_SA{81:100, "Wt_Calc"}, 'y~-1+x1');
m_BBSAF2 = lm_BBSAF2.Coefficients.Estimate(1);
fitY_BBSAF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_BB_SA{81:100, "Wt_Calc"}), max(Table_BB_SA{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_SAF2 = lm_BBSAF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBSAF2(i) = m_BBSAF2 * nVal_F2(i);
end

subplot(1, 2, 1)
plot(Table_BB_SA{81:100, "Wt_Meas"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Particle Surface Area.')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_BBSAF2, '--r')
plot(nVal_F2, 1.3*nVal_F2, ':k')
plot(nVal_F2, 0.7*nVal_F2, ':k')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSAF2, r_sq_SAF2), '', '',  'location', 'best');
set(gca,'YLim', [0.001, nMax_F2*1.1] )
set(gca,'XLim', [0.001, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_BBProjF2 = fitlm(Table_BB_Proj{81:100, "Wt_Meas"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'y~-1+x1');
m_BBProjF2 = lm_BBProjF2.Coefficients.Estimate(1);
fitY_BBProjF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_BB_Proj{81:100, "Wt_Calc"}), max(Table_BB_Proj{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_ProjF2 = lm_BBProjF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBProjF2(i) = m_BBProjF2 * nVal_F2(i);
end

subplot(1, 2, 2)
plot(Table_BB_Proj{81:100, "Wt_Meas"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Projected Area of Volume Equivalent Sphere.')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_BBProjF2, '--r')
plot(nVal_F2, 1.3*nVal_F2, ':k')
plot(nVal_F2, 0.7*nVal_F2, ':k')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProjF2, r_sq_ProjF2), '', '', 'location', 'best');
set(gca,'YLim', [0.007, nMax_F2*1.1] )
set(gca,'XLim', [0.007, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_MeasVsCalc_FitF2.jpg', 'Resolution', 300);

%% C E) Plot film separately with fitted model

% Fit linear model through the intercept: SA
lm_BBSAF1 = fitlm(Table_BB_SA{101:140, "Wt_Meas"}, Table_BB_SA{101:140, "Wt_Calc"}, 'y~-1+x1');
m_BBSAF1 = lm_BBSAF1.Coefficients.Estimate(1);
fitY_BBSAF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_BB_SA{101:140, "Wt_Calc"}), max(Table_BB_SA{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_SAF1 = lm_BBSAF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBSAF1(i) = m_BBSAF1 * nVal_F1(i);
end

subplot(1, 2, 1)
plot(Table_BB_SA{101:140, "Wt_Meas"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Particle Surface Area.')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_BBSAF1, '--g')
plot(nVal_F1, 1.3*nVal_F1, ':k')
plot(nVal_F1, 0.7*nVal_F1, ':k')
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSAF1, r_sq_SAF1), '', '', 'location', 'best');
set(gca,'YLim', [0.0001, nMax_F1*1.1] )
set(gca,'XLim', [0.0001, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_BBProjF1 = fitlm(Table_BB_Proj{101:140, "Wt_Meas"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'y~-1+x1');
m_BBProjF1 = lm_BBProjF1.Coefficients.Estimate(1);
fitY_BBProjF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_BB_Proj{101:140, "Wt_Calc"}), max(Table_BB_Proj{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.001, nMax_F1, 1000);
r_sq_ProjF1 = lm_BBProjF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_BBProjF1(i) = m_BBProjF1 * nVal_F1(i);
end

subplot(1, 2, 2)
plot(Table_BB_Proj{101:140, "Wt_Meas"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Bagheri and Bonadonna (2016): Using Projected Area of Volume Equivalent Sphere.')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_BBProjF1, '--g')
plot(nVal_F1, 1.3*nVal_F1, ':k')
plot(nVal_F1, 0.7*nVal_F1, ':k')
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProjF1, r_sq_ProjF1), 'location', 'best');
set(gca,'YLim', [0.003, nMax_F1*1.1] )
set(gca,'XLim', [0.003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_MeasVsCalc_FitF1.jpg', 'Resolution', 300);

%% Combine all m and r_sq values into the error table: Projected Area
Error_table = readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_Proj.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

BB_rsq_proj = [r_sq_Proj; r_sq_ProjF3; r_sq_ProjF2; r_sq_ProjF1];
BB_m_proj = [m_BBProj; m_BBProjF3; m_BBProjF2; m_BBProjF1];

BB_Proj_Table = array2table([BB_m_proj BB_rsq_proj]);
BB_Proj_Table.Properties.VariableNames = Col_names;
BB_Proj_Table.Properties.RowNames = Row_names;

Error_table_Proj = [Error_table BB_Proj_Table];

writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Bagheri/BagheriFinalTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Bagheri/BagheriFinalTableVM_Proj.xls', 'WriteRowNames', true);

%% Combine all m and r_sq values into the error table: Surface Area
Error_table = readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriErrorTableVM_SA.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

BB_rsq_SA = [r_sq_SA; r_sq_SAF3; r_sq_SAF2; r_sq_SAF1];
BB_m_SA = [m_BBSA; m_BBSAF3; m_BBSAF2; m_BBSAF1];

BB_SA_Table = array2table([BB_m_SA BB_rsq_SA]);
BB_SA_Table.Properties.VariableNames = Col_names;
BB_SA_Table.Properties.RowNames = Row_names;

Error_table_SA = [Error_table BB_SA_Table];

writetable(Error_table_SA, './DragModelsTest/Output/20220621/Bagheri/BagheriFinalTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_SA, './DragModelsTest/Output/20220621/Bagheri/BagheriFinalTableVM_SA.xls', 'WriteRowNames', true);
%% D) Re against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_BB_SA.('Re_Meas'), Table_BB_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_SA.('Re_Calc'), Table_BB_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Bagheri Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_BB_Proj.('Re_Meas'), Table_BB_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_Proj.('Re_Calc'), Table_BB_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ReVsCd.jpg', 'Resolution', 300);

%% D2) Re against Cd (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('Re_Meas'), Table_BB_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_SA{1:80, "Re_Calc"}, Table_BB_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "Re_Calc"}, Table_BB_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "Re_Calc"}, Table_BB_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model. Using Particle Surface Area.')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_BB_Proj.('Re_Meas'), Table_BB_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_Proj{1:80, "Re_Calc"}, Table_BB_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "Re_Calc"}, Table_BB_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "Re_Calc"}, Table_BB_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ReVsCd_Shapes.jpg', 'Resolution', 300);

%% E1) ESD against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Bagheri Model: Using Particle Surface Area.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ESDVsCd.jpg', 'Resolution', 300);

%% E2) ESD against Cd (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_SA{1:80, "ESD"}, Table_BB_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "ESD"}, Table_BB_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "ESD"}, Table_BB_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Particle Surface Area.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_BB_Proj{1:80, "ESD"}, Table_BB_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "ESD"}, Table_BB_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "ESD"}, Table_BB_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model: Using Projected Area of Volume Equivalent Sphere.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Bagheri/BagheriVM_ESDVsCd_Shapes.jpg', 'Resolution', 300);