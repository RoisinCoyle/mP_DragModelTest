%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DioguardiScript: VM
% Date created: 23.04.22
% Date last mostified: 02.03.23
% Purpose: To test the implementation of the Dioguardi drag model on a range of
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

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_Dio(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
    hold on
    plot(timesec(1:(FinalStep_Dio(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_Dio = zeros(140, 11);

for i=1:140
    Results_Dio(i, 1) = d_equi(i);
    Results_Dio(i, 2) = CSF(i);
    Results_Dio(i, 3) = wtFinal_Dio(i);
    Results_Dio(i, 4) = wvel_meas(i);
    Results_Dio(i, 5) = FinalTime_Dio(i);
    Results_Dio(i, 6) = DistTot_Dio(i);
    Results_Dio(i, 7) = timestep;
    Results_Dio(i, 8) = Reynolds(i);
    Results_Dio(i, 9) = ReFinal_Dio(i);
    Results_Dio(i, 10) = Cd_meas(i);
    Results_Dio(i, 11) = CdFinal_Dio(i);
end 

Table_Dio_SA = array2table(Results_Dio, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});


Table_Dio_SA = [VM_Dataset.Shape Table_Dio_SA];
Table_Dio_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Dio_SA, './DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dio_SA, './DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_SA.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error(i) = ((residual(i)) / wvel_meas(i))*100;
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i))^2;
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
    residual_F3(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i)) / wvel_meas(i))*100;
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Abs_AE_Sum_F3 = Abs_AE_Sum_F3 + abs(Percentage_Error_F3(i));
    Percentage_Error_sq_F3(i) = (Percentage_Error_F3(i))^2;
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
    residual_F2(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i)) / wvel_meas(i))*100;
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i));
    Percentage_Error_sq_F2(i) = (Percentage_Error_F2(i))^2;
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
    residual_F1(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F1(i) = ((residual_F1(i)) / wvel_meas(i))*100;
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i);
    Abs_AE_Sum_F1 = Abs_AE_Sum_F1 + abs(Percentage_Error_F1(i));
    Percentage_Error_sq_F1(i) = (Percentage_Error_F1(i))^2;
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

writetable(Error_table, './DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_SA.xls', 'WriteRowNames', true);

%% Dioguardi Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area in the calculation of the drag force.

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

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_Dio(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
    hold on
    plot(timesec(1:(FinalStep_Dio(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_Dio = zeros(140, 11);

for i=1:140
    Results_Dio(i, 1) = d_equi(i);
    Results_Dio(i, 2) = CSF(i);
    Results_Dio(i, 3) = wtFinal_Dio(i);
    Results_Dio(i, 4) = wvel_meas(i);
    Results_Dio(i, 5) = FinalTime_Dio(i);
    Results_Dio(i, 6) = DistTot_Dio(i);
    Results_Dio(i, 7) = timestep;
    Results_Dio(i, 8) = Reynolds(i);
    Results_Dio(i, 9) = ReFinal_Dio(i);
    Results_Dio(i, 10) = Cd_meas(i);
    Results_Dio(i, 11) = CdFinal_Dio(i);
end 

Table_Dio_Proj = array2table(Results_Dio, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Results_Dioguardi = zeros(140, 5);

Table_Dio_Proj = [VM_Dataset.Shape Table_Dio_Proj];
Table_Dio_Proj.Properties.VariableNames(1) = {'Shape'};

% writetable(Table_Dio_Proj, './DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dio_Proj, './DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_Proj.xls', 'WriteRowNames', true);

%% Distance assumption calculation and plot

DistConst_Dio = zeros(140, 1);

for i=1:140
    DistConst_Dio(i) = FinalTime_Dio(i) * wtFinal_Dio(i);
end

% Fit linear model through the intercept: SA
lm_DioDist = fitlm(DistTot_Dio, DistConst_Dio, 'y~-1+x1');
m_DioDist = lm_DioDist.Coefficients.Estimate(1);
fitY_DioDist = zeros(1000, 1);
% Generate data using linear model:
n1=[max(DistTot_Dio), max(DistConst_Dio)] ;
nMax = max(n1);

nVal=linspace(0.00001, nMax, 1000);
r_sq_Dist = lm_DioDist.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioDist(i) = m_DioDist * nVal(i);
end

subplot(1, 2, 2)
plot(DistTot_Dio, DistConst_Dio, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[1, 1, 0]')
ylabel('Distance travelled at constant velocity (m)')
xlabel('Distance travelled in attaining terminal velocity (m)')
title(sprintf('Figure showing that distance travelled in attaining terminal settling velocity is approximately \n\r equal to the distance travelled in the same time interval at constant velocity'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using particle projection area as effective area.'))
hold on
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_DioDist, '--k', 'LineWidth', 1)
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 2, 'Color', [.7 .7 .7])
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 2, 'Color', [.7 .7 .7])
legend(' ', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioDist, r_sq_Dist), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.00001, nMax*1.3] )
set(gca,'XLim', [0.00001, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'Maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230301/Distance/Dio_DistanceProj.jpg', 'Resolution', 1200)

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error(i) = ((residual(i)) / wvel_meas(i))*100;
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
    residual_F3(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i)) / wvel_meas(i))*100;
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Abs_AE_Sum_F3 = Abs_AE_Sum_F3 + abs(Percentage_Error_F3(i));
    Percentage_Error_sq_F3(i) = (Percentage_Error_F3(i))^2;
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
    residual_F2(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i)) / wvel_meas(i))*100;
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i));
    Percentage_Error_sq_F2(i) = (Percentage_Error_F2(i))^2;
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
    residual_F1(i) = (wtFinal_Dio(i) - wvel_meas(i));
    Percentage_Error_F1(i) = ((residual_F1(i)) / wvel_meas(i))*100;
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

writetable(Error_table, './DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Dioguardi output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Dio_SA= readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_SA.txt", "Delimiter", ",");
Table_Dio_Proj= readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_Proj.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

% Method 1: All
subplot(1, 2, 1)
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dioguardi Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

% Method 2: All
subplot(1, 2, 2)
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dioguardi Model: Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_SA{1:80, "ESD"}, Table_Dio_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_SA{81:100, "ESD"}, Table_Dio_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_SA{101:140, "ESD"}, Table_Dio_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_Proj{1:80, "ESD"}, Table_Dio_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_Proj{81:100, "ESD"}, Table_Dio_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "ESD"}, Table_Dio_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model: Using Projected Area of Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Dio_SA.('CSF'), Table_Dio_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_SA.('CSF'), Table_Dio_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dioguardi Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Dio_Proj.('CSF'), Table_Dio_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_Proj.('CSF'), Table_Dio_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dioguardi Model: Using Projected Area of Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Dio_SA.('CSF'), Table_Dio_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_SA{1:80, "CSF"}, Table_Dio_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_SA{81:100, "CSF"}, Table_Dio_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_SA{101:140, "CSF"}, Table_Dio_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model: Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Dio_Proj.('CSF'), Table_Dio_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dio_Proj{1:80, "CSF"}, Table_Dio_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_Proj{81:100, "CSF"}, Table_Dio_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "CSF"}, Table_Dio_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model: Using Projected Area of Equivalent Sphere.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_CSFVsW_Shapes.jpg', 'Resolution', 300);


%% C) wt against wt measured using Matlab fitlm function
% ========================================================

% C A) All shapes

% Fit linear model through the intercept: SA
lm_DioSA = fitlm(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'y~-1+x1');
m_DioSA = lm_DioSA.Coefficients.Estimate(1);
fitY_DioSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Dio_SA.Wt), max(Table_Dio_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_DioSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioSA(i) = m_DioSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r particle surface area as effective area.'))
hold on
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_DioSA, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSA, r_sq_SA), 'y= x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.0003, nMax*1.1] )
set(gca,'XLim', [0.0003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_DioProj = fitlm(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Dio_Proj.Wt), max(Table_Dio_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.001, nMax, 1000);
r_sq_Proj = lm_DioProj.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Bagheri and Bonadonna (2016) using \n\r projection area of volume equivalent sphere as effective area.'))
hold on
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_DioProj, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProj, r_sq_Proj), 'y= x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.1] )
set(gca,'XLim', [0.003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Dioguardi/DioguardiVM_MeasVsCalc_Fit.jpg', 'Resolution', 1200);

%% C B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_DioSA = fitlm(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'y~-1+x1');
m_DioSA = lm_DioSA.Coefficients.Estimate(1);
fitY_DioSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Dio_SA.Wt), max(Table_Dio_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_DioSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioSA(i) = m_DioSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Dio_SA{1:80, "Wt_Meas"}, Table_Dio_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r particle surface area as effective area.'))
hold on
plot(Table_Dio_SA{81:100, "Wt_Meas"}, Table_Dio_SA{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_SA{101:140, "Wt_Meas"}, Table_Dio_SA{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_DioSA, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSA, r_sq_SA), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.0003, nMax*1.3] )
set(gca,'XLim', [0.0003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_DioProj = fitlm(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Dio_Proj.Wt), max(Table_Dio_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_Proj = lm_DioProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Dio_Proj{1:80, "Wt_Meas"}, Table_Dio_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r using projection area of volume equivalent sphere as effective area.'))
hold on
plot(Table_Dio_Proj{81:100, "Wt_Meas"}, Table_Dio_Proj{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "Wt_Meas"}, Table_Dio_Proj{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_DioProj, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProj, r_sq_Proj), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Dioguardi/DioVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 1200);

%% C C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_DioSAF3 = fitlm(Table_Dio_SA{1:80, "Wt_Meas"}, Table_Dio_SA{1:80, "Wt"}, 'y~-1+x1');
m_DioSAF3 = lm_DioSAF3.Coefficients.Estimate(1);
fitY_DioSAF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Dio_SA{1:80, "Wt"}), max(Table_Dio_SA{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0, nMax_F3, 1000);
r_sq_SAF3 = lm_DioSAF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioSAF3(i) = m_DioSAF3 * nVal_F3(i);
end

subplot(1, 2, 1)
plot(Table_Dio_SA{1:80, "Wt_Meas"}, Table_Dio_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fragment terminal settling velocity to \n\r mP fragment terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r particle surface area as effective area.'))
hold on
plot(nVal_F3, nVal_F3, '-k', 'LineWidth', 1)
plot(nVal_F3, fitY_DioSAF3, '--b', 'LineWidth', 1)
plot(nVal_F3, 1.3*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F3, 0.7*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSAF3, r_sq_SAF3), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.0003, nMax_F3*1.1] )
set(gca,'XLim', [0.0003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_DioProjF3 = fitlm(Table_Dio_Proj{1:80, "Wt_Meas"}, Table_Dio_Proj{1:80, "Wt"}, 'y~-1+x1');
m_DioProjF3 = lm_DioProjF3.Coefficients.Estimate(1);
fitY_DioProjF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Dio_Proj{1:80, "Wt"}), max(Table_Dio_Proj{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_ProjF3 = lm_DioProjF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioProjF3(i) = m_DioProjF3 * nVal_F3(i);
end

subplot(1, 2, 2)
plot(Table_Dio_Proj{1:80, "Wt_Meas"}, Table_Dio_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fragment terminal settling velocity to \n\r mP fragment terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r projection area of volume equivalent sphere as effective area.'))
hold on
plot(nVal_F3, nVal_F3, '-k', 'LineWidth', 1)
plot(nVal_F3, fitY_DioProjF3, '--b', 'LineWidth', 1)
plot(nVal_F3, 1.3*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F3, 0.7*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProjF3, r_sq_ProjF3), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Dioguardi/DioVM_MeasVsCalc_FitF3.jpg', 'Resolution', 1200);

%% C D) Plot fibres separately with fitted model

% Fit linear model through the intercept: SA
lm_DioSAF2 = fitlm(Table_Dio_SA{81:100, "Wt_Meas"}, Table_Dio_SA{81:100, "Wt"}, 'y~-1+x1');
m_DioSAF2 = lm_DioSAF2.Coefficients.Estimate(1);
fitY_DioSAF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Dio_SA{81:100, "Wt"}), max(Table_Dio_SA{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_SAF2 = lm_DioSAF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioSAF2(i) = m_DioSAF2 * nVal_F2(i);
end

subplot(1, 2, 1)
plot(Table_Dio_SA{81:100, "Wt_Meas"}, Table_Dio_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fibre terminal settling velocity to \n\r mP fibre terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r particle surface area as effective area.'))
hold on
plot(nVal_F2, nVal_F2, '-k', 'LineWidth', 1)
plot(nVal_F2, fitY_DioSAF2, '--r', 'LineWidth', 1)
plot(nVal_F2, 1.3*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F2, 0.7*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSAF2, r_sq_SAF2), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.0003, nMax_F2*1.1] )
set(gca,'XLim', [0.0003, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

hold off

% Fit linear model through the intercept: Projected area
lm_DioProjF2 = fitlm(Table_Dio_Proj{81:100, "Wt_Meas"}, Table_Dio_Proj{81:100, "Wt"}, 'y~-1+x1');
m_DioProjF2 = lm_DioProjF2.Coefficients.Estimate(1);
fitY_DioProjF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Dio_Proj{81:100, "Wt"}), max(Table_Dio_Proj{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_ProjF2 = lm_DioProjF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioProjF2(i) = m_DioProjF2 * nVal_F2(i);
end

subplot(1, 2, 2)
plot(Table_Dio_Proj{81:100, "Wt_Meas"}, Table_Dio_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fibre terminal settling velocity to \n\r mP fibre terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r projection area of volume equivalent sphere as effective area.'))
hold on
plot(nVal_F2, nVal_F2, '-k', 'LineWidth', 1)
plot(nVal_F2, fitY_DioProjF2, '--r', 'LineWidth', 1)
plot(nVal_F2, 1.3*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F2, 0.7*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProjF2, r_sq_ProjF2), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F2*1.1] )
set(gca,'XLim', [0.003, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Dioguardi/DioVM_MeasVsCalc_FitF2.jpg', 'Resolution', 1200);

%% C E) Plot film separately with fitted model

% Fit linear model through the intercept: SA
lm_DioSAF1 = fitlm(Table_Dio_SA{101:140, "Wt_Meas"}, Table_Dio_SA{101:140, "Wt"}, 'y~-1+x1');
m_DioSAF1 = lm_DioSAF1.Coefficients.Estimate(1);
fitY_DioSAF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Dio_SA{101:140, "Wt"}), max(Table_Dio_SA{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_SAF1 = lm_DioSAF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioSAF1(i) = m_DioSAF1 * nVal_F1(i);
end

subplot(1, 2, 1)
plot(Table_Dio_SA{101:140, "Wt_Meas"}, Table_Dio_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP film terminal settling velocity to \n\r mP film terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r particle surface area as effective area.'))
hold on
plot(nVal_F1, nVal_F1, '-k', 'LineWidth', 1)
plot(nVal_F1, fitY_DioSAF1, '--g', 'LineWidth', 1)
plot(nVal_F1, 1.3*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F1, 0.7*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSAF1, r_sq_SAF1), 'y = x +/- 30%', '',  'location', 'best');
set(gca,'YLim', [0.0003, nMax_F1*1.1] )
set(gca,'XLim', [0.0003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_DioProjF1 = fitlm(Table_Dio_Proj{101:140, "Wt_Meas"}, Table_Dio_Proj{101:140, "Wt"}, 'y~-1+x1');
m_DioProjF1 = lm_DioProjF1.Coefficients.Estimate(1);
fitY_DioProjF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Dio_Proj{101:140, "Wt"}), max(Table_Dio_Proj{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_ProjF1 = lm_DioProjF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_DioProjF1(i) = m_DioProjF1 * nVal_F1(i);
end

subplot(1, 2, 2)
plot(Table_Dio_Proj{101:140, "Wt_Meas"}, Table_Dio_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP film terminal settling velocity to \n\r mP film terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Dioguardi et al (2018) using \n\r projection area of volume equivalent sphere as effective area.'))
hold on
plot(nVal_F1, nVal_F1, '-k', 'LineWidth', 1)
plot(nVal_F1, fitY_DioProjF1, '--g', 'LineWidth', 1)
plot(nVal_F1, 1.3*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F1, 0.7*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProjF1, r_sq_ProjF1), 'y = x +/- 30%', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F1*1.1] )
set(gca,'XLim', [0.003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Dioguardi/DioVM_MeasVsCalc_FitF1.jpg', 'Resolution', 1200);

%% Combine all m and r_sq values into the error table: Projected Area
Error_table = readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_Proj.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

Dio_rsq_proj = [r_sq_Proj; r_sq_ProjF3; r_sq_ProjF2; r_sq_ProjF1];
Dio_m_proj = [m_DioProj; m_DioProjF3; m_DioProjF2; m_DioProjF1];

Dio_Proj_Table = array2table([Dio_m_proj Dio_rsq_proj]);
Dio_Proj_Table.Properties.VariableNames = Col_names;
Dio_Proj_Table.Properties.RowNames = Row_names;

Error_table_Proj = [Error_table Dio_Proj_Table];

writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Dioguardi/DioFinalTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Dioguardi/DioFinalTableVM_Proj.xls', 'WriteRowNames', true);

%% Combine all m and r_sq values into the error table: Surface Area
Error_table = readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiErrorTableVM_SA.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

Dio_rsq_SA = [r_sq_SA; r_sq_SAF3; r_sq_SAF2; r_sq_SAF1];
Dio_m_SA = [m_DioSA; m_DioSAF3; m_DioSAF2; m_DioSAF1];

Dio_SA_Table = array2table([Dio_m_SA Dio_rsq_SA]);
Dio_SA_Table.Properties.VariableNames = Col_names;
Dio_SA_Table.Properties.RowNames = Row_names;

Error_table_SA = [Error_table Dio_SA_Table];

writetable(Error_table_SA, './DragModelsTest/Output/20220621/Dioguardi/DioFinalTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_SA, './DragModelsTest/Output/20220621/Dioguardi/DioFinalTableVM_SA.xls', 'WriteRowNames', true);
%% D1) Re against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Dio_SA.('Re_Meas'), Table_Dio_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_SA.('Re_Calc'), Table_Dio_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Dioguardi Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Dio_Proj.('Re_Meas'), Table_Dio_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_Proj.('Re_Calc'), Table_Dio_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Dioguardi Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ReVsCd.jpg', 'Resolution', 300);

%% D2) Re against Cd (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Dio_SA.('Re_Meas'), Table_Dio_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_SA{1:80, "Re_Calc"}, Table_Dio_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_SA{81:100, "Re_Calc"}, Table_Dio_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_SA{101:140, "Re_Calc"}, Table_Dio_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Dio_Proj.('Re_Meas'), Table_Dio_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_Proj{1:80, "Re_Calc"}, Table_Dio_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_Proj{81:100, "Re_Calc"}, Table_Dio_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "Re_Calc"}, Table_Dio_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ReVsCd_Shapes.jpg', 'Resolution', 300);

%% E) ESD against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Dioguardi Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Dioguardi Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ESDVsCd.jpg', 'Resolution', 300);

%% E2) wt against CSF (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Dio_SA.('ESD'), Table_Dio_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_SA{1:80, "ESD"}, Table_Dio_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_SA{81:100, "ESD"}, Table_Dio_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_SA{101:140, "ESD"}, Table_Dio_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Dio_Proj.('ESD'), Table_Dio_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Dio_Proj{1:80, "ESD"}, Table_Dio_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio_Proj{81:100, "ESD"}, Table_Dio_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "ESD"}, Table_Dio_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dioguardi Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dioguardi/DioguardiVM_ESDVsCd_Shapes.jpg', 'Resolution', 300);

