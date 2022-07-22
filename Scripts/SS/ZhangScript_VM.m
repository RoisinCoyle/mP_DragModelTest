%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: ZhangScript: VM
% Date created: 23.04.22
% Date last mostified: 22.07.22
% Purpose: To test the implementation of the Zhang drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
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

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_ZC(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
    hold on
    plot(timesec(1:(FinalStep_ZC(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_ZC = zeros(140, 12);

for i=1:140
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = shape_ASF(i);
    Results_ZC(i, 4) = wtFinal_ZC(i);
    Results_ZC(i, 5) = wvel_meas(i);
    Results_ZC(i, 6) = FinalTime_ZC(i);
    Results_ZC(i, 7) = DistTot_ZC(i);
    Results_ZC(i, 8) = timestep;
    Results_ZC(i, 9) = Reynolds(i);
    Results_ZC(i, 10) = ReFinal_ZC(i);
    Results_ZC(i, 11) = Cd_meas(i);
    Results_ZC(i, 12) = CdFinal_ZC(i);
    
end 

Table_ZC_SA = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'ASF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_ZC_SA = [VM_Dataset.Shape Table_ZC_SA];
Table_ZC_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_SA.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error(i) = ((residual(i) / wvel_meas(i))*100);
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
    residual_F3(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i) / wvel_meas(i))*100);
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
    residual_F2(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i) / wvel_meas(i))*100);
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
    residual_F1(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error_F1(i) = ((residual_F1(i) / wvel_meas(i))*100);
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

writetable(Error_table, './DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_SA.xls', 'WriteRowNames', true);

%% Zhang Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area, using Newtons Drag formula

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

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:140
    subplot(14, 10, n)
    meas = zeros(FinalStep_ZC(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
    hold on
    plot(timesec(1:(FinalStep_ZC(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_ZC = zeros(140, 12);

for i=1:140
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = shape_ASF(i);
    Results_ZC(i, 4) = wtFinal_ZC(i);
    Results_ZC(i, 5) = wvel_meas(i);
    Results_ZC(i, 6) = FinalTime_ZC(i);
    Results_ZC(i, 7) = DistTot_ZC(i);
    Results_ZC(i, 8) = timestep;
    Results_ZC(i, 9) = Reynolds(i);
    Results_ZC(i, 10) = ReFinal_ZC(i);
    Results_ZC(i, 11) = Cd_meas(i);
    Results_ZC(i, 12) = CdFinal_ZC(i);
    
end 

Table_ZC_Proj = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'ASF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_ZC_Proj = [VM_Dataset.Shape Table_ZC_Proj];
Table_ZC_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_Proj.xls', 'WriteRowNames', true);

%% Distance assumption calculation and plot

DistConst_ZC = zeros(140, 1);

for i=1:140
    DistConst_ZC(i) = FinalTime_ZC(i) * wtFinal_ZC(i);
end

% Fit linear model through the intercept: SA
lm_ZCDist = fitlm(DistTot_ZC, DistConst_ZC, 'y~-1+x1');
m_ZCDist = lm_ZCDist.Coefficients.Estimate(1);
fitY_ZCDist = zeros(1000, 1);
% Generate data using linear model:
n1=[max(DistTot_ZC), max(DistConst_ZC)] ;
nMax = max(n1);
nVal=linspace(0.00001, nMax, 1000);
r_sq_Dist = lm_ZCDist.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCDist(i) = m_ZCDist * nVal(i);
end

subplot(1, 2, 2)
plot(DistTot_ZC, DistConst_ZC, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[1, 1, 0]')
ylabel('Distance travelled at constant velocity (m)')
xlabel('Distance travelled in attaining terminal velocity (m)')
title('Zhang and Choi (2021): Using Particle Projected Area.')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCDist, '--k')
plot(nVal, 0.7*nVal, ':k')
plot(nVal, 1.3*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCDist, r_sq_Dist), '', '', 'location', 'best');
set(gca,'YLim', [0.00005, nMax*1.3] )
set(gca,'XLim', [0.00005, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'Maximized')
exportgraphics(gcf, './DragModelsTest/Output/20220621/Distance/ZC_DistanceProj.jpg', 'Resolution', 300)
%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error(i) = ((residual(i) / wvel_meas(i))*100);
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
    residual_F3(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error_F3(i) = ((residual_F3(i) / wvel_meas(i))*100);
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
    residual_F2(i) = (wtFinal_ZC(i) - wvel_meas(i));
    Percentage_Error_F2(i) = ((residual_F2(i) / wvel_meas(i))*100);
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
    residual_F1(i) = (wtFinal_ZC(i) - wvel_meas(i));
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

writetable(Error_table, './DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Zhang output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Zhang_SA= readtable("./DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_SA.txt", "Delimiter", ",");
Table_Zhang_Proj= readtable("./DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_Proj.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

% Method 1: All
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

% Method 2: All
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ESDVsWt.jpg', 'Resolution', 300)

%% A2) wt against ESD: Shapes separately
% =========================================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_SA{1:80, "ESD"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "ESD"}, Table_Zhang_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "ESD"}, Table_Zhang_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_Proj{1:80, "ESD"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "ESD"}, Table_Zhang_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "ESD"}, Table_Zhang_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Zhang_SA.('CSF'), Table_Zhang_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_SA.('CSF'), Table_Zhang_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('CSF'), Table_Zhang_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_Proj.('CSF'), Table_Zhang_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF: Shapes separate
% ======================================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Zhang_SA.('CSF'), Table_Zhang_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_SA{1:80, "CSF"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "CSF"}, Table_Zhang_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "CSF"}, Table_Zhang_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('CSF'), Table_Zhang_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Zhang_Proj{1:80, "CSF"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "CSF"}, Table_Zhang_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "CSF"}, Table_Zhang_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest_SA(1) = max(Table_Zhang_SA.Wt);
Highest_SA(2) = max(Table_Zhang_SA.Wt_Meas);
Highest_Proj(1) =  max(Table_Zhang_Proj.Wt);
Highest_Proj(2) = max(Table_Zhang_Proj.Wt_Meas);
MaxW_SA = max(Highest_SA);
MaxW_Proj = max(Highest_Proj);
yx_SA=linspace(0, MaxW_SA, 100);
yx_Proj=linspace(0, MaxW_Proj, 100);

% Method 1: Plot shapes separately
subplot(1, 2, 1)
plot(yx_SA, yx_SA)
hold on
plot(Table_Zhang_SA{1:80, "Wt_Meas"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "Wt_Meas"}, Table_Zhang_SA{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "Wt_Meas"}, Table_Zhang_SA{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Zhang Model. Using Particle Surface Area.')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('Measured=Calculated', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_SA*1.1] )
set(gca,'XLim', [0, MaxW_SA*1.1] )
hold off

% Method 2: Plot shapes separately
subplot(1, 2, 2)
plot(yx_Proj, yx_Proj)
hold on
plot(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Zhang Model. Estimated projection area using max CSA.')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('Measured=Calculated', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_Proj*1.1] )
set(gca,'XLim', [0, MaxW_Proj*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
% ===============================================

subplot(1, 2, 1)
plot(Table_Zhang_SA.('Wt_Meas'), Table_Zhang_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_SA, yx_SA, '-k')
p=polyfit(Table_Zhang_SA.('Wt_Meas'), Table_Zhang_SA.('Wt'), 1);
px=[min(Table_Zhang_SA.('Wt_Meas')) max(Table_Zhang_SA.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.5*px(2), 1.2*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Zhang_SA.("Wt_Meas")\Table_Zhang_SA.("Wt");
mx = m*Table_Zhang_SA.("Wt_Meas");
plot(Table_Zhang_SA.('Wt_Meas'), mx, '-g');
text(0.5*px(2), 1.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Zhang Model. Using Particle Surface Area.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_SA])
set(gca, 'Xlim', [0, 1.1*MaxW_SA])
hold off

subplot(1, 2, 2)
plot(Table_Zhang_Proj.('Wt_Meas'), Table_Zhang_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_Proj, yx_Proj, '-k')
p=polyfit(Table_Zhang_Proj.('Wt_Meas'), Table_Zhang_Proj.('Wt'), 1);
px=[min(Table_Zhang_Proj.('Wt_Meas')) max(Table_Zhang_Proj.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.4*px(2), 1.1*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Zhang_Proj.("Wt_Meas")\Table_Zhang_Proj.("Wt");
mx = m*Table_Zhang_Proj.("Wt_Meas");
plot(Table_Zhang_Proj.('Wt_Meas'), mx, '-g');
text(0.4*px(2), max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Zhang Model. Estimated projection area using max CSA.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_Proj])
set(gca, 'Xlim', [0, 1.1*MaxW_Proj])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% D2) wt against wt measured using Matlab fitlm function
% ========================================================

% D2 A) All shapes

% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model. Using particle Surface Area')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCSA, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq_SA), 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.1] )
set(gca,'XLim', [0.003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_Proj = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model. Estimated projection area using max CSA.')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCProj, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq_Proj), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.1] )
set(gca,'XLim', [0.003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_FitAll.jpg', 'Resolution', 300);

%% D2 B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_SA = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA{1:80, "Wt_Meas"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model. Using particle Surface Area.')
hold on
plot(Table_Zhang_SA{81:100, "Wt_Meas"}, Table_Zhang_SA{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "Wt_Meas"}, Table_Zhang_SA{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCSA, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq_SA), 'location', 'best');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq_Proj = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model. Estimated projection area using max CSA.')
hold on
plot(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCProj, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq_Proj), 'location', 'best');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 300);

%% D2 C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_ZCSAF3 = fitlm(Table_Zhang_SA{1:80, "Wt_Meas"}, Table_Zhang_SA{1:80, "Wt"}, 'y~-1+x1');
m_ZCSAF3 = lm_ZCSAF3.Coefficients.Estimate(1);
fitY_ZCSAF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Zhang_SA{1:80, "Wt"}), max(Table_Zhang_SA{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_SAF3 = lm_ZCSAF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCSAF3(i) = m_ZCSAF3 * nVal_F3(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA{1:80, "Wt_Meas"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Using Particle Surface Area.')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_ZCSAF3, '--b')
plot(nVal_F3, 1.3*nVal_F3, ':k')
plot(nVal_F3, 0.7*nVal_F3, ':k')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSAF3, r_sq_SAF3), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProjF3 = fitlm(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'y~-1+x1');
m_ZCProjF3 = lm_ZCProjF3.Coefficients.Estimate(1);
fitY_ZCProjF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Zhang_Proj{1:80, "Wt"}), max(Table_Zhang_Proj{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_ProjF3 = lm_ZCProjF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCProjF3(i) = m_ZCProjF3 * nVal_F3(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Estimated projection area using max CSA.')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_ZCProjF3, '--b')
plot(nVal_F3, 1.3*nVal_F3, ':k')
plot(nVal_F3, 0.7*nVal_F3, ':k')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProjF3, r_sq_ProjF3), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_FitF3.jpg', 'Resolution', 300);

%% D2 D) Plot fibres separately with fitted model

% Fit linear model through the intercept: SA
lm_ZCSAF2 = fitlm(Table_Zhang_SA{81:100, "Wt_Meas"}, Table_Zhang_SA{81:100, "Wt"}, 'y~-1+x1');
m_ZCSAF2 = lm_ZCSAF2.Coefficients.Estimate(1);
fitY_ZCSAF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Zhang_SA{81:100, "Wt"}), max(Table_Zhang_SA{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_SAF2 = lm_ZCSAF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCSAF2(i) = m_ZCSAF2 * nVal_F2(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA{81:100, "Wt_Meas"}, Table_Zhang_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Using Particle Surface Area.')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_ZCSAF2, '--r')
plot(nVal_F2, 1.3*nVal_F2, ':k')
plot(nVal_F2, 0.7*nVal_F2, ':k')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSAF2, r_sq_SAF2), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F2*1.1] )
set(gca,'XLim', [0.003, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProjF2 = fitlm(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'y~-1+x1');
m_ZCProjF2 = lm_ZCProjF2.Coefficients.Estimate(1);
fitY_ZCProjF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Zhang_Proj{81:100, "Wt"}), max(Table_Zhang_Proj{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_ProjF2 = lm_ZCProjF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCProjF2(i) = m_ZCProjF2 * nVal_F2(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Estimated projection area using max CSA.')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_ZCProjF2, '--r')
plot(nVal_F2, 1.3*nVal_F2, ':k')
plot(nVal_F2, 0.7*nVal_F2, ':k')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProjF2, r_sq_ProjF2), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F2*1.1] )
set(gca,'XLim', [0.003, nMax_F2*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_FitF2.jpg', 'Resolution', 300);

%% D2 E) Plot film separately with fitted model

% Fit linear model through the intercept: SA
lm_ZCSAF1 = fitlm(Table_Zhang_SA{101:140, "Wt_Meas"}, Table_Zhang_SA{101:140, "Wt"}, 'y~-1+x1');
m_ZCSAF1 = lm_ZCSAF1.Coefficients.Estimate(1);
fitY_ZCSAF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Zhang_SA{101:140, "Wt"}), max(Table_Zhang_SA{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_SAF1 = lm_ZCSAF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCSAF1(i) = m_ZCSAF1 * nVal_F1(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA{101:140, "Wt_Meas"}, Table_Zhang_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Using Particle Surface Area.')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_ZCSAF1, '--g')
plot(nVal_F1, 1.3*nVal_F1, ':k')
plot(nVal_F1, 0.7*nVal_F1, ':k')
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSAF1, r_sq_SAF1), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F1*1.1] )
set(gca,'XLim', [0.003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProjF1 = fitlm(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, 'y~-1+x1');
m_ZCProjF1 = lm_ZCProjF1.Coefficients.Estimate(1);
fitY_ZCProjF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Zhang_Proj{101:140, "Wt"}), max(Table_Zhang_Proj{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_ProjF1 = lm_ZCProjF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_ZCProjF1(i) = m_ZCProjF1 * nVal_F1(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model: Estimated projection area using max CSA.')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_ZCProjF1, '--g')
plot(nVal_F1, 1.3*nVal_F1, ':k')
plot(nVal_F1, 0.7*nVal_F1, ':k')
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProjF1, r_sq_ProjF1), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F1*1.1] )
set(gca,'XLim', [0.003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_MeasVsCalc_FitF1.jpg', 'Resolution', 300);

%% Combine all m and r_sq values into the error table: Projected Area
Error_table = readtable("./DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_Proj.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

ZC_rsq_proj = [r_sq_Proj; r_sq_ProjF3; r_sq_ProjF2; r_sq_ProjF1];
ZC_m_proj = [m_ZCProj; m_ZCProjF3; m_ZCProjF2; m_ZCProjF1];

ZC_Proj_Table = array2table([ZC_m_proj ZC_rsq_proj]);
ZC_Proj_Table.Properties.VariableNames = Col_names;
ZC_Proj_Table.Properties.RowNames = Row_names;

Error_table_Proj = [Error_table ZC_Proj_Table];

writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Zhang/ZhangFinalTableVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_Proj, './DragModelsTest/Output/20220621/Zhang/ZhangFinalTableVM_Proj.xls', 'WriteRowNames', true);

%% Combine all m and r_sq values into the error table: Surface Area
Error_table = readtable("./DragModelsTest/Output/20220621/Zhang/ZhangErrorTableVM_SA.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

ZC_rsq_SA = [r_sq_SA; r_sq_SAF3; r_sq_SAF2; r_sq_SAF1];
ZC_m_SA = [m_ZCSA; m_ZCSAF3; m_ZCSAF2; m_ZCSAF1];

ZC_SA_Table = array2table([ZC_m_SA ZC_rsq_SA]);
ZC_SA_Table.Properties.VariableNames = Col_names;
ZC_SA_Table.Properties.RowNames = Row_names;

Error_table_SA = [Error_table ZC_SA_Table];

writetable(Error_table_SA, './DragModelsTest/Output/20220621/Zhang/ZhangFinalTableVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table_SA, './DragModelsTest/Output/20220621/Zhang/ZhangFinalTableVM_SA.xls', 'WriteRowNames', true);
%% E1) Re against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Zhang_SA.('Re_Meas'), Table_Zhang_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_SA.('Re_Calc'), Table_Zhang_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Zhang Model. Using Particle Surface Area.')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('Re_Meas'), Table_Zhang_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_Proj.('Re_Calc'), Table_Zhang_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Zhang Model. Estimated projection area using max CSA.')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ReVsCd.jpg', 'Resolution', 300);


%% E2) Cd against Re (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Zhang_SA.('Re_Meas'), Table_Zhang_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_SA{1:80, "Re_Calc"}, Table_Zhang_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "Re_Calc"}, Table_Zhang_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "Re_Calc"}, Table_Zhang_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('Re_Meas'), Table_Zhang_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_Proj{1:80, "Re_Calc"}, Table_Zhang_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "Re_Calc"}, Table_Zhang_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "Re_Calc"}, Table_Zhang_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ReVsCd_Shapes.jpg', 'Resolution', 300);

%% F1) ESD against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Zhang Model. Using Particle Surface Area.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')

hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Zhang Model. Estimated projection area using max CSA.')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ESDVsCd.jpg', 'Resolution', 300);


%% F2) ESD against Re (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ESD'), Table_Zhang_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_SA{1:80, "ESD"}, Table_Zhang_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "ESD"}, Table_Zhang_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "ESD"}, Table_Zhang_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ESD'), Table_Zhang_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Zhang_Proj{1:80, "ESD"}, Table_Zhang_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "ESD"}, Table_Zhang_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "ESD"}, Table_Zhang_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('ESD (m)')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ESDVsCd_Shapes.jpg', 'Resolution', 300);

%% G1) wt against ASF
% ====================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ASF'), Table_Zhang_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Zhang_SA.('ASF'), Table_Zhang_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('ASF')
set(gca, 'XScale', 'log')
hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ASF'), Table_Zhang_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Zhang_Proj.('ASF'), Table_Zhang_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('ASF')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ASFVsW.jpg', 'Resolution', 300);

%% G2) wt against ASF: Shapes separate
% ======================================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Zhang_SA.('ASF'), Table_Zhang_SA.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Zhang_SA{1:80, "ASF"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_SA{81:100, "ASF"}, Table_Zhang_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "ASF"}, Table_Zhang_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Using Particle Surface Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('ASF')
set(gca, 'XScale', 'log')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Zhang_Proj.('ASF'), Table_Zhang_Proj.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Zhang_Proj{1:80, "ASF"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "ASF"}, Table_Zhang_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "ASF"}, Table_Zhang_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Zhang and Choi Model. Estimated projection area using max CSA.')
ylabel('Terminal settling velocity (m/s)')
xlabel('ASF')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Zhang/ZhangVM_ASFVsW_Shapes.jpg', 'Resolution', 300);