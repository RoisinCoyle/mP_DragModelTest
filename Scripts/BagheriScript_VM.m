%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: BagheriScript: VM
% Date created: 23.04.22
% Date last mostified: 17.05.22
% Purpose: To test the implementation of the Bagheri drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file

% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
VM_Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

rho_p = table2array(VM_Dataset(:, "ParticleDensity"));
rho_f = table2array(VM_Dataset(:, "FluidDensity"));
vis_dyn = table2array(VM_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(VM_Dataset(:, "KinematicVisvosity"));

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
		
		Cd_BB(i,t) = ((24.0*Correction_S(i))/Re_BB(i,t))*(1+ 0.125*(((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))*(2.0/3.0))) ...
			+ (0.46*Correction_N(i))/(1 + (5330/((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))));
	
		Fd_BB(i,t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_BB(i,t)^2.0)*Cd_BB(i,t);
	
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

writetable(Table_BB_SA, './DragModelsTest/Output/20220517/BagheriOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_SA, './DragModelsTest/Output/20220517/BagheriOutputVM_SA.xls', 'WriteRowNames', true);


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
		
		Cd_BB(i,t) = ((24.0*Correction_S(i))/Re_BB(i,t))*(1+ 0.125*(((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))*(2.0/3.0))) ...
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

writetable(Table_BB_Proj, './DragModelsTest/Output/20220517/BagheriOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_Proj, './DragModelsTest/Output/20220517/BagheriOutputVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Bagheri output
% <<<<<<<<<<<<<<<<<<<
Table_BB_SA= readtable("./DragModelsTest/Output/20220517/BagheriOutputVM_SA.txt", "Delimiter", ",");
Table_BB_Proj= readtable("./DragModelsTest/Output/20220517/BagheriOutputVM_Proj.txt", "Delimiter", ",");

%% A1) wt against ESD (ALL)
% ==========================

% Method 1: All
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

% Method 2: All
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('ESD'), Table_BB_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_SA{1:80, "ESD"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "ESD"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "ESD"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_BB_Proj.('ESD'), Table_BB_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_Proj{1:80, "ESD"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "ESD"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "ESD"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF (ALL)
% ==========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_BB_SA.('CSF'), Table_BB_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_BB_SA.('CSF'), Table_BB_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
plot(Table_BB_Proj.('CSF'), Table_BB_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB_Proj{1:80, "CSF"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "CSF"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "CSF"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest_SA(1) = max(Table_BB_SA.Wt_Calc);
Highest_SA(2) = max(Table_BB_SA.Wt_Meas);
Highest_Proj(1) =  max(Table_BB_Proj.Wt_Calc);
Highest_Proj(2) = max(Table_BB_Proj.Wt_Meas);
MaxW_SA = max(Highest_SA);
MaxW_Proj = max(Highest_Proj);
yx_SA=linspace(0, MaxW_SA, 100);
yx_Proj=linspace(0, MaxW_Proj, 100);

% Method 1: Plot shapes separately
subplot(1, 2, 1)
plot(yx_SA, yx_SA)
hold on
plot(Table_BB_SA{1:80, "Wt_Meas"}, Table_BB_SA{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_SA{81:100, "Wt_Meas"}, Table_BB_SA{81:100, "Wt_Calc"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_SA{101:140, "Wt_Meas"}, Table_BB_SA{101:140, "Wt_Calc"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Bagheri Model. Using Particle Surface Area')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_SA*1.1] )
set(gca,'XLim', [0, MaxW_SA*1.1] )
hold off

% Method 2: Plot shapes separately
subplot(1, 2, 2)
plot(yx_Proj, yx_Proj)
hold on
plot(Table_BB_Proj{1:80, "Wt_Meas"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB_Proj{81:100, "Wt_Meas"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "Wt_Meas"}, Table_BB_Proj{101:140, "Wt_Calc"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_Proj*1.1] )
set(gca,'XLim', [0, MaxW_Proj*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
% ===============================================

subplot(1, 2, 1)
plot(Table_BB_SA.('Wt_Meas'), Table_BB_SA.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_SA, yx_SA, '-k')
p=polyfit(Table_BB_SA.('Wt_Meas'), Table_BB_SA.('Wt_Calc'), 1);
px=[min(Table_BB_SA.('Wt_Meas')) max(Table_BB_SA.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.6*px(2), py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_BB_SA.("Wt_Meas")\Table_BB_SA.("Wt_Calc");
mx = m*Table_BB_SA.("Wt_Meas");
plot(Table_BB_SA.('Wt_Meas'), mx, '-g');
text(0.75*px(2), 0.7*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Bagheri Model. Using Particle Surface Area')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_SA])
set(gca, 'Xlim', [0, 1.1*MaxW_SA])
hold off

subplot(1, 2, 2)
plot(Table_BB_Proj.('Wt_Meas'), Table_BB_Proj.('Wt_Calc'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_Proj, yx_Proj, '-k')
p=polyfit(Table_BB_Proj.('Wt_Meas'), Table_BB_Proj.('Wt_Calc'), 1);
px=[min(Table_BB_Proj.('Wt_Meas')) max(Table_BB_Proj.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.75*px(2), 0.7*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_BB_Proj.("Wt_Meas")\Table_BB_Proj.("Wt_Calc");
mx = m*Table_BB_Proj.("Wt_Meas");
plot(Table_BB_Proj.('Wt_Meas'), mx, '-g');
text(0.6*px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_Proj])
set(gca, 'Xlim', [0, 1.1*MaxW_Proj])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% E1) Re against Cd (ALL)
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
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_ReVsCd.jpg', 'Resolution', 300);

%% B2) wt against CSF (SHAPES)
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
title('Bagheri Model. Using Particle Surface Area')
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
title('Bagheri Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/BagheriVM_ReVsCd_Shapes.jpg', 'Resolution', 300);