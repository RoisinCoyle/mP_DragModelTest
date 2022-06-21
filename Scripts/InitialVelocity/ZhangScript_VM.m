%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: ZhangScript: VM
% Date created: 23.04.22
% Date last mostified: 17.05.22
% Purpose: To test the implementation of the Zhang' drag model on a range of
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

%% Zhang and Choi's method 1
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

% Set timestep
timestep = 0.0002;

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
		
		Re_ZC(i, t) = abs((rho_p(i) * wvel_ZC(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_ZC(i,t) = (58.58*(shape_ASF(i)^0.1936))/(Re_ZC(i)^0.8273);
	
		Fd_ZC(i,t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_ZC(i,t)^2.0)*Cd_ZC(i,t);
	
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
		    CdFinal_ZC(i) = (24.0/ReFinal_ZC(i))*(((1.0-shape_del(i))/(ReFinal_ZC(i)+1.0))^0.25) ...
			     + (24.0/ReFinal_ZC(i))*0.1806*(ReFinal_ZC(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_ZC(i)^0.08))) ...
			     + 0.4251/(1.0+((6880.95/ReFinal_ZC(i))*(shape_del(i)^5.05)));
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
Results_ZC = zeros(140, 10);

for i=1:140
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = wtFinal_ZC(i);
    Results_ZC(i, 4) = wvel_meas(i);
    Results_ZC(i, 5) = FinalTime_ZC(i);
    Results_ZC(i, 6) = DistTot_ZC(i);
    Results_ZC(i, 7) = timestep;
    Results_ZC(i, 8) = Reynolds(i);
    Results_ZC(i, 9) = ReFinal_ZC(i);
    Results_ZC(i, 10) = Cd_meas(i);
    Results_ZC(i, 11) = CdFinal_ZC(i);
    
end 

Table_ZC_SA = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_ZC_SA = [VM_Dataset.Shape Table_ZC_SA];
Table_ZC_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_SA, './DragModelsTest/Output/20220517/ZhangOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_SA, './DragModelsTest/Output/20220517/ZhangOutputVM_SA.xls', 'WriteRowNames', true);


%% Zhang Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area, using Newtons Drag formula which assumes a spherical shape

% Set timestep
timestep = 0.0002;

% Set up variable arrays
Re_ZC = zeros(140, 10000);
shape_ASF = zeros(140, 1);
beta2 = zeros(140,1);
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
		
		Re_ZC(i, t) = abs((rho_p(i) * wvel_ZC(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_ZC(i,t) = (58.58*(shape_ASF(i)^0.1936))/(Re_ZC(i)^0.8273);
	
		Fd_ZC(i,t) = 0.5*rho_f(i)*ProjA_ESD(i)*(wvel_ZC(i,t)^2.0)*Cd_ZC(i,t);
	
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
		    CdFinal_ZC(i) = (24.0/ReFinal_ZC(i))*(((1.0-shape_del(i))/(ReFinal_ZC(i)+1.0))^0.25) ...
			     + (24.0/ReFinal_ZC(i))*0.1806*(ReFinal_ZC(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_ZC(i)^0.08))) ...
			     + 0.4251/(1.0+((6880.95/ReFinal_ZC(i))*(shape_del(i)^5.05)));
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
Results_ZC = zeros(140, 11);

for i=1:140
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = wtFinal_ZC(i);
    Results_ZC(i, 4) = wvel_meas(i);
    Results_ZC(i, 5) = FinalTime_ZC(i);
    Results_ZC(i, 6) = DistTot_ZC(i);
    Results_ZC(i, 7) = timestep;
    Results_ZC(i, 8) = Reynolds(i);
    Results_ZC(i, 9) = ReFinal_ZC(i);
    Results_ZC(i, 10) = Cd_meas(i);
    Results_ZC(i, 11) = CdFinal_ZC(i);
    
end 

Table_ZC_Proj = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_ZC_Proj = [VM_Dataset.Shape Table_ZC_Proj];
Table_ZC_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_Proj, './DragModelsTest/Output/20220517/ZhangOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_Proj, './DragModelsTest/Output/20220517/ZhangOutputVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Zhang output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Zhang_SA= readtable("./DragModelsTest/Output/20220517/ZhangOutputVM_SA.txt", "Delimiter", ",");
Table_Zhang_Proj= readtable("./DragModelsTest/Output/20220517/ZhangOutputVM_Proj.txt", "Delimiter", ",");

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
title('Zhang and Choi Model. Using Particle Surface Area')
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
title('Zhang and Choi Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

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
title('Zhang and Choi Model. Using Particle Surface Area')
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
title('Zhang and Choi Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

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
title('Zhang and Choi Model. Using Particle Surface Area')
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
title('Zhang and Choi Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

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
title('Zhang and Choi Model. Using Particle Surface Area')
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
title('Zhang and Choi Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

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
title('Zhang Model. Using Particle Surface Area')
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
plot(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Zhang Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_Proj*1.1] )
set(gca,'XLim', [0, MaxW_Proj*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_MeasVsCalc.jpg', 'Resolution', 300);

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
title('Zhang Model. Using Particle Surface Area')
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
title('Zhang Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_Proj])
set(gca, 'Xlim', [0, 1.1*MaxW_Proj])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% D2) wt against wt measured using Matlab fitlm function
% ========================================================

% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_MeasVsCalc_Fit.jpg', 'Resolution', 300);

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
title('ZCguardi Model. Using Particle Surface Area')
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
title('Zhang Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_ReVsCd.jpg', 'Resolution', 300);

%% E2) wt against CSF (SHAPES)
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
exportgraphics(gcf, './DragModelsTest/Output/20220517/ZhangVM_ReVsCd_Shapes.jpg', 'Resolution', 300);

