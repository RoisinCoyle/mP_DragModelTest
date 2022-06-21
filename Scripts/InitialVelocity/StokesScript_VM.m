%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: StokesScript: VM
% Date created: 23.04.22
% Date last mostified: 17.05.22
% Purpose: To test the implementation of the Stokes drag model on a range of
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

%% Stokes' method 1
% <<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

% Set timestep
timestep = 0.0002;

% Set initial velocity and timestep
 wvel_Stokes = zeros(140, 10000);
 wvel_Stokes(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
% Set up variable arrays
Cd_Stokes = zeros(140, 10000);
CdFinal_Stokes = zeros(140, 1);
Re_Stokes = zeros(140, 10000);
ReFinal_Stokes = zeros(140, 1);
Fd_Stokes = zeros(140, 10000);
Fg_Stokes = zeros(140, 10000);
Fb_Stokes = zeros(140, 10000);
Fnet_Stokes = zeros(140, 10000);
Dist1_Stokes = zeros(140, 10000);
Acc_Stokes = zeros(140, 10000);
wtFinal_Stokes = zeros(140, 1);
DistTot_Stokes = zeros(140, 1);
FinalStep_Stokes = zeros(140, 1);
FinalTime_Stokes = zeros(140, 1);

% Begin calculation
for i=1:140
    for t=1:10000
        
        Re_Stokes(i, t) = (rho_f(i)*wvel_Stokes(i, t)*d_equi(i))/vis_dyn(i);
        Cd_Stokes(t) = 24.0/Re_Stokes(i, t);
	
		Fd_Stokes(t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_Stokes(i,t)^2.0)*Cd_Stokes(t);
	
		Fg_Stokes(t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Stokes(t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Stokes(t) = Fg_Stokes(t) - Fb_Stokes(t) - Fd_Stokes(t);
	
		wvel_Stokes(i, t+1) = ((Fnet_Stokes(t)/Mass_mP(i))*timestep)+wvel_Stokes(i, t);
	
		Dist1_Stokes(i, t) = wvel_Stokes(i, t) * timestep;
		DistTot_Stokes(i) = DistTot_Stokes(i) + Dist1_Stokes(i, t);
		Acc_Stokes(t) = (wvel_Stokes(i, t+1) - wvel_Stokes(i, t))/timestep;
        
        if (Acc_Stokes(t)< 0.001)
			FinalStep_Stokes(i) = (t+1);
            FinalTime_Stokes(i) = (t+1)*timestep;
            wtFinal_Stokes(i)=wvel_Stokes(i, t+1);
            Dist1_Stokes(i, t+1) = wvel_Stokes(i, t+1) * timestep;
            DistTot_Stokes(i) = DistTot_Stokes(i) + Dist1_Stokes(i, t+1);
            ReFinal_Stokes(i) = abs((rho_p(i) * wvel_Stokes(i, t+1) * d_equi(i))/ vis_dyn(i));
		    CdFinal_Stokes(i) = (24.0/ReFinal_Stokes(i))*(((1.0-shape_del(i))/(ReFinal_Stokes(i)+1.0))^0.25) ...
			     + (24.0/ReFinal_Stokes(i))*0.1806*(ReFinal_Stokes(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_Stokes(i)^0.08))) ...
			     + 0.4251/(1.0+((6880.95/ReFinal_Stokes(i))*(shape_del(i)^5.05)));
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
    meas = zeros(FinalStep_Stokes(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
    hold on
    plot(timesec(1:(FinalStep_Stokes(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_Stokes = zeros(140, 11);

for i=1:140
    Results_Stokes(i, 1) = d_equi(i);
    Results_Stokes(i, 2) = CSF(i);
    Results_Stokes(i, 3) = wtFinal_Stokes(i);
    Results_Stokes(i, 4) = wvel_meas(i);
    Results_Stokes(i, 5) = FinalTime_Stokes(i);
    Results_Stokes(i, 6) = DistTot_Stokes(i);
    Results_Stokes(i, 7) = timestep;
    Results_Stokes(i, 8) = Reynolds(i);
    Results_Stokes(i, 9) = ReFinal_Stokes(i);
    Results_Stokes(i, 10) = Cd_meas(i);
    Results_Stokes(i, 11) = CdFinal_Stokes(i);
end 

Table_Stokes_SA = array2table(Results_Stokes, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_Stokes_SA = [VM_Dataset.Shape Table_Stokes_SA];
Table_Stokes_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Stokes_SA, './DragModelsTest/Output/20220517/StokesOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Stokes_SA, './DragModelsTest/Output/20220517/StokesOutputVM_SA.xls', 'WriteRowNames', true);


%% Stokes Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area in the calculation of the drag force.

% Set timestep
timestep = 0.0005;

% Set initial velocity and timestep
 wvel_Stokes = zeros(140, 10000);
 wvel_Stokes(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
% Set up variable arrays
Cd_Stokes = zeros(140, 10000);
CdFinal_Stokes = zeros(140, 1);
Re_Stokes = zeros(140, 10000);
ReFinal_Stokes = zeros(140, 1);
Fd_Stokes = zeros(140, 10000);
Fg_Stokes = zeros(140, 10000);
Fb_Stokes = zeros(140, 10000);
Fnet_Stokes = zeros(140, 10000);
Dist1_Stokes = zeros(140, 10000);
Acc_Stokes = zeros(140, 10000);
wtFinal_Stokes = zeros(140, 1);
DistTot_Stokes = zeros(140, 1);
FinalStep_Stokes = zeros(140, 1);
FinalTime_Stokes = zeros(140, 1);

for i=1:140
    for t=1:10000
        
        Re_Stokes(i, t) = (rho_f(i)*wvel_Stokes(i, t)*d_equi(i))/vis_dyn(i);
        Cd_Stokes(t) = 24.0/Re_Stokes(i, t);
	
		Fd_Stokes(t) = 0.5*rho_f(i)*ProjA_ESD(i)*(wvel_Stokes(i,t)^2.0)*Cd_Stokes(t);
	
		Fg_Stokes(t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Stokes(t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Stokes(t) = Fg_Stokes(t) - Fb_Stokes(t) - Fd_Stokes(t);
	
		wvel_Stokes(i, t+1) = ((Fnet_Stokes(t)/Mass_mP(i))*timestep)+wvel_Stokes(i, t);
	
		Dist1_Stokes(i, t) = wvel_Stokes(i, t) * timestep;
		DistTot_Stokes(i) = DistTot_Stokes(i) + Dist1_Stokes(i);
		Acc_Stokes(t) = (wvel_Stokes(i, t+1) - wvel_Stokes(i, t))/timestep;
        
        if (Acc_Stokes(t)< 0.001)
			FinalStep_Stokes(i) = (t+1);
            FinalTime_Stokes(i) = (t+1)*timestep;
            wtFinal_Stokes(i)=wvel_Stokes(i, t+1);
            Dist1_Stokes(i, t+1) = wvel_Stokes(i, t+1) * timestep;
            DistTot_Stokes(i) = DistTot_Stokes(i) + Dist1_Stokes(i, t+1);
            ReFinal_Stokes(i) = abs((rho_p(i) * wvel_Stokes(i, t+1) * d_equi(i))/ vis_dyn(i));
		    CdFinal_Stokes(i) = (24.0/ReFinal_Stokes(i))*(((1.0-shape_del(i))/(ReFinal_Stokes(i)+1.0))^0.25) ...
			     + (24.0/ReFinal_Stokes(i))*0.1806*(ReFinal_Stokes(i)^0.6459)*(shape_del(i)^(-1.0*(ReFinal_Stokes(i)^0.08))) ...
			     + 0.4251/(1.0+((6880.95/ReFinal_Stokes(i))*(shape_del(i)^5.05)));
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
    meas = zeros(FinalStep_Stokes(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
    hold on
    plot(timesec(1:(FinalStep_Stokes(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_Stokes_Proj = zeros(140, 11);

for i=1:140
    Results_Stokes_Proj(i, 1) = d_equi(i);
    Results_Stokes_Proj(i, 2) = CSF(i);
    Results_Stokes_Proj(i, 3) = wtFinal_Stokes(i);
    Results_Stokes_Proj(i, 4) = wvel_meas(i);
    Results_Stokes_Proj(i, 5) = FinalTime_Stokes(i);
    Results_Stokes_Proj(i, 6) = DistTot_Stokes(i);
    Results_Stokes_Proj(i, 7) = timestep;
    Results_Stokes_Proj(i, 8) = Reynolds(i);
    Results_Stokes_Proj(i, 9) = ReFinal_Stokes(i);
    Results_Stokes_Proj(i, 10) = Cd_meas(i);
    Results_Stokes_Proj(i, 11) = CdFinal_Stokes(i);
end 

Table_Stokes_Proj = array2table(Results_Stokes_Proj, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time', ...
    'Distance', 'Timestep', 'Re_Meas', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_Stokes_Proj = [VM_Dataset.Shape Table_Stokes_Proj];
Table_Stokes_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Stokes_Proj, './DragModelsTest/Output/20220517/StokesOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Stokes_Proj, './DragModelsTest/Output/20220517/StokesOutputVM_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Stokes output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Stokes_SA= readtable("./DragModelsTest/Output/20220517/StokesOutputVM_SA.txt", "Delimiter", ",");
Table_Stokes_Proj= readtable("./DragModelsTest/Output/20220517/StokesOutputVM_Proj.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

% Method 1: All
subplot(1, 2, 1)
plot(Table_Stokes_SA.('ESD'), Table_Stokes_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_SA.('ESD'), Table_Stokes_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Stokes Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

% Method 2: All
subplot(1, 2, 2)
plot(Table_Stokes_Proj.('ESD'), Table_Stokes_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_Proj.('ESD'), Table_Stokes_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Stokes_SA.('ESD'), Table_Stokes_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_SA{1:80, "ESD"}, Table_Stokes_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_SA{81:100, "ESD"}, Table_Stokes_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_SA{101:140, "ESD"}, Table_Stokes_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Stokes_Proj.('ESD'), Table_Stokes_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_Proj{1:80, "ESD"}, Table_Stokes_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_Proj{81:100, "ESD"}, Table_Stokes_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_Proj{101:140, "ESD"}, Table_Stokes_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Stokes_SA.('CSF'), Table_Stokes_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_SA.('CSF'), Table_Stokes_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Stokes Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Plotting all
subplot(1, 2, 2)
plot(Table_Stokes_Proj.('CSF'), Table_Stokes_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_Proj.('CSF'), Table_Stokes_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Stokes_SA.('CSF'), Table_Stokes_SA.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_SA{1:80, "CSF"}, Table_Stokes_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_SA{81:100, "CSF"}, Table_Stokes_SA{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_SA{101:140, "CSF"}, Table_Stokes_SA{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Stokes_Proj.('CSF'), Table_Stokes_Proj.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes_Proj{1:80, "CSF"}, Table_Stokes_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_Proj{81:100, "CSF"}, Table_Stokes_Proj{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_Proj{101:140, "CSF"}, Table_Stokes_Proj{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest_SA(1) = max(Table_Stokes_SA.Wt);
Highest_SA(2) = max(Table_Stokes_SA.Wt_Meas);
Highest_Proj(1) =  max(Table_Stokes_Proj.Wt);
Highest_Proj(2) = max(Table_Stokes_Proj.Wt_Meas);
MaxW_SA = max(Highest_SA);
MaxW_Proj = max(Highest_Proj);
yx_SA=linspace(0, MaxW_SA, 100);
yx_Proj=linspace(0, MaxW_Proj, 100);

% Method 1: Plot shapes separately
subplot(1, 2, 1)
plot(yx_SA, yx_SA)
hold on
plot(Table_Stokes_SA{1:80, "Wt_Meas"}, Table_Stokes_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_SA{81:100, "Wt_Meas"}, Table_Stokes_SA{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_SA{101:140, "Wt_Meas"}, Table_Stokes_SA{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Stokes Model. Using Particle Surface Area')
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
plot(Table_Stokes_Proj{1:80, "Wt_Meas"}, Table_Stokes_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_Proj{81:100, "Wt_Meas"}, Table_Stokes_Proj{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_Proj{101:140, "Wt_Meas"}, Table_Stokes_Proj{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_Proj*1.1] )
set(gca,'XLim', [0, MaxW_Proj*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_MeasVsCalc.jpg', 'Resolution', 300);

%% D1) wt against wt measured with fitted lines
% ===============================================

subplot(1, 2, 1)
plot(Table_Stokes_SA.('Wt_Meas'), Table_Stokes_SA.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_SA, yx_SA, '-k')
p=polyfit(Table_Stokes_SA.('Wt_Meas'), Table_Stokes_SA.('Wt'), 1);
px=[min(Table_Stokes_SA.('Wt_Meas')) max(Table_Stokes_SA.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Stokes_SA.("Wt_Meas")\Table_Stokes_SA.("Wt");
mx = m*Table_Stokes_SA.("Wt_Meas");
plot(Table_Stokes_SA.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Stokes Model. Using Particle Surface Area')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_SA])
set(gca, 'Xlim', [0, 1.1*MaxW_SA])
hold off

subplot(1, 2, 2)
plot(Table_Stokes_Proj.('Wt_Meas'), Table_Stokes_Proj.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx_Proj, yx_Proj, '-k')
p=polyfit(Table_Stokes_Proj.('Wt_Meas'), Table_Stokes_Proj.('Wt'), 1);
px=[min(Table_Stokes_Proj.('Wt_Meas')) max(Table_Stokes_Proj.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Stokes_Proj.("Wt_Meas")\Table_Stokes_Proj.("Wt");
mx = m*Table_Stokes_Proj.("Wt_Meas");
plot(Table_Stokes_Proj.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Stokes Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW_Proj])
set(gca, 'Xlim', [0, 1.1*MaxW_Proj])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% D2) wt against wt measured using Matlab fitlm function
% ========================================================

% Fit linear model through the intercept: SA
lm_StokesSA = fitlm(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'y~-1+x1');
m_StokesSA = lm_StokesSA.Coefficients.Estimate(1);
fitY_StokesSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_SA.Wt), max(Table_Stokes_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesSA(i) = m_StokesSA * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

% Fit linear model through the intercept: Projected area
lm_StokesProj = fitlm(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'y~-1+x1');
m_StokesProj = lm_StokesProj.Coefficients.Estimate(1);
fitY_StokesProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_Proj.Wt), max(Table_Stokes_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesProj(i) = m_StokesProj * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_MeasVsCalc_Fit.jpg', 'Resolution', 300);
%% E1) Re against Cd (ALL)
% =========================

% Method 1: Plotting all 
subplot(1, 2, 1)
plot(Table_Stokes_SA.('Re_Meas'), Table_Stokes_SA.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Stokes_SA.('Re_Calc'), Table_Stokes_SA.('Cd_Calc'), 's', ...
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
plot(Table_Stokes_Proj.('Re_Meas'), Table_Stokes_Proj.('Cd_Meas'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Stokes_Proj.('Re_Calc'), Table_Stokes_Proj.('Cd_Calc'), 's', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Cd', 'Calculated Cd', 'location', 'best')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_ReVsCd.jpg', 'Resolution', 300);

%% E2) wt against CSF (SHAPES)
% =============================

% Method 1: Shapes Plotted Separately
subplot(1, 2, 1)
plot(Table_Stokes_SA.('Re_Meas'), Table_Stokes_SA.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Stokes_SA{1:80, "Re_Calc"}, Table_Stokes_SA{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_SA{81:100, "Re_Calc"}, Table_Stokes_SA{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_SA{101:140, "Re_Calc"}, Table_Stokes_SA{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Particle Surface Area')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%set(gca, 'Xlim', [0.01, 10000])
%set(gca, 'Ylim', [0.01, 10000])
hold off

% Method 2: Shapes plotted separately
subplot(1, 2, 2)
plot(Table_Stokes_Proj.('Re_Meas'), Table_Stokes_Proj.('Cd_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7 .7 .7]')
hold on
plot(Table_Stokes_Proj{1:80, "Re_Calc"}, Table_Stokes_Proj{1:80, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes_Proj{81:100, "Re_Calc"}, Table_Stokes_Proj{81:100, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_Proj{101:140, "Re_Calc"}, Table_Stokes_Proj{101:140, "Cd_Calc"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Cd', 'Calculated Cd, Fragment', 'Calculated Cd, Fibre', ...
       'Calculated Cd, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Cd')
xlabel('Re')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/StokesVM_ReVsCd_Shapes.jpg', 'Resolution', 300);

