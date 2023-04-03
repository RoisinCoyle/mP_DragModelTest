%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DioguardiScript: VM
% Date created: 23.04.22
% Date last mostified: 01.03.23
% Purpose: To test the impact of initial velocity on the terminal velocity
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
VM_Dataset = readtable("SettlingVelocity calc\VelocityTestTableNew.txt");

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
Powers = table2array(VM_Dataset(:, "Powers"));


% Set up and calculate additional variables:
SA_mP = zeros(36, 1);
SA_EqSph = zeros(36, 1);
Vol_mP = zeros(36, 1);
Mass_mP = zeros(36, 1);
CSF = zeros(36, 1);
rho_rel = zeros(36, 1);
ProjA_ESD = zeros(36, 1);
g=9.81;

for i=1:36
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dioguardi Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area in the calculation of the draf force.

Re_Dio = zeros(36, 10000);
ReFinal_Dio = zeros(36, 1);
Cd_Dio = zeros(36, 10000);
CdFinal_Dio = zeros(36, 1);
Fd_Dio = zeros(36, 10000);
Fg_Dio = zeros(36, 10000);
Fb_Dio = zeros(36, 10000);
Fnet_Dio = zeros(36, 10000);
Dist1_Dio = zeros(36, 10000);
DistTot_Dio=zeros(36, 1);
Acc_Dio = zeros(36, 10000);	
FinalTime_Dio = zeros(36, 1);
FinalStep_Dio = zeros(36, 1);
wtFinal_Dio = zeros(36, 1);
wvel_Dio = zeros(36, 10000);

% Begin calculation
for i=1:36
    % Set initial velocity and timestep
    wvel_Dio(i, 1) = VM_Dataset.Initial_w(i);
    timestep = VM_Dataset.Initial_w(i)*0.1;
end
                              
% Begin calculation
for i=1:36
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

for n = 1:36
    subplot(6, 6, n)
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end

% Store output in one array
Results_Dio = zeros(36, 9);

for i=1:36
    Results_Dio(i, 1) = d_equi(i);
    Results_Dio(i, 2) = CSF(i);
    Results_Dio(i, 3) = wtFinal_Dio(i);
    Results_Dio(i, 4) = wvel_Dio(i, 1);
    Results_Dio(i, 5) = FinalTime_Dio(i);
    Results_Dio(i, 6) = DistTot_Dio(i);
    Results_Dio(i, 7) = timestep;
    Results_Dio(i, 8) = ReFinal_Dio(i);
    Results_Dio(i, 9) = CdFinal_Dio(i);
end 

Table_Dio_Proj = array2table(Results_Dio, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Initial_W', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Results_Dioguardi = zeros(36, 5);

Table_Dio_Proj = [VM_Dataset.Shape Table_Dio_Proj];
Table_Dio_Proj.Properties.VariableNames(1) = {'Shape'};

% writetable(Table_Dio_Proj, './DragModelsTest/Output/20220621/Velocity/DioguardiVelOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dio_Proj, './DragModelsTest/Output/20220621/Velocity/DioguardiVelOutputVM_Proj.xls', 'WriteRowNames', true);

%% Read in data

Table_Dio_Proj = readtable('./DragModelsTest/Output/20220621/Velocity/DioguardiVelOutputVM_Proj.txt', 'Delimiter', ',')

%% Plot the 6 plots

subplot(3, 2, 1)
plot(timesec(1:(FinalStep_Dio(1))), wvel_Dio(1, (1:(FinalStep_Dio(1)))))
hold on
for n=2:6
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Dio(1, 1)), sprintf('%4.2e', wvel_Dio(2, 1)), ...
    sprintf('%4.2e', wvel_Dio(3, 1)), sprintf('%4.2e', wvel_Dio(4, 1)),...
    sprintf('%4.2e', wvel_Dio(5, 1)), sprintf('%4.2e', wvel_Dio(6, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 2)
plot(timesec(1:(FinalStep_Dio(7))), wvel_Dio(7, (1:(FinalStep_Dio(7)))))
hold on
for n=8:12
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.012])
lgnd = legend(sprintf('%4.2e', wvel_Dio(7, 1)), sprintf('%4.2e', wvel_Dio(8, 1)), ...
    sprintf('%4.2e', wvel_Dio(9, 1)), sprintf('%4.2e', wvel_Dio(10, 1)),...
    sprintf('%4.2e', wvel_Dio(11, 1)), sprintf('%4.2e', wvel_Dio(12, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 3)
plot(timesec(1:(FinalStep_Dio(13))), wvel_Dio(13, (1:(FinalStep_Dio(13)))))
hold on
for n=14:18
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Dio(13, 1)), sprintf('%4.2e', wvel_Dio(14, 1)), ...
    sprintf('%4.2e', wvel_Dio(15, 1)), sprintf('%4.2e', wvel_Dio(16, 1)),...
    sprintf('%4.2e', wvel_Dio(17, 1)), sprintf('%4.2e', wvel_Dio(18, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 4)
plot(timesec(1:(FinalStep_Dio(19))), wvel_Dio(19, (1:(FinalStep_Dio(19)))))
hold on
for n=20:24
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Dio(19, 1)), sprintf('%4.2e', wvel_Dio(20, 1)), ...
    sprintf('%4.2e', wvel_Dio(21, 1)), sprintf('%4.2e', wvel_Dio(22, 1)),...
    sprintf('%4.2e', wvel_Dio(23, 1)), sprintf('%4.2e', wvel_Dio(24, 1)), 'Location', 'east', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 5)
plot(timesec(1:(FinalStep_Dio(25))), wvel_Dio(25, (1:(FinalStep_Dio(25)))))
hold on
for n=26:30
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Dio(25, 1)), sprintf('%4.2e', wvel_Dio(26, 1)), ...
    sprintf('%4.2e', wvel_Dio(27, 1)), sprintf('%4.2e', wvel_Dio(28, 1)),...
    sprintf('%4.2e', wvel_Dio(29, 1)), sprintf('%4.2e', wvel_Dio(30, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 6)
plot(timesec(1:(FinalStep_Dio(31))), wvel_Dio(31, (1:(FinalStep_Dio(31)))))
hold on
for n=32:36
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.012])
lgnd = legend(sprintf('%4.2e', wvel_Dio(31, 1)), sprintf('%4.2e', wvel_Dio(32, 1)), ...
    sprintf('%4.2e', wvel_Dio(33, 1)), sprintf('%4.2e', wvel_Dio(34, 1)),...
    sprintf('%4.2e', wvel_Dio(35, 1)), sprintf('%4.2e', wvel_Dio(36, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

sgtitle(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Dioguardi et al (2018) using particle projection area as the effective area.'), 'FontWeight', 'Bold');


set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230301/Velocity/Diox6.jpeg', 'Resolution', 1200)

%% New plot:
% All particles on one plot, initial velocity on x axis (log) and terminal
% settling velocity on the y axis (not log).

PlotColor = {[0 0 1] [0 0 1] [1 0 0] [1 0 0] [0 1 0] [0 1 0]};

PlotColor = num2cell(PlotColor, 3);

ColorOrder_RC = [0 0 1; 1 0 0; 0 1 0];

FaceOrder_RC = [0 0 1; 0 0 1; 1 0 0; 1 0 0; 0 1 0; 0 1 0];

LineStyleOrder_RC = ["-o"; "-^"];

for i = 1:12:36
    plot(Table_Dio_Proj.Initial_W(i:i+5), Table_Dio_Proj.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
for i = 7:12:36
    plot(Table_Dio_Proj.Initial_W(i:i+5), Table_Dio_Proj.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
hold off
title(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Dioguardi et al (2018) using particle projection area as the effective area.'), 'FontWeight', 'Bold');
ylabel(sprintf('Modelled Terminal Settling \n Velocity (m/s)'))
xlabel('Initial Settling Velocity specified (m/s)')
ax=gca;
ax.ColorOrder = ColorOrder_RC;
ax.LineStyleOrder = LineStyleOrder_RC;
set(gca, 'XScale', 'log')
set (gca, 'Xlim', [3e-6 2e-3])
lgnd = legend(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(1), '\rho', VM_Dataset.ParticleDensity(1)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(13), '\rho', VM_Dataset.ParticleDensity(13)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(25), '\rho', VM_Dataset.ParticleDensity(25)), ...
    sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(7), '\rho', VM_Dataset.ParticleDensity(7)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(19), '\rho', VM_Dataset.ParticleDensity(19)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_Dio_Proj.ESD(31), '\rho', VM_Dataset.ParticleDensity(31)), ...
    'Location', 'southoutside', 'NumColumns', 2);
title(lgnd, 'Particle Properties', 'FontWeight', 'bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230403/Velocity/DioxAll_Terminal_Init.jpeg', 'Resolution', 1200)
