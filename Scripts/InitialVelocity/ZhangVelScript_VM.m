%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: ZhangScript: VM
% Date created: 23.04.22
% Date last mostified: 01.03.22
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
Zhang_EstVol = zeros(140, 1);
Zhang_deq = zeros(140, 1);
Zhang_Deqv = zeros(140, 1);
ZhangProjA = zeros(140, 1);
Zhang_EstMass = zeros(140, 1);
g=9.81;

for i=1:36
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

% Begin calculation
for i=1:36
    % Set initial velocity and timestep
    wvel_ZC(i, 1) = VM_Dataset.Initial_w(i);
    timestep = VM_Dataset.Initial_w(i)*0.1;
end

% Begin calculation
for i=1:36
    
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

for n = 1:36
    subplot(6, 6, n)
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end

% Store output in one array
Results_ZC = zeros(36, 9);

for i=1:36
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = wtFinal_ZC(i);
    Results_ZC(i, 4) = wvel_ZC(i, 1);
    Results_ZC(i, 5) = FinalTime_ZC(i);
    Results_ZC(i, 6) = DistTot_ZC(i);
    Results_ZC(i, 7) = timestep;
    Results_ZC(i, 8) = ReFinal_ZC(i);
    Results_ZC(i, 9) = CdFinal_ZC(i);
end 

Table_ZC_SA = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Initial_W', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_ZC_SA = [VM_Dataset.Shape Table_ZC_SA];
Table_ZC_SA.Properties.VariableNames(1) = {'Shape'};

% writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Velocity/ZhangVelOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Velocity/ZhangOutputVM_SA.xls', 'WriteRowNames', true);

%% Read in data
Table_ZC_SA = readtable('./DragModelsTest/Output/20220621/Velocity/ZhangVelOutputVM_SA.txt', 'Delimiter', ',');

%% Plot the 3 plots

subplot(3, 2, 1)
plot(timesec(1:(FinalStep_ZC(1))), wvel_ZC(1, (1:(FinalStep_ZC(1)))))
hold on
for n=2:6
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.0035])
lgnd = legend(sprintf('%4.2e', wvel_ZC(1, 1)), sprintf('%4.2e', wvel_ZC(2, 1)), ...
    sprintf('%4.2e', wvel_ZC(3, 1)), sprintf('%4.2e', wvel_ZC(4, 1)),...
    sprintf('%4.2e', wvel_ZC(5, 1)), sprintf('%4.2e', wvel_ZC(6, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 2)
plot(timesec(1:(FinalStep_ZC(7))), wvel_ZC(7, (1:(FinalStep_ZC(7)))))
hold on
for n=8:12
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.009])
lgnd = legend(sprintf('%4.2e', wvel_ZC(7, 1)), sprintf('%4.2e', wvel_ZC(8, 1)), ...
    sprintf('%4.2e', wvel_ZC(9, 1)), sprintf('%4.2e', wvel_ZC(10, 1)),...
    sprintf('%4.2e', wvel_ZC(11, 1)), sprintf('%4.2e', wvel_ZC(12, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 3)
plot(timesec(1:(FinalStep_ZC(13))), wvel_ZC(13, (1:(FinalStep_ZC(13)))))
hold on
for n=14:18
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.0125])
lgnd = legend(sprintf('%4.2e', wvel_ZC(13, 1)), sprintf('%4.2e', wvel_ZC(14, 1)), ...
    sprintf('%4.2e', wvel_ZC(15, 1)), sprintf('%4.2e', wvel_ZC(16, 1)),...
    sprintf('%4.2e', wvel_ZC(17, 1)), sprintf('%4.2e', wvel_ZC(18, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 4)
plot(timesec(1:(FinalStep_ZC(19))), wvel_ZC(19, (1:(FinalStep_ZC(19)))))
hold on
for n=20:24
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_ZC(19, 1)), sprintf('%4.2e', wvel_ZC(20, 1)), ...
    sprintf('%4.2e', wvel_ZC(21, 1)), sprintf('%4.2e', wvel_ZC(22, 1)),...
    sprintf('%4.2e', wvel_ZC(23, 1)), sprintf('%4.2e', wvel_ZC(24, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 5)
plot(timesec(1:(FinalStep_ZC(25))), wvel_ZC(25, (1:(FinalStep_ZC(25)))))
hold on
for n=26:30
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.007])
lgnd = legend(sprintf('%4.2e', wvel_ZC(25, 1)), sprintf('%4.2e', wvel_ZC(26, 1)), ...
    sprintf('%4.2e', wvel_ZC(27, 1)), sprintf('%4.2e', wvel_ZC(28, 1)),...
    sprintf('%4.2e', wvel_ZC(29, 1)), sprintf('%4.2e', wvel_ZC(30, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 6)
plot(timesec(1:(FinalStep_ZC(31))), wvel_ZC(31, (1:(FinalStep_ZC(31)))))
hold on
for n=32:36
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.008])
lgnd = legend(sprintf('%4.2e', wvel_ZC(31, 1)), sprintf('%4.2e', wvel_ZC(32, 1)), ...
    sprintf('%4.2e', wvel_ZC(33, 1)), sprintf('%4.2e', wvel_ZC(34, 1)),...
    sprintf('%4.2e', wvel_ZC(35, 1)), sprintf('%4.2e', wvel_ZC(36, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');
sgtitle(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Zhang and Choi (2021) using particle surface area as the effective area.'), 'FontWeight', 'Bold');


set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230301/Velocity/ZCx6SA.jpeg', 'Resolution', 1200)

%% New plot:
% All particles on one plot, initial velocity on x axis (log) and terminal
% settling velocity on the y axis (not log).

PlotColor = {[0 0 1] [0 0 1] [1 0 0] [1 0 0] [0 1 0] [0 1 0]};

PlotColor = num2cell(PlotColor, 3);

ColorOrder_RC = [0 0 1; 1 0 0; 0 1 0];

FaceOrder_RC = [0 0 1; 0 0 1; 1 0 0; 1 0 0; 0 1 0; 0 1 0];

LineStyleOrder_RC = ["-o"; "-^"];

for i = 1:12:36
    plot(Table_ZC_SA.Initial_W(i:i+5), Table_ZC_SA.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
for i = 7:12:36
    plot(Table_ZC_SA.Initial_W(i:i+5), Table_ZC_SA.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
hold off
title(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Zhang and Choi (2021) using particle surface area as the effective area.'), 'FontWeight', 'Bold');
ylabel(sprintf('Modelled Terminal Settling \n Velocity (m/s)'))
xlabel('Initial Settling Velocity specified (m/s)')
ax=gca;
ax.ColorOrder = ColorOrder_RC;
ax.LineStyleOrder = LineStyleOrder_RC;
set(gca, 'XScale', 'log')
set (gca, 'Xlim', [3e-6 2e-3])
lgnd = legend(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(1), '\rho', VM_Dataset.ParticleDensity(1)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(13), '\rho', VM_Dataset.ParticleDensity(13)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(25), '\rho', VM_Dataset.ParticleDensity(25)), ...
    sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(7), '\rho', VM_Dataset.ParticleDensity(7)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(19), '\rho', VM_Dataset.ParticleDensity(19)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_SA.ESD(31), '\rho', VM_Dataset.ParticleDensity(31)), ...
    'Location', 'southoutside', 'NumColumns', 2);
title(lgnd, 'Particle Properties', 'FontWeight', 'bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230403/Velocity/ZC_SAxAll_Terminal_Init.jpeg', 'Resolution', 1200)

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
for i=1:36
    % Set initial velocity and timestep
    wvel_ZC(i, 1) = VM_Dataset.Initial_w(i);
    timestep = VM_Dataset.Initial_w(i)*0.1;
end

% Begin calculation
for i=1:36
    
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
for n = 1:36
    subplot(6, 6, n)
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end

% Store output in one array
Results_ZC = zeros(36, 9);

for i=1:36
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = wtFinal_ZC(i);
    Results_ZC(i, 4) = wvel_ZC(i, 1);
    Results_ZC(i, 5) = FinalTime_ZC(i);
    Results_ZC(i, 6) = DistTot_ZC(i);
    Results_ZC(i, 7) = timestep;
    Results_ZC(i, 8) = ReFinal_ZC(i);
    Results_ZC(i, 9) = CdFinal_ZC(i);
end 

Table_ZC_Proj = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Initial_W', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_ZC_Proj = [VM_Dataset.Shape Table_ZC_Proj];
Table_ZC_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Velocity/ZhangVelOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Velocity/ZhangVelOutputVM_Proj.xls', 'WriteRowNames', true);

%% Read in data
%clear
Table_ZC_Proj = readtable('./DragModelsTest/Output/20220621/Velocity/ZhangVelOutputVM_Proj.txt', 'Delimiter', ',');

%% Plot the 3 plots

subplot(3, 2, 1)
plot(timesec(1:(FinalStep_ZC(1))), wvel_ZC(1, (1:(FinalStep_ZC(1)))))
hold on
for n=2:6
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_ZC(1, 1)), sprintf('%4.2e', wvel_ZC(2, 1)), ...
    sprintf('%4.2e', wvel_ZC(3, 1)), sprintf('%4.2e', wvel_ZC(4, 1)),...
    sprintf('%4.2e', wvel_ZC(5, 1)), sprintf('%4.2e', wvel_ZC(6, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 2)
plot(timesec(1:(FinalStep_ZC(7))), wvel_ZC(7, (1:(FinalStep_ZC(7)))))
hold on
for n=8:12
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.015])
lgnd = legend(sprintf('%4.2e', wvel_ZC(7, 1)), sprintf('%4.2e', wvel_ZC(8, 1)), ...
    sprintf('%4.2e', wvel_ZC(9, 1)), sprintf('%4.2e', wvel_ZC(10, 1)),...
    sprintf('%4.2e', wvel_ZC(11, 1)), sprintf('%4.2e', wvel_ZC(12, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 3)
plot(timesec(1:(FinalStep_ZC(13))), wvel_ZC(13, (1:(FinalStep_ZC(13)))))
hold on
for n=14:18
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.017])
lgnd = legend(sprintf('%4.2e', wvel_ZC(13, 1)), sprintf('%4.2e', wvel_ZC(14, 1)), ...
    sprintf('%4.2e', wvel_ZC(15, 1)), sprintf('%4.2e', wvel_ZC(16, 1)),...
    sprintf('%4.2e', wvel_ZC(17, 1)), sprintf('%4.2e', wvel_ZC(18, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 4)
plot(timesec(1:(FinalStep_ZC(19))), wvel_ZC(19, (1:(FinalStep_ZC(19)))))
hold on
for n=20:24
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'Ylim', [0 0.025])
lgnd = legend(sprintf('%4.2e', wvel_ZC(19, 1)), sprintf('%4.2e', wvel_ZC(20, 1)), ...
    sprintf('%4.2e', wvel_ZC(21, 1)), sprintf('%4.2e', wvel_ZC(22, 1)),...
    sprintf('%4.2e', wvel_ZC(23, 1)), sprintf('%4.2e', wvel_ZC(24, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 5)
plot(timesec(1:(FinalStep_ZC(25))), wvel_ZC(25, (1:(FinalStep_ZC(25)))))
hold on
for n=26:30
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'Ylim', [0 0.01])
lgnd = legend(sprintf('%4.2e', wvel_ZC(25, 1)), sprintf('%4.2e', wvel_ZC(26, 1)), ...
    sprintf('%4.2e', wvel_ZC(27, 1)), sprintf('%4.2e', wvel_ZC(28, 1)),...
    sprintf('%4.2e', wvel_ZC(29, 1)), sprintf('%4.2e', wvel_ZC(30, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 6)
plot(timesec(1:(FinalStep_ZC(31))), wvel_ZC(31, (1:(FinalStep_ZC(31)))))
hold on
for n=32:36
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Modelled Settling \n Velocity (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.012])
lgnd = legend(sprintf('%4.2e', wvel_ZC(31, 1)), sprintf('%4.2e', wvel_ZC(32, 1)), ...
    sprintf('%4.2e', wvel_ZC(33, 1)), sprintf('%4.2e', wvel_ZC(34, 1)),...
    sprintf('%4.2e', wvel_ZC(35, 1)), sprintf('%4.2e', wvel_ZC(36, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');
sgtitle(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Zhang and Choi (2021) using particle projection area as the effective area.'), 'FontWeight', 'Bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230301/Velocity/ZCx6Proj.jpeg', 'Resolution', 1200)

%% New plot:
% All particles on one plot, initial velocity on x axis (log) and terminal
% settling velocity on the y axis (not log).

PlotColor = {[0 0 1] [0 0 1] [1 0 0] [1 0 0] [0 1 0] [0 1 0]};

PlotColor = num2cell(PlotColor, 3);

ColorOrder_RC = [0 0 1; 1 0 0; 0 1 0];

FaceOrder_RC = [0 0 1; 0 0 1; 1 0 0; 1 0 0; 0 1 0; 0 1 0];

LineStyleOrder_RC = ["-o"; "-^"];

for i = 1:12:36
    plot(Table_ZC_Proj.Initial_W(i:i+5), Table_ZC_Proj.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
for i = 7:12:36
    plot(Table_ZC_Proj.Initial_W(i:i+5), Table_ZC_Proj.Wt(i:i+5),...
        'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', FaceOrder_RC(((i+5)/6), :))
    hold on
end
hold off
title(sprintf('Graphs demonstrating that the specified initial velocity has negligible impact on the modelled terminal settling velocity. \r\n Model applied: Zhang and Choi (2021) using particle projection area as the effective area.'), 'FontWeight', 'Bold');
ylabel(sprintf('Modelled Terminal Settling \n Velocity (m/s)'))
xlabel('Initial Settling Velocity specified (m/s)')
ax=gca;
ax.ColorOrder = ColorOrder_RC;
ax.LineStyleOrder = LineStyleOrder_RC;
set(gca, 'XScale', 'log')
set (gca, 'Xlim', [3e-6 2e-3])
lgnd = legend(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(1), '\rho', VM_Dataset.ParticleDensity(1)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(13), '\rho', VM_Dataset.ParticleDensity(13)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(25), '\rho', VM_Dataset.ParticleDensity(25)), ...
    sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(7), '\rho', VM_Dataset.ParticleDensity(7)), ...
    sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(19), '\rho', VM_Dataset.ParticleDensity(19)), ...
    sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', Table_ZC_Proj.ESD(31), '\rho', VM_Dataset.ParticleDensity(31)), ...
    'Location', 'southoutside', 'NumColumns', 2);
title(lgnd, 'Particle Properties', 'FontWeight', 'bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20230403/Velocity/ZC_ProjxAll_Terminal_Init.jpeg', 'Resolution', 1200)
