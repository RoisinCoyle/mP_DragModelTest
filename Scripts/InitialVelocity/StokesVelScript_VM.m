%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: StokesScript: VM
% Date created: 23.04.22
% Date last mostified: 17.05.22
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

%% Stokes' method 1
% <<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

% Set up variable arrays
Cd_Stokes = zeros(36, 10000);
CdFinal_Stokes = zeros(36, 1);
Re_Stokes = zeros(36, 10000);
ReFinal_Stokes = zeros(36, 1);
Fd_Stokes = zeros(36, 10000);
Fg_Stokes = zeros(36, 10000);
Fb_Stokes = zeros(36, 10000);
Fnet_Stokes = zeros(36, 10000);
Dist1_Stokes = zeros(36, 10000);
Acc_Stokes = zeros(36, 10000);
wtFinal_Stokes = zeros(36, 1);
DistTot_Stokes = zeros(36, 1);
FinalStep_Stokes = zeros(36, 1);
FinalTime_Stokes = zeros(36, 1);
wvel_Stokes = zeros(36, 10000);

% Begin calculation
for i=1:36
    % Set initial velocity and timestep
    wvel_Stokes(i, 1) = VM_Dataset.Initial_w(i);
    timestep = VM_Dataset.Initial_w(i)*0.1;
end

for i=1:36
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
		    CdFinal_Stokes(i) = (24.0/ReFinal_Stokes(i));
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
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end

% Store output in one array
Results_Stokes = zeros(36, 10);

for i=1:36
    Results_Stokes(i, 1) = d_equi(i);
    Results_Stokes(i, 2) = CSF(i);
    Results_Stokes(i, 3) = wtFinal_Stokes(i);
    Results_Stokes(i, 4) = wvel_Stokes(i, 1);
    Results_Stokes(i, 5) = FinalTime_Stokes(i);
    Results_Stokes(i, 6) = DistTot_Stokes(i);
    Results_Stokes(i, 7) = timestep;
    Results_Stokes(i, 8) = ReFinal_Stokes(i);
    Results_Stokes(i, 9) = CdFinal_Stokes(i);
end 

Table_Stokes_SA = array2table(Results_Stokes, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'Initial_W', ...
    'Time', 'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Meas', 'Cd_Calc'});

Table_Stokes_SA = [VM_Dataset.Shape Table_Stokes_SA];
Table_Stokes_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Stokes_SA, './DragModelsTest/Output/20220621/Stokes/StokesVelOutputVM_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Stokes_SA, './DragModelsTest/Output/20220621/Stokes/StokesVelOutputVM_SA.xls', 'WriteRowNames', true);

%% Plot the 3 plots

subplot(3, 2, 1)
plot(timesec(1:(FinalStep_Stokes(1))), wvel_Stokes(1, (1:(FinalStep_Stokes(1)))))
hold on
for n=2:6
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 3.5e-3])
lgnd = legend(sprintf('%4.2e', wvel_Stokes(1, 1)), sprintf('%4.2e', wvel_Stokes(2, 1)), ...
    sprintf('%4.2e', wvel_Stokes(3, 1)), sprintf('%4.2e', wvel_Stokes(4, 1)),...
    sprintf('%4.2e', wvel_Stokes(5, 1)), sprintf('%4.2e', wvel_Stokes(6, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 2)
plot(timesec(1:(FinalStep_Stokes(7))), wvel_Stokes(7, (1:(FinalStep_Stokes(7)))))
hold on
for n=8:12
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Stokes(7, 1)), sprintf('%4.2e', wvel_Stokes(8, 1)), ...
    sprintf('%4.2e', wvel_Stokes(9, 1)), sprintf('%4.2e', wvel_Stokes(10, 1)),...
    sprintf('%4.2e', wvel_Stokes(11, 1)), sprintf('%4.2e', wvel_Stokes(12, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 3)
plot(timesec(1:(FinalStep_Stokes(13))), wvel_Stokes(13, (1:(FinalStep_Stokes(13)))))
hold on
for n=14:18
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Stokes(13, 1)), sprintf('%4.2e', wvel_Stokes(14, 1)), ...
    sprintf('%4.2e', wvel_Stokes(15, 1)), sprintf('%4.2e', wvel_Stokes(16, 1)),...
    sprintf('%4.2e', wvel_Stokes(17, 1)), sprintf('%4.2e', wvel_Stokes(18, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 4)
plot(timesec(1:(FinalStep_Stokes(19))), wvel_Stokes(19, (1:(FinalStep_Stokes(19)))))
hold on
for n=20:24
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 3.5e-2])
lgnd = legend(sprintf('%4.2e', wvel_Stokes(19, 1)), sprintf('%4.2e', wvel_Stokes(20, 1)), ...
    sprintf('%4.2e', wvel_Stokes(21, 1)), sprintf('%4.2e', wvel_Stokes(22, 1)),...
    sprintf('%4.2e', wvel_Stokes(23, 1)), sprintf('%4.2e', wvel_Stokes(24, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 5)
plot(timesec(1:(FinalStep_Stokes(25))), wvel_Stokes(25, (1:(FinalStep_Stokes(25)))))
hold on
for n=26:30
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 2e-3])
lgnd = legend(sprintf('%4.2e', wvel_Stokes(25, 1)), sprintf('%4.2e', wvel_Stokes(26, 1)), ...
    sprintf('%4.2e', wvel_Stokes(27, 1)), sprintf('%4.2e', wvel_Stokes(28, 1)),...
    sprintf('%4.2e', wvel_Stokes(29, 1)), sprintf('%4.2e', wvel_Stokes(30, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 6)
plot(timesec(1:(FinalStep_Stokes(31))), wvel_Stokes(31, (1:(FinalStep_Stokes(31)))))
hold on
for n=32:36
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(n, (1:(FinalStep_Stokes(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_Stokes(31, 1)), sprintf('%4.2e', wvel_Stokes(32, 1)), ...
    sprintf('%4.2e', wvel_Stokes(33, 1)), sprintf('%4.2e', wvel_Stokes(34, 1)),...
    sprintf('%4.2e', wvel_Stokes(35, 1)), sprintf('%4.2e', wvel_Stokes(36, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

sgtitle(sprintf('Impact of initial velocity on calculated terminal settling velocity. \r\n Stokes (1851): Using Particle Surface Area'), 'FontWeight', 'Bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20220621/Velocity/Stokesx6.jpeg', 'Resolution', 300)