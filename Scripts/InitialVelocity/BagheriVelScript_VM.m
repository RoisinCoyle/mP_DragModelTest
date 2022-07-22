%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: BagheriScript: VM_InitialVel
% Date created: 26.05.22
% Date last mostified: 26.05.22
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

%% Bagheri Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective area
% in the calculation of the drag force.

% Set up variable arrays
Re_BB = zeros(36, 10000);
ReFinal_BB = zeros(36, 1);
CdFinal_BB = zeros(36, 1);
FormFactor_S = zeros(36, 1);
FormFactor_N = zeros(36, 1);
Correction_S = zeros(36, 1);
Correction_N = zeros(36, 1);
ratio_density = zeros(36,1);
alpha2= zeros(36,1);
beta2 = zeros(36,1);
wvel_BB = zeros(36, 10000);
Cd_BB = zeros(36, 10000);
Fd_BB = zeros(36, 10000);
Fg_BB = zeros(36, 10000);
Fb_BB = zeros(36, 10000);
Fnet_BB = zeros(36, 10000);
Dist1_BB = zeros(36, 10000);
DistTot_BB = zeros(36, 1);
Acc_BB = zeros(36, 10000);
FinalTime_BB = zeros(36, 1);
FinalStep_BB = zeros(36, 1);
wtFinal_BB = zeros(36, 1);

% Begin calculation
for i=1:36
    % Set initial velocity and timestep
    wvel_BB(i, 1) = VM_Dataset.Initial_w(i);
    timestep = VM_Dataset.Initial_w(i)*0.1;
end
     
for i=1:36	
    
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

for n = 1:36
    subplot(6, 6, n)
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end

% Store output in one array
Results_BB_10 = zeros(36, 9);

for i=1:36
    Results_BB_10(i, 1) = d_equi(i);
    Results_BB_10(i, 2) = CSF(i);
    Results_BB_10(i, 3) = wtFinal_BB(i);
    Results_BB_10(i, 4) = wvel_BB(i, 1);
    Results_BB_10(i, 5) = FinalTime_BB(i);
    Results_BB_10(i, 6) = DistTot_BB(i);
    Results_BB_10(i, 7) = timestep;
    Results_BB_10(i, 8) = ReFinal_BB(i);
    Results_BB_10(i, 9) = CdFinal_BB(i);
end 

Table_BB_Proj = array2table(Results_BB_10, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Initial_W', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_BB_Proj = [VM_Dataset.Shape Table_BB_Proj];
Table_BB_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Velocity/BagheriVelOutputVM_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Velocity/BagheriVelOutputVM_Proj.xls', 'WriteRowNames', true);

%% Plot the 3 plots

subplot(3, 2, 1)
plot(timesec(1:(FinalStep_BB(1))), wvel_BB(1, (1:(FinalStep_BB(1)))))
hold on
for n=2:6
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.007])
lgnd = legend(sprintf('%4.2e', wvel_BB(1, 1)), sprintf('%4.2e', wvel_BB(2, 1)), ...
    sprintf('%4.2e', wvel_BB(3, 1)), sprintf('%4.2e', wvel_BB(4, 1)),...
    sprintf('%4.2e', wvel_BB(5, 1)), sprintf('%4.2e', wvel_BB(6, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 2)
plot(timesec(1:(FinalStep_BB(7))), wvel_BB(7, (1:(FinalStep_BB(7)))))
hold on
for n=8:12
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Fragment, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_BB(7, 1)), sprintf('%4.2e', wvel_BB(8, 1)), ...
    sprintf('%4.2e', wvel_BB(9, 1)), sprintf('%4.2e', wvel_BB(10, 1)),...
    sprintf('%4.2e', wvel_BB(11, 1)), sprintf('%4.2e', wvel_BB(12, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 3)
plot(timesec(1:(FinalStep_BB(13))), wvel_BB(13, (1:(FinalStep_BB(13)))))
hold on
for n=14:18
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.017])
lgnd = legend(sprintf('%4.2e', wvel_BB(13, 1)), sprintf('%4.2e', wvel_BB(14, 1)), ...
    sprintf('%4.2e', wvel_BB(15, 1)), sprintf('%4.2e', wvel_BB(16, 1)),...
    sprintf('%4.2e', wvel_BB(17, 1)), sprintf('%4.2e', wvel_BB(18, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 4)
plot(timesec(1:(FinalStep_BB(19))), wvel_BB(19, (1:(FinalStep_BB(19)))))
hold on
for n=20:24
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Fibre, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_BB(19, 1)), sprintf('%4.2e', wvel_BB(20, 1)), ...
    sprintf('%4.2e', wvel_BB(21, 1)), sprintf('%4.2e', wvel_BB(22, 1)),...
    sprintf('%4.2e', wvel_BB(23, 1)), sprintf('%4.2e', wvel_BB(24, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 5)
plot(timesec(1:(FinalStep_BB(25))), wvel_BB(25, (1:(FinalStep_BB(25)))))
hold on
for n=26:30
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
lgnd = legend(sprintf('%4.2e', wvel_BB(25, 1)), sprintf('%4.2e', wvel_BB(26, 1)), ...
    sprintf('%4.2e', wvel_BB(27, 1)), sprintf('%4.2e', wvel_BB(28, 1)),...
    sprintf('%4.2e', wvel_BB(29, 1)), sprintf('%4.2e', wvel_BB(30, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

subplot(3, 2, 6)
plot(timesec(1:(FinalStep_BB(31))), wvel_BB(31, (1:(FinalStep_BB(31)))))
hold on
for n=32:36
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
hold off
title(sprintf('Film, ESD=%3.2e m, %s_{mP} = %4.1f kg/m^3', VM_Dataset.ParticleSize(n), '\rho', rho_p(n)))
ylabel(sprintf('Settling Velocity \n (m/s)'))
xlabel('Time (sec)')
set(gca, 'YLim', [0 0.012])
lgnd = legend(sprintf('%4.2e', wvel_BB(31, 1)), sprintf('%4.2e', wvel_BB(32, 1)), ...
    sprintf('%4.2e', wvel_BB(33, 1)), sprintf('%4.2e', wvel_BB(34, 1)),...
    sprintf('%4.2e', wvel_BB(35, 1)), sprintf('%4.2e', wvel_BB(36, 1)), 'Location', 'best', 'NumColumns', 2);
title(lgnd, 'Initial Velocity (m/s)', 'FontWeight', 'bold');

sgtitle(sprintf('Impact of initial velocity on calculated terminal settling velocity. \r\n Bagheri et al (2016): Using Particle Projection Area'), 'FontWeight', 'Bold');

set(gcf, 'WindowState', 'maximized')
exportgraphics(gcf, './DragModelsTest/Output/20220621/Velocity/BBx6.jpeg', 'Resolution', 300)