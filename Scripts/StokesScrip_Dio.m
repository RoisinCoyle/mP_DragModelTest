%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: StokesScript: Dioguardi Dataset
% Date created: 24.04.22
% Date last mostified: 23.04.22
% Purpose: To test the implementation of the Stokes drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file

% Dioguardi (2018) DOI: 10.1002/2017JB014926
% ====================================================
Dio_Dataset = readtable("SettlingVelocity calc\DioguardiSIDataSet.txt");

rho_p = table2array(Dio_Dataset(1:200, "ParticleDensity"));
rho_f = table2array(Dio_Dataset(1:200, "FluidDensity"));
vis_dyn = table2array(Dio_Dataset(1:200, "DynamicViscosity"));

d_equi = table2array(Dio_Dataset(1:200, "ParticleSize"));
size_a = table2array(Dio_Dataset(1:200, "a"));
size_b = table2array(Dio_Dataset(1:200, "b"));
size_c = table2array(Dio_Dataset(1:200, "c"));

shape_flt = table2array(Dio_Dataset(1:200, "flatness"));
shape_eln = table2array(Dio_Dataset(1:200, "elongation"));
shape_del = table2array(Dio_Dataset(1:200, "Dellino"));
shape_sph = table2array(Dio_Dataset(1:200, "Sphericity"));
shape_cir = table2array(Dio_Dataset(1:200, "Circularity"));
Reynolds = table2array(Dio_Dataset(1:200, "Re"));

wvel_meas = table2array(Dio_Dataset(:, "Wmeasured"));

% Set up and calculate additional variables:
SA_mP = zeros(200, 1);
SA_EqSph = zeros(200, 1);
Vol_mP = zeros(200, 1);
Mass_mP = zeros(200, 1);
CSF = zeros(200, 1);
rho_rel = zeros(200, 1);
vis_kin = zeros(200, 1);
ProjA_ESD = zeros(200, 1);
g=9.81;

for i=1:200
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    vis_kin(i) = vis_dyn(i) / rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Stokes' method 1
% <<<<<<<<<<<<<<<<<
% Method 1: Computing Drag force using surface area as the effective area

% Set timestep
timestep = 0.0002;

% Set initial velocity
 wvel_Stokes = zeros(200, 10000);
 wvel_Stokes(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
% Set up variable arrays
Cd_Stokes = zeros(200, 10000);
Re_Stokes = zeros(200, 10000);
Fd_Stokes = zeros(200, 10000);
Fg_Stokes = zeros(200, 10000);
Fb_Stokes = zeros(200, 10000);
Fnet_Stokes = zeros(200, 10000);
Dist1_Stokes = zeros(200, 10000);
Acc_Stokes = zeros(200, 10000);
wtFinal_Stokes = zeros(200, 1);
DistTot_Stokes = 0.0;
FinalStep_Stokes = zeros(200, 1);
FinalTime_Stokes = zeros(200, 1);
timestep_var = zeros(200, 1);

% Begin calculation
for i=1:200
    for t=1:10000

        Re_Stokes(i, t) = (rho_f(i)*wvel_Stokes(i, t)*d_equi(i))/vis_dyn(i);
        Cd_Stokes(t) = 24.0/Re_Stokes(i, t);
	
		Fd_Stokes(t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_Stokes(i,t)^2.0)*Cd_Stokes(t);
	
		Fg_Stokes(t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Stokes(t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Stokes(t) = Fg_Stokes(t) - Fb_Stokes(t) - Fd_Stokes(t);
	
		wvel_Stokes(i, t+1) = ((Fnet_Stokes(t)/Mass_mP(i))*timestep)+wvel_Stokes(i, t);
	
		Dist1_Stokes(i) = wvel_Stokes(i, t) * timestep_var(i);
		DistTot_Stokes = DistTot_Stokes + Dist1_Stokes(i);
		Acc_Stokes(t) = (wvel_Stokes(i, t+1) - wvel_Stokes(i, t))/timestep;
        
        if (Acc_Stokes(t)< 0.001)
			FinalStep_Stokes(i) = (t+1);
            FinalTime_Stokes(i) = (t+1)*timestep;
            wtFinal_Stokes(i)=wvel_Stokes(i, t+1);
            break
        end
    end
end

timesec = zeros(10000, 1);
for t=1:10000
        timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:100
    subplot(10, 10, n)
    meas = zeros(FinalStep_Stokes(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(1:(FinalStep_Stokes(n))))
     hold on
     plot(timesec(1:(FinalStep_Stokes(n))), meas(:, 1))
     hold off
end
for n = 101:200
    subplot(10, 10, n)
    meas = zeros(FinalStep_Stokes(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Stokes(n))), wvel_Stokes(1:(FinalStep_Stokes(n))))
     hold on
     plot(timesec(1:(FinalStep_Stokes(n))), meas(:, 1))
     hold off
end


% Store output in one array
Results_Stokes = zeros(200, 7);

for i=1:200
    Results_Stokes(i, 1) = d_equi(i);
    Results_Stokes(i, 2) = rho_p(i);
    Results_Stokes(i, 3) = CSF(i);
    Results_Stokes(i, 4) = timestep;
    Results_Stokes(i, 5) = wtFinal_Stokes(i);
    Results_Stokes(i, 6) = FinalTime_Stokes(i);
    Results_Stokes(i, 7) = wvel_meas(i);
end 

Table_Stokes_SA = array2table(Results_Stokes, "VariableNames", ...
    {'ESD', 'ParticleDensity', 'CSF', 'Timestep', ...
     'Wt', 'Time', 'Wt_Meas'});

writetable(Table_Stokes_SA, './DragModelsTest/Output/StokesOutputDio_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Stokes_SA, './DragModelsTest/Output/StokesOutputDio_SA.xls', 'WriteRowNames', true);


%% Stokes Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area, using Newtons Drag formula which assumes a spherical shape

% Set timestep
timestep = 0.0005;

% Set initial velocity and timestep
 wvel_Stokes2 = zeros(200, 10000);
 wvel_Stokes2(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
% Set up variable arrays
Cd_Stokes2 = zeros(200, 10000);
Re_Stokes2 = zeros(200, 10000);
Fd_Stokes2 = zeros(200, 10000);
Fg_Stokes2 = zeros(200, 10000);
Fb_Stokes2 = zeros(200, 10000);
Fnet_Stokes2 = zeros(200, 10000);
Dist1_Stokes2 = zeros(200, 10000);
Acc_Stokes2 = zeros(200, 10000);
wtFinal_Stokes2 = zeros(200, 1);
DistTot_Stokes2 = 0.0;
FinalStep_Stokes2 = zeros(200, 1);
FinalTime_Stokes2 = zeros(200, 1);

for i=1:200
    for t=1:10000
        
        Re_Stokes2(i, t) = (rho_f(i)*wvel_Stokes2(i, t)*d_equi(i))/vis_dyn(i);
        Cd_Stokes2(t) = 24.0/Re_Stokes2(i, t);
	
		Fd_Stokes2(t) = 0.5*rho_f(i)*ProjA_ESD(i)*(wvel_Stokes2(i,t)^2.0)*Cd_Stokes2(t);
	
		Fg_Stokes2(t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Stokes2(t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Stokes2(t) = Fg_Stokes2(t) - Fb_Stokes2(t) - Fd_Stokes2(t);
	
		wvel_Stokes2(i, t+1) = ((Fnet_Stokes2(t)/Mass_mP(i))*timestep)+wvel_Stokes2(i, t);
	
		Dist1_Stokes2(i) = wvel_Stokes2(i, t) * timestep;
		DistTot_Stokes2 = DistTot_Stokes2 + Dist1_Stokes2(i);
		Acc_Stokes2(t) = (wvel_Stokes2(i, t+1) - wvel_Stokes2(i, t))/timestep;
        
        if (Acc_Stokes2(t)< 0.001)
			FinalStep_Stokes2(i) = (t+1);
            FinalTime_Stokes2(i) = (t+1)*timestep;
            wtFinal_Stokes2(i)=wvel_Stokes2(i, t+1);
            break
        end
    end
end

timesec = zeros(10000, 1);
for t=1:10000
    timesec(t) = t*timestep;
end

%Check that the timestep is ok by plotting graph of w against time

for n = 1:200
    subplot(14, 10, n)
    meas = zeros(FinalStep_Stokes2(n), 1);
    meas(:, 1) = wvel_meas(n);
    plot(timesec(1:(FinalStep_Stokes2(n))), wvel_Stokes2(n, (1:(FinalStep_Stokes2(n)))))
    hold on
    plot(timesec(1:(FinalStep_Stokes2(n))), meas(:, 1))
    hold off
end

% Store output in one array
Results_Stokes_Proj = zeros(200, 5);

for i=1:200
    Results_Stokes_Proj(i, 1) = d_equi(i);
    Results_Stokes_Proj(i, 2) = CSF(i);
    Results_Stokes_Proj(i, 3) = wtFinal_Stokes2(i);
    Results_Stokes_Proj(i, 4) = FinalTime_Stokes2(i);
    Results_Stokes_Proj(i, 5) = wvel_meas(i);
end 

Table_Stokes_Proj = array2table(Results_Stokes_Proj, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', ...
    'Time', 'Wt_Meas'});

writetable(Table_Stokes_Proj, './DragModelsTest/Output/StokesOutputDio_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Stokes_Proj, './DragModelsTest/Output/StokesOutputDio_Proj.xls', 'WriteRowNames', true);

%% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:200

%% Plot Stokes output
% <<<<<<<<<<<<<<<<<<<

Table_Stokes_SA= readtable("./DragModelsTest/Output/StokesOutputDio_SA.txt", "Delimiter", ",");
Table_Stokes_Proj= readtable("./DragModelsTest/Output/StokesOutputDio_Proj.txt", "Delimiter", ",");

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
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_ESDVsW.jpg', 'Resolution', 300)

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
plot(Table_Stokes_SA{101:200, "ESD"}, Table_Stokes_SA{101:200, "Wt"}, 'og', ...
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
plot(Table_Stokes_Proj{101:200, "ESD"}, Table_Stokes_Proj{101:200, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_ESDVsW_Shapes.jpg', 'Resolution', 300)

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
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_CSFVsW.jpg', 'Resolution', 300);

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
plot(Table_Stokes_SA{101:200, "CSF"}, Table_Stokes_SA{101:200, "Wt"}, 'og', ...
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
plot(Table_Stokes_Proj{101:200, "CSF"}, Table_Stokes_Proj{101:200, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_CSFVsW_Shapes.jpg', 'Resolution', 300);

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
plot(Table_Stokes_SA{101:200, "Wt_Meas"}, Table_Stokes_SA{101:200, "Wt"}, 'og',...
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
plot(Table_Stokes_Proj{101:200, "Wt_Meas"}, Table_Stokes_Proj{101:200, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Stokes Model. Using Projected Area of Equivalent Sphere')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW_Proj*1.1] )
set(gca,'XLim', [0, MaxW_Proj*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
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
exportgraphics(gcf, './DragModelsTest/Output/StokesDio_MeasVsCalc_Eqn.jpg', 'Resolution', 300);