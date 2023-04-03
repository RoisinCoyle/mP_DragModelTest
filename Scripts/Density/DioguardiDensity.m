%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: Density assumption test: Dioguardi
% Date created: 26.06.22
% Date last mostified: 02.03.23
% Purpose: To test the model by Dioguardi satisfies the density and initial
% velocity assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
Density_Dataset = readtable("SettlingVelocity calc\DensityTestTableNew.txt");

rho_p = table2array(Density_Dataset(:, "ParticleDensity"));
rho_f = table2array(Density_Dataset(:, "FluidDensity"));
vis_dyn = table2array(Density_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(Density_Dataset(:, "KinematicViscosity"));

d_equi = table2array(Density_Dataset(:, "ParticleSize"));
size_a = table2array(Density_Dataset(:, "a"));
size_b = table2array(Density_Dataset(:, "b"));
size_c = table2array(Density_Dataset(:, "c"));
shape = table2array(Density_Dataset(:, "Shape"));

shape_flt = table2array(Density_Dataset(:, "Flatness"));
shape_eln = table2array(Density_Dataset(:, "elongation"));
shape_del = table2array(Density_Dataset(:, "Dellino"));
shape_sph = table2array(Density_Dataset(:, "Sphericity"));
shape_cir = table2array(Density_Dataset(:, "Circularity"));


% Set up and calculate additional variables:
SA_mP = zeros(54, 1);
SA_EqSph = zeros(54, 1);
Vol_mP = zeros(54, 1);
Mass_mP = zeros(54, 1);
CSF = zeros(54, 1);
rho_rel = zeros(54, 1);
ProjA_ESD = zeros(54, 1);
g=9.81;

for i=1:54
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dioguardi Method 2: Projected Area
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area in the calculation of the drag force.

Re_Dio = zeros(54, 10000);
ReFinal_Dio = zeros(54, 1);
Cd_Dio = zeros(54, 10000);
CdFinal_Dio = zeros(54, 1);
Fd_Dio = zeros(54, 10000);
Fg_Dio = zeros(54, 10000);
Fb_Dio = zeros(54, 10000);
Fnet_Dio = zeros(54, 10000);
Dist1_Dio = zeros(54, 10000);
DistTot_Dio=zeros(54, 1);
Acc_Dio = zeros(54, 10000);	
FinalTime_Dio = zeros(54, 1);
FinalStep_Dio = zeros(54, 1);
wtFinal_Dio = zeros(54, 1);

% Set timestep
timestep = 0.0002;

% Set initial velocity and timestep
wvel_Dio = zeros(54, 10000);
wvel_Dio(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
                              
% Begin calculation
for i=1:54
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

for n = 1:54
    subplot(6, 9, n)
    plot(timesec(1:(FinalStep_Dio(n))), wvel_Dio(n, (1:(FinalStep_Dio(n)))))
end

% Store output in one array
Results_Dio = zeros(54, 10);

for i=1:54
    Results_Dio(i, 1) = d_equi(i);
    Results_Dio(i, 2) = CSF(i);
    Results_Dio(i, 3) = rho_f(i);
    Results_Dio(i, 4) = rho_p(i);
    Results_Dio(i, 5) = wtFinal_Dio(i);
    Results_Dio(i, 6) = FinalTime_Dio(i);
    Results_Dio(i, 7) = DistTot_Dio(i);
    Results_Dio(i, 8) = timestep;
    Results_Dio(i, 9) = ReFinal_Dio(i);
    Results_Dio(i, 10) = CdFinal_Dio(i);
end 

Table_Dio = array2table(Results_Dio, "VariableNames", ...
    {'ESD', 'CSF', 'rho_f', 'rho_p', 'Wt', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_Dio = [Density_Dataset.Shape Table_Dio];
Table_Dio.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Dio, './DragModelsTest/Output/20220621/Density/DioguardiOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dio, './DragModelsTest/Output/20220621/Density/DioguardiOutputDensity.xls', 'WriteRowNames', true);

%% Read data in to plot
% <<<<<<<<<<<<<<<<<<<<<<<
clear
Table_Dio= readtable("./DragModelsTest/Output/20220621/Density/DioguardiOutputDensity.txt", "Delimiter", ",");

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_Dio.rho_p(i), Table_Dio.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 && i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_Dio.rho_p(i), Table_Dio.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 && i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_Dio.rho_p(i), Table_Dio.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_Dio label_t];

%% Plot boxplot

colors = {[0.6980 1 0.4118] [0.6980 1 0.4118] [0.6980 1 0.4118] ...
          [1 0.6000 0.6000] [1 0.6000 0.6000] [1 0.6000 0.6000] ...
          [0.4000 0.6980 1] [0.4000 0.6980 1] [0.4000 0.6980 1] };

fig = figure
hold on
boxplot(new_table.Wt, new_table.label_m, 'position', (1:6:54), 'widths', 5, 'boxstyle', 'outline', 'Colors', 'k')
ylim([0 0.036])
boxes = fig.Children.Children(1,1).Children(19:27)
for j = 1:length(boxes) % draw a colored patch behind each bar
        patch(boxes(j).XData, boxes(j).YData, colors{j},'FaceAlpha',.5,'EdgeAlpha',0.3);
end
ylabel('Modelled Terminal Settling Velocity (m/s)')
title(sprintf('Boxplots showing the range of modelled terminal settling velocity attained when fluid density varies from %s_f = %5.2f to %5.2f kg/m^3.', '\rho', Table_Dio.rho_f(1), Table_Dio.rho_f(6)))
subtitle('Model applied: Dioguardi et al (2018) using particle projection area as the effective area.')
set(gcf, 'WindowState', 'maximized');

exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Dio_Boxplot.jpg', 'Resolution', 1200)
%% Range table

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_Dio.Wt(i:i+5)) - min(Table_Dio.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_Dio.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_Dio.Shape(i) Table_Dio.ESD(i) Table_Dio.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/DioguardiRangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/DioguardiRangeWDensity.xls', 'WriteRowNames', true);

%% A1) wt against Density
% =====================

plot(Table_Dio{1:6, 'rho_f'}, Table_Dio{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Dio{7:12, 'rho_f'}, Table_Dio{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio{13:18, 'rho_f'}, Table_Dio{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio{19:24, 'rho_f'}, Table_Dio{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{25:30, 'rho_f'}, Table_Dio{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{31:36, 'rho_f'}, Table_Dio{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{37:42, 'rho_f'}, Table_Dio{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Dio{43:48, 'rho_f'}, Table_Dio{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Dio{49:54, 'rho_f'}, Table_Dio{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(1), Table_Dio.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(7), Table_Dio.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(13), Table_Dio.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(19), Table_Dio.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(25), Table_Dio.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(31), Table_Dio.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(37), Table_Dio.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(43), Table_Dio.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Dio.rho_p(49), Table_Dio.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title(sprintf("The impact of fluid density (%s_f) on modelled terminal settling velocity of six mP particles selected randomly from Van Melkebeke et al (2020)'s dataset.", '\rho'))
subtitle('Model applied: Dioguardi et al (2018) using particle projection area as the effective area.')
ylabel('Modelled Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Dioguardi_Density.jpg', 'Resolution', 1200)