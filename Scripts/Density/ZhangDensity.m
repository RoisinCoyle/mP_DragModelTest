%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: Density assumption test: Zhang
% Date created: 26.06.22
% Date last mostified: 02.03.23
% Purpose: To test the model by Zhang satisfies the density and initial
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
Re_ZC = zeros(54, 10000);
ReFinal_ZC = zeros(54, 1);
shape_ASF = zeros(54, 1);
wvel_ZC = zeros(54, 10000);
Cd_ZC = zeros(54, 10000);
CdFinal_ZC = zeros(54, 1);
Fd_ZC = zeros(54, 10000);
Fg_ZC = zeros(54, 10000);
Fb_ZC = zeros(54, 10000);
Fnet_ZC = zeros(54, 10000);
Dist1_ZC = zeros(54, 10000);
DistTot_ZC = zeros(54, 1);
Acc_ZC = zeros(54, 10000);
FinalTime_ZC = zeros(54, 1);
FinalStep_ZC = zeros(54, 1);
wtFinal_ZC = zeros(54, 1);

% Set initial velocity and timestep
wvel_ZC(:, 1) = 0.0001; % Note that earlier tests have shown that the
                        % terminal velocity is independent of the
                        % initial velocity.

% Begin calculation
for i=1:54	
    
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

for n = 1:54
    subplot(6, 9, n)
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end

% Store output in one array
Results_ZC = zeros(54, 10);

for i=1:54
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = rho_f(i);
    Results_ZC(i, 4) = rho_p(i);
    Results_ZC(i, 5) = wtFinal_ZC(i);
    Results_ZC(i, 6) = FinalTime_ZC(i);
    Results_ZC(i, 7) = DistTot_ZC(i);
    Results_ZC(i, 8) = timestep;
    Results_ZC(i, 9) = ReFinal_ZC(i);
    Results_ZC(i, 10) = CdFinal_ZC(i);
end 

Table_ZC_SA = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'rho_f', 'rho_p', 'Wt', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_ZC_SA = [Density_Dataset.Shape Table_ZC_SA];
Table_ZC_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Density/ZhangOutputDensity_SA.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_SA, './DragModelsTest/Output/20220621/Density/ZhangOutputDensity_SA.xls', 'WriteRowNames', true);

%% Zhang Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective
% area, using Newtons Drag formula

% Set timestep
timestep = 0.00015;

% Set up variable arrays
Re_ZC = zeros(54, 10000);
shape_ASF = zeros(54, 1);
wvel_ZC = zeros(54, 10000);
Cd_ZC = zeros(54, 10000);
Fd_ZC = zeros(54, 10000);
Fg_ZC = zeros(54, 10000);
Fb_ZC = zeros(54, 10000);
Fnet_ZC = zeros(54, 10000);
Dist1_ZC = zeros(54, 10000);
DistTot_ZC = zeros(54, 1);
Acc_ZC = zeros(54, 10000);
FinalTime_ZC = zeros(54, 1);
FinalStep_ZC = zeros(54, 1);
wtFinal_ZC = zeros(54, 1);

% Set initial velocity and timestep
wvel_ZC(:, 1) = 0.0001; % Note that earlier tests have shown that the
                        % terminal velocity is independent of the
                        % initial velocity.

% Begin calculation
for i=1:54	
    
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

for n = 1:54
    subplot(6, 9, n)
    plot(timesec(1:(FinalStep_ZC(n))), wvel_ZC(n, (1:(FinalStep_ZC(n)))))
end

% Store output in one array
Results_ZC = zeros(54, 10);

for i=1:54
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = rho_f(i);
    Results_ZC(i, 4) = rho_p(i);
    Results_ZC(i, 5) = wtFinal_ZC(i);
    Results_ZC(i, 6) = FinalTime_ZC(i);
    Results_ZC(i, 7) = DistTot_ZC(i);
    Results_ZC(i, 8) = timestep;
    Results_ZC(i, 9) = ReFinal_ZC(i);
    Results_ZC(i, 10) = CdFinal_ZC(i);
end 

Table_ZC_Proj = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'rho_f', 'rho_p', 'Wt', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_ZC_Proj = [Density_Dataset.Shape Table_ZC_Proj];
Table_ZC_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Density/ZhangOutputDensity_Proj.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_ZC_Proj, './DragModelsTest/Output/20220621/Density/ZhangOutputDensity_Proj.xls', 'WriteRowNames', true);
%% Read in data
% <<<<<<<<<<<<<<<<<<<<<
clear
Table_ZC_SA = readtable("./DragModelsTest/Output/20220621/Density/ZhangOutputDensity_SA.txt", "Delimiter", ",");
Table_ZC_Proj = readtable("./DragModelsTest/Output/20220621/Density/ZhangOutputDensity_Proj.txt", "Delimiter", ",");

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_ZC_SA.rho_p(i), Table_ZC_SA.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 & i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_ZC_SA.rho_p(i), Table_ZC_SA.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 & i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_ZC_SA.rho_p(i), Table_ZC_SA.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_ZC_SA label_t];

%% Plot boxplot: SA
colors = {[0.6980 1 0.4118] [0.6980 1 0.4118] [0.6980 1 0.4118] ...
          [1 0.6000 0.6000] [1 0.6000 0.6000] [1 0.6000 0.6000] ...
          [0.4000 0.6980 1] [0.4000 0.6980 1] [0.4000 0.6980 1] };

fig = figure
hold on
boxplot(new_table.Wt, new_table.label_m, 'position', (1:6:54), 'widths', 5, 'boxstyle', 'outline', 'Colors', 'k')
ylim([0 0.05])
boxes = fig.Children.Children(1,1).Children(19:27)
for j = 1:length(boxes) % draw a colored patch behind each bar
        patch(boxes(j).XData, boxes(j).YData, colors{j},'FaceAlpha',.5,'EdgeAlpha',0.3);
end
ylabel('Modelled Terminal Settling Velocity (m/s)')
title(sprintf('Boxplots showing the range of modelled terminal settling velocity attained when fluid density varies from %s_f = %5.2f to %5.2f kg/m^3.', '\rho', Table_ZC_SA.rho_f(1), Table_ZC_SA.rho_f(6)))
subtitle('Model applied: Zhang and Choi (2021) using particle surface area as the effective area.')

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/ZhangSA_Boxplot.jpg', 'Resolution', 1200)

%% Range table: SA

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_ZC_SA.Wt(i:i+5)) - min(Table_ZC_SA.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_ZC_SA.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_ZC_SA.Shape(i) Table_ZC_SA.ESD(i) Table_ZC_SA.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/ZhangSARangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/ZhangSARangeWDensity.xls', 'WriteRowNames', true);

%% Read data: Proj

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_ZC_Proj.rho_p(i), Table_ZC_Proj.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 & i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_ZC_Proj.rho_p(i), Table_ZC_Proj.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 & i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_ZC_Proj.rho_p(i), Table_ZC_Proj.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_ZC_Proj label_t];

%% Plot boxplot: Proj
colors = {[0.6980 1 0.4118] [0.6980 1 0.4118] [0.6980 1 0.4118] ...
          [1 0.6000 0.6000] [1 0.6000 0.6000] [1 0.6000 0.6000] ...
          [0.4000 0.6980 1] [0.4000 0.6980 1] [0.4000 0.6980 1] };

fig = figure
hold on
boxplot(new_table.Wt, new_table.label_m, 'position', (1:6:54), 'widths', 5, 'boxstyle', 'outline', 'Colors', 'k')
ylim([0 0.08])
boxes = fig.Children.Children(1,1).Children(19:27)
for j = 1:length(boxes) % draw a colored patch behind each bar
        patch(boxes(j).XData, boxes(j).YData, colors{j},'FaceAlpha',.5,'EdgeAlpha',0.3);
end
ylabel('Modelled Terminal Settling Velocity (m/s)')
title(sprintf('Boxplots showing the range of modelled terminal settling velocity attained when fluid density varies from %s_f = %5.2f to %5.2f kg/m^3.', '\rho', Table_ZC_Proj.rho_f(1), Table_ZC_Proj.rho_f(6)))
subtitle('Model applied: Zhang and Choi (2021) using particle projection area as the effective area.')

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/ZhangProj_Boxplot.jpg', 'Resolution', 1200)
%% Range table: Projected Area

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_ZC_Proj.Wt(i:i+5)) - min(Table_ZC_Proj.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_ZC_Proj.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_ZC_Proj.Shape(i) Table_ZC_Proj.ESD(i) Table_ZC_Proj.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/ZhangProjRangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/ZhangProjRangeWDensity.xls', 'WriteRowNames', true);

%% A1) wt against Density: SA
% ============================

plot(Table_ZC_SA{1:6, 'rho_f'}, Table_ZC_SA{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_ZC_SA{7:12, 'rho_f'}, Table_ZC_SA{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC_SA{13:18, 'rho_f'}, Table_ZC_SA{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC_SA{19:24, 'rho_f'}, Table_ZC_SA{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_SA{25:30, 'rho_f'}, Table_ZC_SA{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_SA{31:36, 'rho_f'}, Table_ZC_SA{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_SA{37:42, 'rho_f'}, Table_ZC_SA{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_ZC_SA{43:48, 'rho_f'}, Table_ZC_SA{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_ZC_SA{49:54, 'rho_f'}, Table_ZC_SA{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(1), Table_ZC_SA.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(7), Table_ZC_SA.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(13), Table_ZC_SA.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(19), Table_ZC_SA.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(25), Table_ZC_SA.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(31), Table_ZC_SA.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(37), Table_ZC_SA.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(43), Table_ZC_SA.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(49), Table_ZC_SA.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title(sprintf("The impact of fluid density (%s_f) on modelled terminal settling velocity of six mP particles selected randomly from Van Melkebeke et al (2020)'s dataset.", '\rho'))
subtitle('Model applied: Zhang and Choi (2021) using particle surface area as the effective area.')
ylabel('Modelled Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Zhang_DensitySA.jpg', 'Resolution', 1200)

%% A2) wt against Density: Proj
% ============================

plot(Table_ZC_Proj{1:6, 'rho_f'}, Table_ZC_Proj{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_ZC_Proj{7:12, 'rho_f'}, Table_ZC_Proj{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC_Proj{13:18, 'rho_f'}, Table_ZC_Proj{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC_Proj{19:24, 'rho_f'}, Table_ZC_Proj{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_Proj{25:30, 'rho_f'}, Table_ZC_Proj{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_Proj{31:36, 'rho_f'}, Table_ZC_Proj{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC_Proj{37:42, 'rho_f'}, Table_ZC_Proj{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_ZC_Proj{43:48, 'rho_f'}, Table_ZC_Proj{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_ZC_Proj{49:54, 'rho_f'}, Table_ZC_Proj{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(1), Table_ZC_SA.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(7), Table_ZC_SA.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(13), Table_ZC_SA.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(19), Table_ZC_SA.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(25), Table_ZC_SA.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(31), Table_ZC_SA.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(37), Table_ZC_SA.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(43), Table_ZC_SA.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_ZC_SA.rho_p(49), Table_ZC_SA.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title(sprintf("The impact of fluid density (%s_f) on modelled terminal settling velocity of six mP particles selected randomly from Van Melkebeke et al (2020)'s dataset.", '\rho'))
subtitle('Model applied: Zhang and Choi (2021) using particle projection area as the effective area.')
ylabel('Modelled Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Zhang_DensityProj.jpg', 'Resolution', 1200)