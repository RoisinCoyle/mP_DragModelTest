%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: Density assumption test: Bagheri
% Date created: 26.06.22
% Date last mostified: 21.07.22
% Purpose: To test the model by Bagheri satisfies the density and initial
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
%% Bagheri Method 2
% <<<<<<<<<<<<<<<<<<<
% Method 2: Computing Drag force using projected area as the effective area
% in the calculation of the drag force.

% Set timestep
timestep = 0.0002;

% Set up variable arrays
Re_BB = zeros(54, 10000);
ReFinal_BB = zeros(54, 1);
CdFinal_BB = zeros(54, 1);
FormFactor_S = zeros(54, 1);
FormFactor_N = zeros(54, 1);
Correction_S = zeros(54, 1);
Correction_N = zeros(54, 1);
ratio_density = zeros(54,1);
alpha2= zeros(54,1);
beta2 = zeros(54,1);
wvel_BB = zeros(54, 10000);
Cd_BB = zeros(54, 10000);
Fd_BB = zeros(54, 10000);
Fg_BB = zeros(54, 10000);
Fb_BB = zeros(54, 10000);
Fnet_BB = zeros(54, 10000);
Dist1_BB = zeros(54, 10000);
DistTot_BB = zeros(54, 1);
Acc_BB = zeros(54, 10000);
FinalTime_BB = zeros(54, 1);
FinalStep_BB = zeros(54, 1);
wtFinal_BB = zeros(54, 1);

% Set initial velocity
wvel_BB(:, 1) = 0.0001;
for i=1:54	
    
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
for n = 1:54
    subplot(6, 9, n)
    plot(timesec(1:(FinalStep_BB(n))), wvel_BB(n, (1:(FinalStep_BB(n)))))
end
% Store output in one array
Results_BB = zeros(54, 10);

for i=1:54
    Results_BB(i, 1) = d_equi(i);
    Results_BB(i, 2) = CSF(i);
    Results_BB(i, 3) = rho_f(i);
    Results_BB(i, 4) = rho_p(i);
    Results_BB(i, 5) = wtFinal_BB(i);
    Results_BB(i, 6) = FinalTime_BB(i);
    Results_BB(i, 7) = DistTot_BB(i);
    Results_BB(i, 8) = timestep;
    Results_BB(i, 9) = ReFinal_BB(i);
    Results_BB(i, 10) = CdFinal_BB(i);
end 

Table_BB_Proj = array2table(Results_BB, "VariableNames", ...
    {'ESD', 'CSF', 'rho_f', 'rho_p', 'Wt', 'Time', ...
    'Distance', 'Timestep', ...
    'Re_Calc', 'Cd_Calc'});

Table_BB_Proj = [Density_Dataset.Shape Table_BB_Proj];
Table_BB_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Density/BagheriOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_BB_Proj, './DragModelsTest/Output/20220621/Density/BagheriOutputDensity.xls', 'WriteRowNames', true);

%% Plot Bagheri output
% <<<<<<<<<<<<<<<<<<<<<
clear
Table_BB= readtable("./DragModelsTest/Output/20220621/Density/BagheriOutputDensity.txt", "Delimiter", ",");

%% Plot boxplot

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_BB.rho_p(i), Table_BB.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 && i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_BB.rho_p(i), Table_BB.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 && i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_BB.rho_p(i), Table_BB.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_BB label_t];

boxplot(new_table.Wt, new_table.label_m)
ylabel('Terminal Settling Velocity (m/s)')
title(sprintf('Bagheri and Bonadonna (2016): Using Particle Surface area \n\r %s_f = %5.2f to %5.2f kg/m^3', '\rho', Table_BB.rho_f(1), Table_BB.rho_f(6)))
set(gcf, 'WindowState', 'maximized');

exportgraphics(gcf, './DragModelsTest/Output/20220621/Density/Bagheri_Boxplot.jpg', 'Resolution', 300)
%% Range table

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_BB.Wt(i:i+5)) - min(Table_BB.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_BB.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_BB.Shape(i) Table_BB.ESD(i) Table_BB.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/BagheriRangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/BagheriRangeWDensity.xls', 'WriteRowNames', true);

%% A1) wt against Density
% =====================

plot(Table_BB{1:6, 'rho_f'}, Table_BB{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_BB{7:12, 'rho_f'}, Table_BB{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB{13:18, 'rho_f'}, Table_BB{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB{19:24, 'rho_f'}, Table_BB{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{25:30, 'rho_f'}, Table_BB{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{31:36, 'rho_f'}, Table_BB{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{37:42, 'rho_f'}, Table_BB{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_BB{43:48, 'rho_f'}, Table_BB{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_BB{49:54, 'rho_f'}, Table_BB{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(1), Table_BB.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(7), Table_BB.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(13), Table_BB.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(19), Table_BB.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(25), Table_BB.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(31), Table_BB.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(37), Table_BB.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(43), Table_BB.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_BB.rho_p(49), Table_BB.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title('Bagheri and Bonadonna (2016): Using Projection Area.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Density/Bagheri_Density.jpg', 'Resolution', 300)