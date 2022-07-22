%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DragModels
% Date created: 21.03.22
% Date last mostified: 21.04.22
% Purpose: To test the implementation of the drag models on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Models tested:
%               Stokes Cd
%               Dietrich Wt
%               Dioguardi et al. (2018) Cd
%               Bagheri and Bonadonna (2019) Cd
%               Francalanci et al. (2021) Wt
%               Zhang and Choi (2021) Cd
%               Yu et al. (2022) Cd and Wt

%% Read in data file

Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

rho_p = table2array(Dataset(:, "ParticleDensity"));
rho_f = table2array(Dataset(:, "FluidDensity"));
vis_dyn = table2array(Dataset(:, "DynamicViscosity"));
vis_kin = table2array(Dataset(:, "KinematicVisvosity"));

d_equi = table2array(Dataset(:, "ParticleSize"));
size_a = table2array(Dataset(:, "a"));
size_b = table2array(Dataset(:, "b"));
size_c = table2array(Dataset(:, "c"));

shape_flt = table2array(Dataset(:, "Flatness"));
shape_eln = table2array(Dataset(:, "elongation"));
shape_del = table2array(Dataset(:, "Dellino"));
shape_sph = table2array(Dataset(:, "Sphericity"));
shape_cir = table2array(Dataset(:, "Circularity"));
Reynolds = table2array(Dataset(:, "Re"));
Powers = table2array(Dataset(:, "Powers"));

wvel_meas = table2array(Dataset(:, "Wmeasured"));

% Set up additional variables
SA_mP = zeros(140, 1);
SA_EqSph = zeros(140, 1);
Vol_mP = zeros(140, 1);
Mass_mP = zeros(140, 1);
CSF = zeros(140, 1);
rho_rel = zeros(140, 1);

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
end

g=9.81;
%% Stokes' method
% <<<<<<<<<<<<<<<<<

% Set initial velocity and timestep
 wvel_Stokes = zeros(140, 10000);
 wvel_Stokes(:, 1) = 0.0001;  % Note that earlier tests have shown that the
                              % terminal velocity is independent of the
                              % initial velocity.
timestep = 0.0002;

% Set up variable arrays
Cd_Stokes = zeros(140, 10000);
Re_Stokes = zeros(140, 10000);
Fd_Stokes = zeros(140, 10000);
Fg_Stokes = zeros(140, 10000);
Fb_Stokes = zeros(140, 10000);
Fnet_Stokes = zeros(140, 10000);
Dist1_Stokes = zeros(140, 10000);
Acc_Stokes = zeros(140, 10000);
wtFinal_Stokes = zeros(140, 1);
DistTot_Stokes = 0.0;
FinalStep_Stokes = zeros(140, 1);
FinalTime_Stokes = zeros(140, 1);

for i=1:140
    for t=1:10000
        
        Re_Stokes(i, t) = (rho_f(i)*wvel_Stokes(i, t)*d_equi(i))/vis_dyn(i);
        Cd_Stokes(t) = 24.0/Re_Stokes(i, t);
	
		Fd_Stokes(t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_Stokes(i,t)^2.0)*Cd_Stokes(t);
	
		Fg_Stokes(t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Stokes(t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Stokes(t) = Fg_Stokes(t) - Fb_Stokes(t) - Fd_Stokes(t);
	
		wvel_Stokes(i, t+1) = ((Fnet_Stokes(t)/Mass_mP(i))*timestep)+wvel_Stokes(i, t);
	
		Dist1_Stokes(i) = wvel_Stokes(i, t) * timestep;
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

% Check that the timestep is ok
% meas = zeros(FinalStep_Stokes(2), 2);
% meas(:, 2) = wvel_meas(2);
% plot(timesec(1:(FinalStep_Stokes(2))), wvel_Stokes(2, (1:(FinalStep_Stokes(2)))))
% hold on
% plot(timesec(1:(FinalStep_Stokes(1))), meas(:, 1))
%set(gca,'YLim', [0, 2.5e-3] )
%set(gca,'XLim', [0, 3.5e-3] )

% Store results required for plotting in one array
Results_Stokes = zeros(140, 5);

for i=1:140
    Results_Stokes(i, 1) = d_equi(i);
    Results_Stokes(i, 2) = CSF(i);
    Results_Stokes(i, 3) = wtFinal_Stokes(i);
    Results_Stokes(i, 4) = FinalTime_Stokes(i);
    Results_Stokes(i, 5) = wvel_meas(i);
end 

Table_Stokes = array2table(Results_Stokes, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', ...
    'Time', 'Wt_Meas'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Stokes output
% <<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_Stokes.('ESD'), Table_Stokes.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes.('ESD'), Table_Stokes.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Stokes Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')

plot(Table_Stokes.('ESD'), Table_Stokes.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes{1:80, "ESD"}, Table_Stokes{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes{81:100, "ESD"}, Table_Stokes{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes{101:140, "ESD"}, Table_Stokes{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Stokes Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% b) wt against CSF

plot(Table_Stokes.('CSF'), Table_Stokes.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Stokes.('CSF'), Table_Stokes.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'locatoin', 'best')
title('Stokes Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% c) wt against wt measured
Max(1) = max(Table_Stokes.Wt);
Max(2) = max(Table_Stokes.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_Stokes.('Wt_Meas'), Table_Stokes.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Stokes.('Wt_Meas'), Table_Stokes.('Wt'), 1);
px=[min(Table_Stokes.('Wt_Meas')) max(Table_Stokes.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Stokes.("Wt_Meas")\Table_Stokes.("Wt");
mx = m*Table_Stokes.("Wt_Meas");
plot(Table_Stokes.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Stokes Cd')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*Maximum])
set(gca, 'Xlim', [0, 1.1*Maximum])
hold off

plot(yx, yx)
hold on
plot(Table_Stokes{1:80, "Wt_Meas"}, Table_Stokes{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Stokes{81:100, "Wt_Meas"}, Table_Stokes{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes{101:140, "Wt_Meas"}, Table_Stokes{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Stokes Cd')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off

%% Dietrich Method
% <<<<<<<<<<<<<<<<<<<
% Dietrich's model cannot be used when CSF<0.2
% This is not an interative procedure, it just calculates the terminal velocity.

% Set up variable arrays
d_dim = zeros(140, 1);
wdim_Dietrich = zeros(140, 1);
wt_Dietrich = zeros(140, 1);
R1 = zeros(140,1);
R2 = zeros(140,1);
R3 = zeros(140,1);
g=9.81;
for i=1:140
        d_dim(i) = ((rho_p(i) - rho_f(i))*g*(d_equi(i)^3.0))/(rho_f(i)*(vis_kin(i)^2.0));
	
	    R1(i) = -3.76715 + 1.92944*(log10(d_dim(i))) - 0.09815*((log10(d_dim(i)))^2.0) ...
		    -0.00575*((log10(d_dim(i)))^3.0) + 0.00056*((log10(d_dim(i)))^4.0);
	    
        if(CSF(i)>0.15)
            R2(i) = (log10(1-((1-CSF(i)))/0.85)) - ((1-CSF(i))^2.3)*tanh(log10(d_dim(i))-4.6)...
		    + 0.3*(0.5-CSF(i))*((1-CSF(i))^2.0)*(log10(d_dim(i))-4.6);
        else
            R2(i)=nan;
        end
	    R3(i) = (0.65-((CSF(i)/2.83)*(tanh(log10(d_dim(i))-4.6))))^(1+((3.5-Powers(i))/2.5));
	    
	    wdim_Dietrich(i) = R3(i) * (10^(R1(i)+R2(i)));
	    
	    wt_Dietrich(i) = ((wdim_Dietrich(i)*(rho_p(i) - rho_f(i))*g*vis_kin(i))/rho_f(i))^(1.0/3.0);
end
%% Róisín, are you using projected area or surface area?? Try testing with whatever one you haven't used too.
% Store results required for plotting in one array
Results_Dietrich = zeros(140, 4);

for i=1:140
    Results_Dietrich(i, 1) = d_equi(i);
    Results_Dietrich(i, 2) = CSF(i);
    Results_Dietrich(i, 3) = wt_Dietrich(i);
    Results_Dietrich(i, 4) = wvel_meas(i);
end 

Table_Dietrich = array2table(Results_Dietrich, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Dietrich output
% <<<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Wt')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

plot(Table_Dietrich(:, "ESD"), Table_Dietrich(:, "Wt_Meas"), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich{1:80, "ESD"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "ESD"}, Table_Dietrich{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "ESD"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
    'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich Wt')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% b) wt against CSF

plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Wt')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% c) wt against wt measured
Max(1) = max(Table_Dietrich.Wt);
Max(2) = max(Table_Dietrich.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_Dietrich.('Wt_Meas'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
% p=polyfit(Table_Dietrich.('Wt_Meas'), Table_Dietrich.('Wt'), 1);
% px=[min(Table_Dietrich.('Wt_Meas')) max(Table_Dietrich.('Wt_Meas'))];
% py=polyval(p, px);
% plot(px, py, '-b')
% text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
%     'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
% m=Table_Dietrich.("Wt_Meas")\Table_Dietrich.("Wt");
% mx = m*Table_Dietrich.("Wt_Meas");
% plot(Table_Dietrich.('Wt_Meas'), mx, '-g');
% text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
%     'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
xlabel('Measured Wt (m/s)')
ylabel('Calcualted Wt (m/s)')
title('Dietrich Wt')
hold off

plot(yx, yx)
hold on
plot(Table_Dietrich{1:80, "Wt_Meas"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "Wt_Meas"}, Table_Dietrich{81:100, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "Wt_Meas"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Dietrich Wt')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off
%% Bagheri and Bonadonna (2016)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This is an iterative procedure

timestep = 0.0002;

Re_BB = zeros(140, 10000);
FormFactor_S = zeros(140, 1);
FormFactor_N = zeros(140, 1);
Correction_S = zeros(140, 1);
Correction_N = zeros(140, 1);
ratio_density = zeros(140,1);
alpha2= zeros(140,1);
beta2 = zeros(140,1);
wvel_BB = zeros(140, 10000);
Cd_BB = zeros(140, 10000);
Fd_BB = zeros(140, 10000);
Fg_BB = zeros(140, 10000);
Fb_BB = zeros(140, 10000);
Fnet_BB = zeros(140, 10000);
Dist1_BB = zeros(140, 10000);
DistTot_BB = 0.0;
Acc_BB = zeros(140, 10000);
FinalTime_BB = zeros(140, 1);
FinalStep_BB = zeros(140, 1);
wtFinal_BB = zeros(140, 1);

% Set initial velocity
wvel_BB(:, 1) = 0.0001;
for i=1:140	
    
    FormFactor_S(i) = shape_flt(i)*(shape_eln(i)^1.3)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
	FormFactor_N(i) = (shape_flt(i)^2.0)*shape_eln(i)*((d_equi(i)^3.0)/(size_a(i)*size_b(i)*size_c(i)));
		
	Correction_S(i) = 0.5*((FormFactor_S(i)^(1.0/3.0))+(FormFactor_S(i)^(-1.0/3.0)));
		
	ratio_density(i) = rho_p(i) / rho_f(i);
	alpha2(i) = 0.45 + (10.0/((exp(2.5*log10(ratio_density(i))))+30));		
    beta2(i) = 1.0 -  (37.0/((exp(3.0*log10(ratio_density(i))))+100));	
	Correction_N(i) = 10.0^(alpha2(i)*((-1.0*log10(FormFactor_N(i)))^beta2(i)));
    
    for t=1:10000
		
		Re_BB(i, t) = abs((rho_p(i) * wvel_BB(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_BB(i,t) = ((24.0*Correction_S(i))/Re_BB(i,t))*(1+ 0.125*(((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))*(2.0/3.0))) ...
			+ (0.46*Correction_N(i))/(1 + (5330/((Re_BB(i,t)*Correction_N(i))/(Correction_S(i)))));
	
		Fd_BB(i,t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_BB(i,t)^2.0)*Cd_BB(i,t);
	
		Fg_BB(i,t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_BB(i,t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_BB(i,t) = Fg_BB(i,t) - Fb_BB(i,t) - Fd_BB(i,t);
	
		wvel_BB(i,t+1) = ((Fnet_BB(i,t)/Mass_mP(i))*timestep)+wvel_BB(i,t);

        Dist1_BB(i,t) = wvel_BB(i,t) * timestep;
		DistTot_BB = DistTot_BB + Dist1_BB(i,t);
		Acc_BB(i,t) = (wvel_BB(i, t+1) - wvel_BB(i,t))/timestep;
		
        if (Acc_BB(i,t)< 0.001)
            FinalTime_BB(i) = (t+1)*timestep;
            FinalStep_BB(i) = t+1;
            wtFinal_BB(i)=wvel_BB(i, t+1);
			break
        end
    end
end

% Store results required for plotting in one array
Results_BB = zeros(140, 4);

for i=1:140
    Results_BB(i, 1) = d_equi(i);
    Results_BB(i, 2) = CSF(i);
    Results_BB(i, 3) = wtFinal_BB(i);
    Results_BB(i, 4) = wvel_meas(i);
    Results_BB(i, 5) = FinalTime_BB(i);
end 

Table_BB = array2table(Results_BB, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Bagheri output
% <<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_BB.('ESD'), Table_BB.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB.('ESD'), Table_BB.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Bagheri and Bonadonna Cd')
hold off

plot(Table_BB.("ESD"), Table_BB.("Wt_Meas"), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB{1:80, "ESD"}, Table_BB{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB{81:100, "ESD"}, Table_BB{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{101:140, "ESD"}, Table_BB{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
    'Calculated Wt, Film', 'NumColumns', 2, 'Location', 'southoutside')
title('Bagheri and Bonadonna Cd')
hold off

% b) wt against CSF

plot(Table_BB.('CSF'), Table_BB.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_BB.('CSF'), Table_BB.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Bagheri and Bonadonna Cd')
hold off

% c) wt against wt measured
Max(1) = max(Table_BB.Wt);
Max(2) = max(Table_BB.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_BB.('Wt_Meas'), Table_BB.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_BB.('Wt_Meas'), Table_BB.('Wt'), 1);
px=[min(Table_BB.('Wt_Meas')) max(Table_BB.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_BB.("Wt_Meas")\Table_BB.("Wt");
mx = m*Table_BB.("Wt_Meas");
plot(Table_BB.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
xlabel('Measured Wt')
ylabel('Calculated Wt')
title('Bagheri and Bonadonna Cd')
hold off

plot(yx, yx)
hold on
plot(Table_BB{1:80, "Wt_Meas"}, Table_BB{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB{81:100, "Wt_Meas"}, Table_BB{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{101:140, "Wt_Meas"}, Table_BB{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Bagheri and Bonadonna Cd')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'southeast')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off	
%% Dioguardi et al. (2018)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This is an iterative procedure

timestep = 0.0002;	

Re_Dio = zeros(140, 10000);
wvel_Dio = zeros(140, 10000);
Cd_Dio = zeros(140, 10000);
Fd_Dio = zeros(140, 10000);
Fg_Dio = zeros(140, 10000);
Fb_Dio = zeros(140, 10000);
Fnet_Dio = zeros(140, 10000);
Dist1_Dio = zeros(140, 10000);
DistTot_Dio=0.0;
Acc_Dio = zeros(140, 10000);	
FinalTime_Dio = zeros(140, 1);
FinalStep_Dio = zeros(140, 1);
wtFinal_Dio = zeros(140, 1);

% Set initial values
wvel_Dio(:, 1) = 0.0001;

for i=1:140
    for t=1:10000
        
        Re_Dio(i,t) = abs((rho_p(i) * wvel_Dio(i, t) * d_equi(i))/ vis_dyn(i));
		
		Cd_Dio(i, t) = (24.0/Re_Dio(i,t))*(((1.0-shape_del(i))/(Re_Dio(i,t)+1.0))^0.25) ...
			     + (24.0/Re_Dio(i,t))*0.1806*(Re_Dio(i,t)^0.6459)*(shape_del(i)^(-1.0*(Re_Dio(i,t)^0.08))) ...
			     + 0.4251/(1.0+((6880.95/Re_Dio(i,t))*(shape_del(i)^5.05)));
	
		Fd_Dio(i,t) = 0.5*rho_f(i)*SA_mP(i)*(wvel_Dio(i,t)^2.0)*Cd_Dio(i,t);
	
		Fg_Dio(i,t) = Vol_mP(i)*rho_p(i)*g;
	
		Fb_Dio(i,t) = Vol_mP(i)*rho_f(i)*g;
	
		Fnet_Dio(i,t) = Fg_Dio(i,t) - Fb_Dio(i,t) - Fd_Dio(i,t);
	
		wvel_Dio(i,t+1) = ((Fnet_Dio(i,t)/Mass_mP(i))*timestep)+wvel_Dio(i,t);
	
		Dist1_Dio(i,t+1) = wvel_Dio(i,t+1) * timestep;
		DistTot_Dio = DistTot_Dio + Dist1_Dio(t);
		Acc_Dio(i,t) = (wvel_Dio(i,t+1) - wvel_Dio(i,t))/timestep;
		
        if (Acc_Dio(i,t)< 0.001)
            FinalTime_Dio(i) = (t+1)*timestep;
            FinalStep_Dio(i) = t+1;
            wtFinal_Dio(i)=wvel_Dio(i, t+1);
            break
        end
    end
end

% Store results required for plotting in one array
Results_Dio = zeros(140, 4);

for i=1:140
    Results_Dio(i, 1) = d_equi(i);
    Results_Dio(i, 2) = CSF(i);
    Results_Dio(i, 3) = wtFinal_Dio(i);
    Results_Dio(i, 4) = wvel_meas(i);
    Results_Dio(i, 5) = FinalTime_Dio(i);
end 

Table_Dio = array2table(Results_Dio, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Dioguardi output
% <<<<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_Dio.('ESD'), Table_Dio.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dio.('ESD'), Table_Dio.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Dioguardi Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

plot(Table_Dio.("ESD"), Table_Dio.("Wt_Meas"), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dio{1:80, "ESD"}, Table_Dio{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio{81:100, "ESD"}, Table_Dio{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{101:140, "ESD"}, Table_Dio{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
    'Calculated Wt, Film', 'Location', 'best')
title('Dioguardi Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% b) wt against CSF

plot(Table_Dio.('CSF'), Table_Dio.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dio.('CSF'), Table_Dio.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Dioguardi Cd')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% c) wt against wt measured
Max(1) = max(Table_Dio.Wt);
Max(2) = max(Table_Dio.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_Dio.('Wt_Meas'), Table_Dio.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Dio.('Wt_Meas'), Table_Dio.('Wt'), 1);
px=[min(Table_Dio.('Wt_Meas')) max(Table_Dio.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Dio.("Wt_Meas")\Table_Dio.("Wt");
mx = m*Table_Dio.("Wt_Meas");
plot(Table_Dio.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
legend('Measured Wt', '', 'Location', 'best')
xlabel('Measured Settling Velocity (m/s)')
ylabel('Calculated Settling Velocity (m/s)')
title('Dioguardi Cd')
set(gca, 'Xlim', [0, inf])
set(gca, 'Ylim', [0, inf])
hold off

plot(yx, yx)
hold on
plot(Table_Dio{1:80, "Wt_Meas"}, Table_Dio{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio{81:100, "Wt_Meas"}, Table_Dio{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{101:140, "Wt_Meas"}, Table_Dio{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Dioguardi Cd')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'southeast')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off

%% Francalanci et al (2021)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This is not an iterative method

Dref_Frn = zeros(140, 1);
Ddim_Frn = zeros(140, 1);
wdim_Frn = zeros(140, 1);
wvel_Frn = zeros(140, 1);

for i=1:140
    Dref_Frn(i) = size_a(i)*((CSF(i))^0.34)*((size_b(i)/size_a(i))^0.5);
    Ddim_Frn(i) = Dref_Frn(i)*((g*rho_rel(i))/(vis_kin(i)^2.0))^(1.0/3.0);
	
	E = size_a(i)*((((size_a(i)^2.0)+(size_b(i)^2.0)+(size_c(i)^2.0))/3.0)^-0.5);
	C1 = 18.0*(E^-0.38);
	C2 = 0.3708*(CSF(i)^-0.1602);
	n = 0.4942*(CSF(i)^-0.059);
	
	wdim_Frn(i) = (Ddim_Frn(i)^2)/(C1+(0.75*C2*(Ddim_Frn(i)^3))^n);
	wvel_Frn(i) = ((wdim_Frn(i)*(rho_p(i) - rho_f(i))*g*vis_kin(i))/rho_f(i))^(1.0/3.0);
end

% Store results required for plotting in one array
Results_Frn = zeros(140, 4);

for i=1:140
    Results_Frn(i, 1) = d_equi(i);
    Results_Frn(i, 2) = CSF(i);
    Results_Frn(i, 3) = wvel_Frn(i);
    Results_Frn(i, 4) = wvel_meas(i);
end 

Table_Frn = array2table(Results_Frn, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Francalanci output
% <<<<<<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_Frn.('ESD'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn.('ESD'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

plot(Table_Frn.("ESD"), Table_Frn.("Wt_Meas"), 'ok')
hold on
plot(Table_Frn{1:80, "ESD"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "ESD"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "ESD"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Fragment', 'Fibre', 'Film', 'Location', 'best')
title('Francalanci Wt')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% b) wt against CSF

plot(Table_Frn.('CSF'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn.('CSF'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

% c) wt against wt measured
Max(1) = max(Table_Frn.Wt);
Max(2) = max(Table_Frn.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_Frn.('Wt_Meas'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Frn.('Wt_Meas'), Table_Frn.('Wt'), 1);
px=[min(Table_Frn.('Wt_Meas')) max(Table_Frn.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Frn.("Wt_Meas")\Table_Frn.("Wt");
mx = m*Table_Frn.("Wt_Meas");
plot(Table_Frn.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
title('Francalanci Wt')
hold off

plot(yx, yx)
hold on
plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Francalanci Wt')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off	

%% Zhang and Choi (2021)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This is an iterative procedure

timestep = 0.0002;

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
DistTot_ZC = 0.0;
Acc_ZC = zeros(140, 10000);
FinalTime_ZC = zeros(140, 1);
FinalStep_ZC = zeros(140, 1);
wtFinal_ZC = zeros(140, 1);

% Set initial velocity
wvel_ZC(:, 1) = 0.0001;
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
		DistTot_ZC = DistTot_ZC + Dist1_ZC(i,t);
		Acc_ZC(i,t) = (wvel_ZC(i, t+1) - wvel_ZC(i,t))/timestep;
		
        if (Acc_ZC(i,t)< 0.001)
            FinalTime_ZC(i) = (t+1)*timestep;
            FinalStep_ZC(i) = t+1;
            wtFinal_ZC(i)=wvel_ZC(i, t+1);
			break
        end
    end
end

% Store results required for plotting in one array
Results_ZC = zeros(140, 4);

for i=1:140
    Results_ZC(i, 1) = d_equi(i);
    Results_ZC(i, 2) = CSF(i);
    Results_ZC(i, 3) = wtFinal_ZC(i);
    Results_ZC(i, 4) = wvel_meas(i);
    Results_ZC(i, 5) = FinalTime_ZC(i);
end 

Table_ZC = array2table(Results_ZC, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas', 'Time'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140

%% Plot Zhang output
% <<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_ZC.('ESD'), Table_ZC.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_ZC.('ESD'), Table_ZC.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Zhang and Choi Cd')
hold off

plot(Table_ZC.("ESD"), Table_ZC.("Wt_Meas"), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_ZC{1:80, "ESD"}, Table_ZC{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC{81:100, "ESD"}, Table_ZC{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC{101:140, "ESD"}, Table_ZC{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
    'Calculated Wt, Film', 'NumColumns', 2, 'Location', 'southoutside')
title('Zhang and Choi Cd')
hold off

% b) wt against CSF

plot(Table_ZC.('CSF'), Table_ZC.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_ZC.('CSF'), Table_ZC.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
title('Zhang and Choi Cd')
hold off

% c) wt against wt measured
Max(1) = max(Table_ZC.Wt);
Max(2) = max(Table_ZC.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_ZC.('Wt_Meas'), Table_ZC.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_ZC.('Wt_Meas'), Table_ZC.('Wt'), 1);
px=[min(Table_ZC.('Wt_Meas')) max(Table_ZC.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_ZC.("Wt_Meas")\Table_ZC.("Wt");
mx = m*Table_ZC.("Wt_Meas");
plot(Table_ZC.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
xlabel('Measured Wt')
ylabel('Calculated Wt')
title('Zhang and Choi Cd')
hold off

plot(yx, yx)
hold on
plot(Table_ZC{1:80, "Wt_Meas"}, Table_ZC{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC{81:100, "Wt_Meas"}, Table_ZC{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC{101:140, "Wt_Meas"}, Table_ZC{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Zhang and Choi Cd')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'southeast')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off	

%% Yu et al. (2022)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This is not an iterative procedure

d_dimYu = zeros(140, 1);
wvel_Yu = zeros(140, 1);
CdSph_Yu = zeros(140, 1);
Cd_Yu = zeros(140, 1);

for i=1:140	
    d_dimYu(i) = (((rho_rel(i)*g)/(vis_kin(i)^2.0))^(1.0/3.0))*d_equi(i);
    CdSph_Yu(i) = (432.0/(d_dimYu(i)^3.0))*((1 + 0.022*(d_dimYu(i)^3.0))^0.54)...
                   + (0.47*(1- exp(-0.15*(d_dimYu(i)^0.45))));
    Cd_Yu(i) = CdSph_Yu(i)/(((d_dimYu(i)^-0.25)*(shape_sph(i)^(d_dimYu(i)^0.03))*(CSF(i)^(d_dimYu(i)^0.33)))^0.25);
    wvel_Yu(i) = ((vis_kin(i)*g*rho_rel(i))^(1.0/3.0))*(((4*d_dimYu(i))/(3*Cd_Yu(i)))^(0.5));
end

% Store results required for plotting in one array
Results_Yu = zeros(140, 4);

for i=1:140
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = wvel_meas(i);
end 

Table_Yu = array2table(Results_Yu, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

% Note that the shapes are in the following rows of the table:
% Fragments: 1:80
% Fibres: 81:100
% Film: 101:140


%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<<<<<<

% a) wt against ESD

plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('ESD'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'Location', 'best')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
title('Yu Wt')
hold off

plot(Table_Yu.("ESD"), Table_Yu.("Wt_Meas"), 'ok')
hold on
plot(Table_Yu{1:80, "ESD"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "ESD"}, Table_Yu{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "ESD"}, Table_Yu{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured', 'Fragment', 'Fibre', 'Film', 'Location', 'best')
title('Yu Wt')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

% b) wt against CSF

plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('CSF'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
title('Yu Wt')
hold off

% c) wt against wt measured
Max(1) = max(Table_Yu.Wt);
Max(2) = max(Table_Yu.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(Table_Yu.('Wt_Meas'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Yu.('Wt_Meas'), Table_Yu.('Wt'), 1);
px=[min(Table_Yu.('Wt_Meas')) max(Table_Yu.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Yu.("Wt_Meas")\Table_Yu.("Wt");
mx = m*Table_Yu.("Wt_Meas");
plot(Table_Yu.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
title('Yu Wt')
hold off

plot(yx, yx)
hold on
plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Yu Wt')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off

%% Plot all models Wt calc Vs Wt meas

% Dietrich
subplot(2, 3, 1)
Max(1) = max(Table_Dietrich.Wt);
Max(2) = max(Table_Dietrich.Wt_Meas);
Maximum = max(Max);

yx=linspace(0, Maximum, 100);
plot(yx, yx, '--k')
hold on
plot(Table_Dietrich{1:80, "Wt_Meas"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "Wt_Meas"}, Table_Dietrich{81:100, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "Wt_Meas"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_Dietrich.('Wt_Meas'), Table_Dietrich.('Wt'), 1);
px=[min(Table_Dietrich.('Wt_Meas')) max(Table_Dietrich.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Dietrich.("Wt_Meas")\Table_Dietrich.("Wt");
mx = m*Table_Dietrich.("Wt_Meas");
plot(Table_Dietrich.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
%legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
title('Dietrich Wt')
hold off

% Bagheri
subplot(2, 3, 2)
Max(1) = max(Table_BB.Wt);
Max(2) = max(Table_BB.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(yx, yx, '--k')
hold on
plot(Table_BB{1:80, "Wt_Meas"}, Table_BB{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_BB{81:100, "Wt_Meas"}, Table_BB{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB{101:140, "Wt_Meas"}, Table_BB{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_BB.('Wt_Meas'), Table_BB.('Wt'), 1);
px=[min(Table_BB.('Wt_Meas')) max(Table_BB.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_BB.("Wt_Meas")\Table_BB.("Wt");
mx = m*Table_BB.("Wt_Meas");
plot(Table_BB.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Bagheri and Bonadonna Cd')
%legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off	

% Dioguardi
subplot(2, 3, 3)

Max(1) = max(Table_Dio.Wt);
Max(2) = max(Table_Dio.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(yx, yx, '--k')
hold on
plot(Table_Dio{1:80, "Wt_Meas"}, Table_Dio{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dio{81:100, "Wt_Meas"}, Table_Dio{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio{101:140, "Wt_Meas"}, Table_Dio{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_Dio.('Wt_Meas'), Table_Dio.('Wt'), 1);
px=[min(Table_Dio.('Wt_Meas')) max(Table_Dio.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Dio.("Wt_Meas")\Table_Dio.("Wt");
mx = m*Table_Dio.("Wt_Meas");
plot(Table_Dio.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Dioguardi Cd')
%legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off

% Francalanci
subplot(2, 3, 4)
Max(1) = max(Table_Frn.Wt);
Max(2) = max(Table_Frn.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(yx, yx, '--k')
hold on
plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_Frn.('Wt_Meas'), Table_Frn.('Wt'), 1);
px=[min(Table_Frn.('Wt_Meas')) max(Table_Frn.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Frn.("Wt_Meas")\Table_Frn.("Wt");
mx = m*Table_Frn.("Wt_Meas");
plot(Table_Frn.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Francalanci Wt')
%legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off	

% Zhang
subplot(2, 3, 5)

Max(1) = max(Table_ZC.Wt);
Max(2) = max(Table_ZC.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(yx, yx, '--k')
hold on
plot(Table_ZC{1:80, "Wt_Meas"}, Table_ZC{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_ZC{81:100, "Wt_Meas"}, Table_ZC{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_ZC{101:140, "Wt_Meas"}, Table_ZC{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_ZC.('Wt_Meas'), Table_ZC.('Wt'), 1);
px=[min(Table_ZC.('Wt_Meas')) max(Table_ZC.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_ZC.("Wt_Meas")\Table_ZC.("Wt");
mx = m*Table_ZC.("Wt_Meas");
plot(Table_ZC.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Zhang and Choi Cd')
%legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'location', 'best')
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off

% Yu
subplot(2, 3, 6)

Max(1) = max(Table_Yu.Wt);
Max(2) = max(Table_Yu.Wt_Meas);
Maximum = max(Max);
yx=linspace(0, Maximum, 100);

plot(yx, yx, '--k')
hold on
plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
p=polyfit(Table_Yu.('Wt_Meas'), Table_Yu.('Wt'), 1);
px=[min(Table_Yu.('Wt_Meas')) max(Table_Yu.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Yu.("Wt_Meas")\Table_Yu.("Wt");
mx = m*Table_Yu.("Wt_Meas");
plot(Table_Yu.('Wt_Meas'), mx, '-g');
text(px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Yu Wt')
legend('y=x', 'Fragment', 'Fibre', 'Film', '', '', 'Position', [0.915798607147828,0.1091079007951,0.064973959090809,0.088541668602544])
set(gca,'YLim', [0, Maximum*1.1] )
set(gca,'XLim', [0, Maximum*1.1] )
hold off


han=axes(gcf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Calculated settling velocity (m/s)');
xlabel(han, 'Measured settling velocity (m/s)')
set(gca, 'XAxisLocation', 'bottom')     
set(gcf, 'WindowState', 'maximized');
%title(han,'Wt meas Vs wt calc, various models');
axesHandles = findobj(get(gcf,'Children'), 'flat');


exportgraphics(gcf, 'C:\Users\roisi\Desktop\mP Model Code\DragModelsTest\Modelsx6.jpeg', 'Resolution', 300)
