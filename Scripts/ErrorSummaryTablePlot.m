%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: All_Error
% Date created: 30.06.22
% Date last mostified: 30.06.22
% Purpose: To compile all the error tables into one
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in tables
clear
BB_SA = readtable('./Output/20220621/Bagheri/BagheriFinalTableVM_SA.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);
BB_Proj = readtable('./Output/20220621/Bagheri/BagheriFinalTableVM_Proj.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Dietrich = readtable('./Output/20220621/Dietrich/DietrichFinalTableVM.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Dio_SA = readtable('./Output/20220621/Dioguardi/DioFinalTableVM_SA.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);
Dio_Proj = readtable('./Output/20220621/Dioguardi/DioFinalTableVM_Proj.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Frn = readtable('./Output/20220621/Francalanci/FrancalanciFinalTableVM.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Stokes_SA = readtable('./Output/20220621/Stokes/StokesFinalTableVM_SA.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);
Stokes_Proj = readtable('./Output/20220621/Stokes/StokesFinalTableVM_Proj.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Yu = readtable('./Output/20220621/Yu/YuFinalTableVM.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

Zhang_SA = readtable('./Output/20220621/Zhang/ZhangFinalTableVM_SA.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);
Zhang_Proj = readtable('./Output/20220621/Zhang/ZhangFinalTableVM_Proj.txt', Delimiter=',', ReadVariableNames=true, ReadRowNames=true);

%% Make an empty table

Col_names = ["Stokes_Proj", "Stokes_SA",  "Dio_Proj", "Dio_SA", ...
             "BB_Proj", "BB_SA", "Zhang_Proj", "Zhang_SA",  "Dietrich", "Francalanci", "Yu"];
Row_names = ["All_AE", "All_AbsAE", "All_RMSE", "All_m", "All_r_sq", ...
             "Fragment_AE", "Fragment_AbsAE", "Fragment_RMSE", "Fragment_m", "Fragment_r_sq", ...
             "Fibre_AE", "Fibre_AbsAE", "Fibre_RMSE", "Fibre_m", "Fibre_r_sq", ...
             "Film_AE", "Film_AbsAE", "Film_RMSE", "Film_m", "Film_r_sq"];
Var_types = ["double","double","double","double","double","double", ...
             "double","double","double","double","double"];

Summary_table = table('Size', [20 11], 'VariableTypes', Var_types);
Summary_table.Properties.VariableNames = Col_names;
Summary_table.Properties.RowNames = Row_names;

%% Populate table

Summary_table("All_AE", "Stokes_Proj") = Stokes_Proj("All", "Error_table_AE");
Summary_table("All_AbsAE", "Stokes_Proj") = Stokes_Proj("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Stokes_Proj") = Stokes_Proj("All", "Error_table_RMSE");
Summary_table("All_m", "Stokes_Proj") = Stokes_Proj("All", "m");
Summary_table("All_r_sq", "Stokes_Proj") = Stokes_Proj("All", "r_sq");
Summary_table("Fragment_AE", "Stokes_Proj") = Stokes_Proj("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Stokes_Proj") = Stokes_Proj("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Stokes_Proj") = Stokes_Proj("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Stokes_Proj") = Stokes_Proj("Fragment", "m");
Summary_table("Fragment_r_sq", "Stokes_Proj") = Stokes_Proj("Fragment", "r_sq");
Summary_table("Fibre_AE", "Stokes_Proj") = Stokes_Proj("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Stokes_Proj") = Stokes_Proj("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Stokes_Proj") = Stokes_Proj("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Stokes_Proj") = Stokes_Proj("Fibre", "m");
Summary_table("Fibre_r_sq", "Stokes_Proj") = Stokes_Proj("Fibre", "r_sq");
Summary_table("Film_AE", "Stokes_Proj") = Stokes_Proj("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Stokes_Proj") = Stokes_Proj("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Stokes_Proj") = Stokes_Proj("Film", "Error_table_RMSE");
Summary_table("Film_m", "Stokes_Proj") = Stokes_Proj("Film", "m");
Summary_table("Film_r_sq", "Stokes_Proj") = Stokes_Proj("Film", "r_sq");

Summary_table("All_AE", "Stokes_SA") = Stokes_SA("All", "Error_table_AE");
Summary_table("All_AbsAE", "Stokes_SA") = Stokes_SA("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Stokes_SA") = Stokes_SA("All", "Error_table_RMSE");
Summary_table("All_m", "Stokes_SA") = Stokes_SA("All", "m");
Summary_table("All_r_sq", "Stokes_SA") = Stokes_SA("All", "r_sq");
Summary_table("Fragment_AE", "Stokes_SA") = Stokes_SA("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Stokes_SA") = Stokes_SA("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Stokes_SA") = Stokes_SA("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Stokes_SA") = Stokes_SA("Fragment", "m");
Summary_table("Fragment_r_sq", "Stokes_SA") = Stokes_SA("Fragment", "r_sq");
Summary_table("Fibre_AE", "Stokes_SA") = Stokes_SA("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Stokes_SA") = Stokes_SA("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Stokes_SA") = Stokes_SA("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Stokes_SA") = Stokes_SA("Fibre", "m");
Summary_table("Fibre_r_sq", "Stokes_SA") = Stokes_SA("Fibre", "r_sq");
Summary_table("Film_AE", "Stokes_SA") = Stokes_SA("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Stokes_SA") = Stokes_SA("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Stokes_SA") = Stokes_SA("Film", "Error_table_RMSE");
Summary_table("Film_m", "Stokes_SA") = Stokes_SA("Film", "m");
Summary_table("Film_r_sq", "Stokes_SA") = Stokes_SA("Film", "r_sq");

Summary_table("All_AE", "Dio_Proj") = Dio_Proj("All", "Error_table_AE");
Summary_table("All_AbsAE", "Dio_Proj") = Dio_Proj("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Dio_Proj") = Dio_Proj("All", "Error_table_RMSE");
Summary_table("All_m", "Dio_Proj") = Dio_Proj("All", "m");
Summary_table("All_r_sq", "Dio_Proj") = Dio_Proj("All", "r_sq");
Summary_table("Fragment_AE", "Dio_Proj") = Dio_Proj("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Dio_Proj") = Dio_Proj("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Dio_Proj") = Dio_Proj("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Dio_Proj") = Dio_Proj("Fragment", "m");
Summary_table("Fragment_r_sq", "Dio_Proj") = Dio_Proj("Fragment", "r_sq");
Summary_table("Fibre_AE", "Dio_Proj") = Dio_Proj("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Dio_Proj") = Dio_Proj("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Dio_Proj") = Dio_Proj("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Dio_Proj") = Dio_Proj("Fibre", "m");
Summary_table("Fibre_r_sq", "Dio_Proj") = Dio_Proj("Fibre", "r_sq");
Summary_table("Film_AE", "Dio_Proj") = Dio_Proj("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Dio_Proj") = Dio_Proj("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Dio_Proj") = Dio_Proj("Film", "Error_table_RMSE");
Summary_table("Film_m", "Dio_Proj") = Dio_Proj("Film", "m");
Summary_table("Film_r_sq", "Dio_Proj") = Dio_Proj("Film", "r_sq");

Summary_table("All_AE", "Dio_SA") = Dio_SA("All", "Error_table_AE");
Summary_table("All_AbsAE", "Dio_SA") = Dio_SA("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Dio_SA") = Dio_SA("All", "Error_table_RMSE");
Summary_table("All_m", "Dio_SA") = Dio_SA("All", "m");
Summary_table("All_r_sq", "Dio_SA") = Dio_SA("All", "r_sq");
Summary_table("Fragment_AE", "Dio_SA") = Dio_SA("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Dio_SA") = Dio_SA("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Dio_SA") = Dio_SA("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Dio_SA") = Dio_SA("Fragment", "m");
Summary_table("Fragment_r_sq", "Dio_SA") = Dio_SA("Fragment", "r_sq");
Summary_table("Fibre_AE", "Dio_SA") = Dio_SA("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Dio_SA") = Dio_SA("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Dio_SA") = Dio_SA("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Dio_SA") = Dio_SA("Fibre", "m");
Summary_table("Fibre_r_sq", "Dio_SA") = Dio_SA("Fibre", "r_sq");
Summary_table("Film_AE", "Dio_SA") = Dio_SA("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Dio_SA") = Dio_SA("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Dio_SA") = Dio_SA("Film", "Error_table_RMSE");
Summary_table("Film_m", "Dio_SA") = Dio_SA("Film", "m");
Summary_table("Film_r_sq", "Dio_SA") = Dio_SA("Film", "r_sq");

Summary_table("All_AE", "BB_Proj") = BB_Proj("All", "Error_table_AE");
Summary_table("All_AbsAE", "BB_Proj") = BB_Proj("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "BB_Proj") = BB_Proj("All", "Error_table_RMSE");
Summary_table("All_m", "BB_Proj") = BB_Proj("All", "m");
Summary_table("All_r_sq", "BB_Proj") = BB_Proj("All", "r_sq");
Summary_table("Fragment_AE", "BB_Proj") = BB_Proj("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "BB_Proj") = BB_Proj("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "BB_Proj") = BB_Proj("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "BB_Proj") = BB_Proj("Fragment", "m");
Summary_table("Fragment_r_sq", "BB_Proj") = BB_Proj("Fragment", "r_sq");
Summary_table("Fibre_AE", "BB_Proj") = BB_Proj("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "BB_Proj") = BB_Proj("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "BB_Proj") = BB_Proj("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "BB_Proj") = BB_Proj("Fibre", "m");
Summary_table("Fibre_r_sq", "BB_Proj") = BB_Proj("Fibre", "r_sq");
Summary_table("Film_AE", "BB_Proj") = BB_Proj("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "BB_Proj") = BB_Proj("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "BB_Proj") = BB_Proj("Film", "Error_table_RMSE");
Summary_table("Film_m", "BB_Proj") = BB_Proj("Film", "m");
Summary_table("Film_r_sq", "BB_Proj") = BB_Proj("Film", "r_sq");

Summary_table("All_AE", "BB_SA") = BB_SA("All", "Error_table_AE");
Summary_table("All_AbsAE", "BB_SA") = BB_SA("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "BB_SA") = BB_SA("All", "Error_table_RMSE");
Summary_table("All_m", "BB_SA") = BB_SA("All", "m");
Summary_table("All_r_sq", "BB_SA") = BB_SA("All", "r_sq");
Summary_table("Fragment_AE", "BB_SA") = BB_SA("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "BB_SA") = BB_SA("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "BB_SA") = BB_SA("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "BB_SA") = BB_SA("Fragment", "m");
Summary_table("Fragment_r_sq", "BB_SA") = BB_SA("Fragment", "r_sq");
Summary_table("Fibre_AE", "BB_SA") = BB_SA("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "BB_SA") = BB_SA("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "BB_SA") = BB_SA("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "BB_SA") = BB_SA("Fibre", "m");
Summary_table("Fibre_r_sq", "BB_SA") = BB_SA("Fibre", "r_sq");
Summary_table("Film_AE", "BB_SA") = BB_SA("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "BB_SA") = BB_SA("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "BB_SA") = BB_SA("Film", "Error_table_RMSE");
Summary_table("Film_m", "BB_SA") = BB_SA("Film", "m");
Summary_table("Film_r_sq", "BB_SA") = BB_SA("Film", "r_sq");

Summary_table("All_AE", "Zhang_Proj") = Zhang_Proj("All", "Error_table_AE");
Summary_table("All_AbsAE", "Zhang_Proj") = Zhang_Proj("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Zhang_Proj") = Zhang_Proj("All", "Error_table_RMSE");
Summary_table("All_m", "Zhang_Proj") = Zhang_Proj("All", "m");
Summary_table("All_r_sq", "Zhang_Proj") = Zhang_Proj("All", "r_sq");
Summary_table("Fragment_AE", "Zhang_Proj") = Zhang_Proj("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Zhang_Proj") = Zhang_Proj("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Zhang_Proj") = Zhang_Proj("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Zhang_Proj") = Zhang_Proj("Fragment", "m");
Summary_table("Fragment_r_sq", "Zhang_Proj") = Zhang_Proj("Fragment", "r_sq");
Summary_table("Fibre_AE", "Zhang_Proj") = Zhang_Proj("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Zhang_Proj") = Zhang_Proj("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Zhang_Proj") = Zhang_Proj("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Zhang_Proj") = Zhang_Proj("Fibre", "m");
Summary_table("Fibre_r_sq", "Zhang_Proj") = Zhang_Proj("Fibre", "r_sq");
Summary_table("Film_AE", "Zhang_Proj") = Zhang_Proj("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Zhang_Proj") = Zhang_Proj("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Zhang_Proj") = Zhang_Proj("Film", "Error_table_RMSE");
Summary_table("Film_m", "Zhang_Proj") = Zhang_Proj("Film", "m");
Summary_table("Film_r_sq", "Zhang_Proj") = Zhang_Proj("Film", "r_sq");

Summary_table("All_AE", "Zhang_SA") = Zhang_SA("All", "Error_table_AE");
Summary_table("All_AbsAE", "Zhang_SA") = Zhang_SA("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Zhang_SA") = Zhang_SA("All", "Error_table_RMSE");
Summary_table("All_m", "Zhang_SA") = Zhang_SA("All", "m");
Summary_table("All_r_sq", "Zhang_SA") = Zhang_SA("All", "r_sq");
Summary_table("Fragment_AE", "Zhang_SA") = Zhang_SA("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Zhang_SA") = Zhang_SA("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Zhang_SA") = Zhang_SA("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Zhang_SA") = Zhang_SA("Fragment", "m");
Summary_table("Fragment_r_sq", "Zhang_SA") = Zhang_SA("Fragment", "r_sq");
Summary_table("Fibre_AE", "Zhang_SA") = Zhang_SA("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Zhang_SA") = Zhang_SA("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Zhang_SA") = Zhang_SA("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Zhang_SA") = Zhang_SA("Fibre", "m");
Summary_table("Fibre_r_sq", "Zhang_SA") = Zhang_SA("Fibre", "r_sq");
Summary_table("Film_AE", "Zhang_SA") = Zhang_SA("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Zhang_SA") = Zhang_SA("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Zhang_SA") = Zhang_SA("Film", "Error_table_RMSE");
Summary_table("Film_m", "Zhang_SA") = Zhang_SA("Film", "m");
Summary_table("Film_r_sq", "Zhang_SA") = Zhang_SA("Film", "r_sq");

nan_m = [nan nan nan nan nan nan];
nan_t = array2table(nan_m);

Summary_table("All_AE", "Dietrich") = Dietrich("All", "Error_table_AE");
Summary_table("All_AbsAE", "Dietrich") = Dietrich("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Dietrich") = Dietrich("All", "Error_table_RMSE");
Summary_table("All_m", "Dietrich") = Dietrich("All", "m");
Summary_table("All_r_sq", "Dietrich") = Dietrich("All", "r_sq");
Summary_table("Fragment_AE", "Dietrich") = nan_t(:,1);
Summary_table("Fragment_AbsAE", "Dietrich") = nan_t(:,1);
Summary_table("Fragment_RMSE", "Dietrich") = nan_t(:,2);
Summary_table("Fragment_m", "Dietrich") = nan_t(:,1);
Summary_table("Fragment_r_sq", "Dietrich") = nan_t(:,2);
Summary_table("Fibre_AE", "Dietrich") = nan_t(:,3);
Summary_table("Fibre_AbsAE", "Dietrich") = nan_t(:,3);
Summary_table("Fibre_RMSE", "Dietrich") = nan_t(:,4);
Summary_table("Fibre_m", "Dietrich") = nan_t(:,3);
Summary_table("Fibre_r_sq", "Dietrich") = nan_t(:,4);
Summary_table("Film_AE", "Dietrich") = nan_t(:,5);
Summary_table("Film_AbsAE", "Dietrich") = nan_t(:,5);
Summary_table("Film_RMSE", "Dietrich") = nan_t(:,6);
Summary_table("Film_m", "Dietrich") = nan_t(:,5);
Summary_table("Film_r_sq", "Dietrich") = nan_t(:,6);

Summary_table("All_AE", "Francalanci") = Frn("All", "Error_table_AE");
Summary_table("All_AbsAE", "Francalanci") = Frn("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Francalanci") = Frn("All", "Error_table_RMSE");
Summary_table("All_m", "Francalanci") = Frn("All", "m");
Summary_table("All_r_sq", "Francalanci") = Frn("All", "r_sq");
Summary_table("Fragment_AE", "Francalanci") = Frn("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Francalanci") = Frn("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Francalanci") = Frn("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Francalanci") = Frn("Fragment", "m");
Summary_table("Fragment_r_sq", "Francalanci") = Frn("Fragment", "r_sq");
Summary_table("Fibre_AE", "Francalanci") = Frn("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Francalanci") = Frn("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Francalanci") = Frn("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Francalanci") = Frn("Fibre", "m");
Summary_table("Fibre_r_sq", "Francalanci") = Frn("Fibre", "r_sq");
Summary_table("Film_AE", "Francalanci") = Frn("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Francalanci") = Frn("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Francalanci") = Frn("Film", "Error_table_RMSE");
Summary_table("Film_m", "Francalanci") = Frn("Film", "m");
Summary_table("Film_r_sq", "Francalanci") = Frn("Film", "r_sq");

Summary_table("All_AE", "Yu") = Yu("All", "Error_table_AE");
Summary_table("All_AbsAE", "Yu") = Yu("All", "Error_table_Abs_AE");
Summary_table("All_RMSE", "Yu") = Yu("All", "Error_table_RMSE");
Summary_table("All_m", "Yu") = Yu("All", "m");
Summary_table("All_r_sq", "Yu") = Yu("All", "r_sq");
Summary_table("Fragment_AE", "Yu") = Yu("Fragment", "Error_table_AE");
Summary_table("Fragment_AbsAE", "Yu") = Yu("Fragment", "Error_table_Abs_AE");
Summary_table("Fragment_RMSE", "Yu") = Yu("Fragment", "Error_table_RMSE");
Summary_table("Fragment_m", "Yu") = Yu("Fragment", "m");
Summary_table("Fragment_r_sq", "Yu") = Yu("Fragment", "r_sq");
Summary_table("Fibre_AE", "Yu") = Yu("Fibre", "Error_table_AE");
Summary_table("Fibre_AbsAE", "Yu") = Yu("Fibre", "Error_table_Abs_AE");
Summary_table("Fibre_RMSE", "Yu") = Yu("Fibre", "Error_table_RMSE");
Summary_table("Fibre_m", "Yu") = Yu("Fibre", "m");
Summary_table("Fibre_r_sq", "Yu") = Yu("Fibre", "r_sq");
Summary_table("Film_AE", "Yu") = Yu("Film", "Error_table_AE");
Summary_table("Film_AbsAE", "Yu") = Yu("Film", "Error_table_Abs_AE");
Summary_table("Film_RMSE", "Yu") = Yu("Film", "Error_table_RMSE");
Summary_table("Film_m", "Yu") = Yu("Film", "m");
Summary_table("Film_r_sq", "Yu") = Yu("Film", "r_sq");

%% Write table

writetable(Summary_table, './Output/20220621/OverallSummaryTable.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Summary_table, './Output/20220621/OverallSummaryTable.xls', 'WriteRowNames', true);

%% Read table in to plot
clear
results_t = readtable('./Output/20220621/OverallSummaryTable.txt', ...
    'Delimiter', ',', 'ReadRowNames', true, 'ReadVariableNames', true);
results_m = table2array(results_t);
%% Make labels for graph

tick_label = ["Stokes:Proj", "Stokes:SA",  "Dio:Proj", "Dio:SA", "BB:Proj", ...
              "BB:SA", "Zhang:Proj", "Zhang:SA",  "Dietrich", "Francalanci", "Yu"];

labelsAE_All = "";
a = num2str(results_m(1, 1), '%7.2f');
labelsAE_All = convertCharsToStrings(a);
labelsAE_All(1);
for i = 2:11
    a = num2str(results_m(1, i), '%7.2f');
    labelsAE_All(i) = convertCharsToStrings(a);
end

labelsAbs_AE_All = "";
a = num2str(results_m(2, 1), '%7.2f');
labelsAbs_AE_All = convertCharsToStrings(a);
labelsAbs_AE_All(1);
for i = 2:11
    a = num2str(results_m(2, i), '%7.2f');
    labelsAbs_AE_All(i) = convertCharsToStrings(a);
end

labelsRMSE_All = "";
a = num2str(results_m(3, 1), '%7.2f');
labelsRMSE_All = convertCharsToStrings(a);
labelsRMSE_All(1);
for i = 2:11
    a = num2str(results_m(3, i), '%7.2f');
    labelsRMSE_All(i) = convertCharsToStrings(a);
end

labelsAE_F3 = "";
a = num2str(results_m(6, 1), '%7.2f');
labelsAE_F3 = convertCharsToStrings(a);
labelsAE_F3(1);
for i = 2:11
    a = num2str(results_m(6, i), '%7.2f');
    labelsAE_F3(i) = convertCharsToStrings(a);
end

labelsAbs_AE_F3 = "";
a = num2str(results_m(7, 1), '%7.2f');
labelsAbs_AE_F3 = convertCharsToStrings(a);
labelsAbs_AE_F3(1);
for i = 2:11
    a = num2str(results_m(7, i), '%7.2f');
    labelsAbs_AE_F3(i) = convertCharsToStrings(a);
end

labelsRMSE_F3 = "";
a = num2str(results_m(8, 1), '%7.2f');
labelsRMSE_F3 = convertCharsToStrings(a);
labelsRMSE_F3(1);
for i = 2:11
    a = num2str(results_m(8, i), '%7.2f');
    labelsRMSE_F3(i) = convertCharsToStrings(a);
end

labelsAE_F1 = "";
a = num2str(results_m(11, 1), '%7.2f');
labelsAE_F1 = convertCharsToStrings(a);
labelsAE_F1(1);
for i = 2:11
    a = num2str(results_m(11, i), '%7.2f');
    labelsAE_F1(i) = convertCharsToStrings(a);
end

labelsAbs_AE_F1 = "";
a = num2str(results_m(12, 1), '%7.2f');
labelsAbs_AE_F1 = convertCharsToStrings(a);
labelsAbs_AE_F1(1);
for i = 2:11
    a = num2str(results_m(12, i), '%7.2f');
    labelsAbs_AE_F1(i) = convertCharsToStrings(a);
end

labelsRMSE_F1 = "";
a = num2str(results_m(13, 1), '%7.2f');
labelsRMSE_F1 = convertCharsToStrings(a);
labelsRMSE_F1(1);
for i = 2:11
    a = num2str(results_m(13, i), '%7.2f');
    labelsRMSE_F1(i) = convertCharsToStrings(a);
end

labelsAE_F2 = "";
a = num2str(results_m(16, 1), '%7.2f');
labelsAE_F2 = convertCharsToStrings(a);
labelsAE_F2(1);
for i = 2:11
    a = num2str(results_m(16, i), '%7.2f');
    labelsAE_F2(i) = convertCharsToStrings(a);
end

labelsAbs_AE_F2 = "";
a = num2str(results_m(17, 1), '%7.2f');
labelsAbs_AE_F2 = convertCharsToStrings(a);
labelsAbs_AE_F2(1);
for i = 2:11
    a = num2str(results_m(17, i), '%7.2f');
    labelsAbs_AE_F2(i) = convertCharsToStrings(a);
end

labelsRMSE_F2 = "";
a = num2str(results_m(18, 1), '%7.2f');
labelsRMSE_F2 = convertCharsToStrings(a);
labelsRMSE_F2(1);
for i = 2:11
    a = num2str(results_m(18, i), '%7.2f');
    labelsRMSE_F2(i) = convertCharsToStrings(a);
end

%% Plot results

%% Average Error: Overall
subplot(3, 3, [1, 6])
[sorted_data, new_indices] = sort(results_m(1, :)); % sorts in *ascending* order
SortlabelsAE_All = labelsAE_All(new_indices); 
tick_labelAE_All = tick_label(new_indices);
b=bar(sorted_data(3:10), 'FaceColor', '[.7, .7, .7]');
hold on
set(gca, 'XTickLabel', tick_labelAE_All(3:10))
set(gca, 'Ylim', [-120 200]);
ylabel("Average Error (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:8
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -15.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAE_All(3:10), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Overall')

subplot(3, 3, 7)
[sorted_data, new_indices] = sort(results_m(6, :)); % sorts in *ascending* order
SortlabelsAE_F3 = labelsAE_F3(new_indices); 
tick_labelAE_F3 = tick_label(new_indices);
b=bar(sorted_data(3:9), 'b');
b(1).BaseValue = 0;
hold on
set(gca, 'XTickLabel', tick_labelAE_F3(3:9))
ylabel("Average Error (%)")
set(gca, 'Ylim', [-120 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAE_F3(3:9), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Fragments Only')

subplot(3, 3, 8)
[sorted_data, new_indices] = sort(results_m(16, :)); % sorts in *ascending* order
SortlabelsAE_F2 = labelsAE_F2(new_indices); 
tick_labelAE_F2 = tick_label(new_indices);
b=bar(sorted_data(3:9), 'g');
hold on
set(gca, 'XTickLabel', tick_labelAE_F2(3:9))
ylabel("Average Error (%)")
set(gca, 'Ylim', [-120 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAE_F2(3:9), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Films Only')

subplot(3, 3, 9)
[sorted_data, new_indices] = sort(results_m(11, :)); % sorts in *ascending* order
SortlabelsAE_F1 = labelsAE_F1(new_indices); 
tick_labelAE_F1 = tick_label(new_indices);
b=bar(sorted_data(3:9), 'r');
hold on
set(gca, 'XTickLabel', tick_labelAE_F1(3:9))
ylabel("Average Error (%)")
set(gca, 'Ylim', [-120 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1, SortlabelsAE_F1(3:9), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Fibres Only')

sgtitle(sprintf('Evaluation of Models: Average Error (AE) \r\n between Measured and Calculated Velocity'));

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './Output/20220621/AE_Bar7.jpg', 'Resolution', 300)

%% Absolute Average Error: Overall
subplot(3, 3, [1, 6])
[sorted_data, new_indices] = sort(results_m(2, :)); % sorts in *ascending* order
SortlabelsAbs_AE_All = labelsAbs_AE_All(new_indices); 
tick_labelAbs_AE_All = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 7 10]), 'FaceColor', '[.7, .7, .7]');
hold on
set(gca, 'XTickLabel', tick_labelAbs_AE_All([1 2 3 4 5 6 7 10]))
set(gca, 'Ylim', [0 200]);
ylabel("Absolute Average Error (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:8
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -15.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAbs_AE_All([1 2 3 4 5 6 7 10]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Overall')

subplot(3, 3, 7)
[sorted_data, new_indices] = sort(results_m(7, :)); % sorts in *ascending* order
SortlabelsAbs_AE_F3 = labelsAbs_AE_F3(new_indices); 
tick_labelAbs_AE_F3 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 9]), 'b');
b(1).BaseValue = 0;
hold on
set(gca, 'XTickLabel', tick_labelAbs_AE_F3([1 2 3 4 5 6 9]))
ylabel("Absolute Average Error (%)")
set(gca, 'Ylim', [0 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAbs_AE_F3([1 2 3 4 5 6 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Fragments Only')

subplot(3, 3, 8)
[sorted_data, new_indices] = sort(results_m(17, :)); % sorts in *ascending* order
SortlabelsAbs_AE_F2 = labelsAbs_AE_F2(new_indices); 
tick_labelAbs_AE_F2 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 9]), 'g');
hold on
set(gca, 'XTickLabel', tick_labelAbs_AE_F2([1 2 3 4 5 6 9]))
ylabel("Absolute Average Error (%)")
set(gca, 'Ylim', [0 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1,SortlabelsAbs_AE_F2([1 2 3 4 5 6 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Films Only')

subplot(3, 3, 9)
[sorted_data, new_indices] = sort(results_m(12, :)); % sorts in *ascending* order
SortlabelsAbs_AE_F1 = labelsAbs_AE_F1(new_indices); 
tick_labelAbs_AE_F1 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 9]), 'r');
hold on
set(gca, 'XTickLabel', tick_labelAbs_AE_F1([1 2 3 4 5 6 9]))
ylabel("Absolute Average Error (%)")
set(gca, 'Ylim', [0 200]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
for i=1:7
    if (ytips1(i) < 0)
        ytips1(i) = ytips1(i) -30.0;
    else
        ytips1(i) = ytips1(i);
    end
end
text(xtips1,ytips1, SortlabelsAbs_AE_F1([1 2 3 4 5 6 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 8)
title('Fibres Only')

sgtitle(sprintf('Evaluation of Models: Absolute Average Error (|AE|) \r\n between Measured and Calculated Velocity'));

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './Output/20220621/Abs_AE_Bar7.jpg', 'Resolution', 300)

%% RMSE All

subplot(3, 3, [1,6])
[sorted_data, new_indices] = sort(results_m(3, :)); % sorts in *ascending* order
SortlabelsRMSE_All = labelsRMSE_All(new_indices); 
tick_labelRMSE_All = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 7 10]), 'FaceColor', '[.7, .7, .7]');
hold on
set(gca, 'XTickLabel', tick_labelRMSE_All([1 2 3 4 5 6 7 10]))
set(gca, 'Ylim', [0 200]);
ylabel("RMSE (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1,ytips1,SortlabelsRMSE_All([1 2 3 4 5 6 7 10]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Overall')

subplot(3, 3, 7)
[sorted_data, new_indices] = sort(results_m(8, :)); % sorts in *ascending* order
SortlabelsRMSE_F3 = labelsRMSE_F3(new_indices); 
tick_labelRMSE_F3 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 8 9]), 'b');
hold on
set(gca, 'XTickLabel', tick_labelRMSE_F3([1 2 3 4 5 8 9]))
ylabel("RMSE (%)")
set(gca, 'Ylim', [0 250]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1,ytips1,SortlabelsRMSE_F3([1 2 3 4 5 8 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 9)
title('Fragments Only')

subplot(3, 3, 8)
[sorted_data, new_indices] = sort(results_m(18, :)); % sorts in *ascending* order
SortlabelsRMSE_F2 = labelsRMSE_F2(new_indices); 
tick_labelRMSE_F2 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 9]), 'g');
hold on
set(gca, 'XTickLabel', tick_labelRMSE_F2([1 2 3 4 5 6 9]))
set(gca, 'Ylim', [0 250]);
ylabel("RMSE (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1,ytips1,SortlabelsRMSE_F2([1 2 3 4 5 6 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 9)
title('Films Only')

subplot(3, 3, 9)
[sorted_data, new_indices] = sort(results_m(13, :)); % sorts in *ascending* order
SortlabelsRMSE_F1 = labelsRMSE_F1(new_indices); 
tick_labelRMSE_F1 = tick_label(new_indices);
b=bar(sorted_data([1 2 3 4 5 6 9]), 'r');
hold on
set(gca, 'XTickLabel', tick_labelRMSE_F1([1 2 3 4 5 6 9]))
set(gca, 'Ylim', [0 250]);
ylabel("RMSE (%)")
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1,ytips1,SortlabelsRMSE_F1([1 2 3 4 5 6 9]), 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 9)
title('Fibres Only')

sgtitle(sprintf('Evaluation of Models: Root Mean Square Error (RMSE) \r\n between Measured and Calculated Velocity'));

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './Output/20220621/RMSE_Bar7.jpg', 'Resolution', 300)
