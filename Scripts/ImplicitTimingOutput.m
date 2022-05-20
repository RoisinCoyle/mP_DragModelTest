%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: ImplictTimes
% Date created: 25.04.22
% Date last mostified: 25.04.22
% Purpose: Output of timing data summary for implicit models
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Table_Stokes_SA= readtable("./DragModelsTest/Output/StokesOutputVM_SA.txt", "Delimiter", ",");
Table_Stokes_Proj= readtable("./DragModelsTest/Output/StokesOutputVM_Proj.txt", "Delimiter", ",");

MaxT_StokesSA = max(Table_Stokes_SA.Time);
MaxT_StokesProj = max(Table_Stokes_Proj.Time);

Table_BB_SA= readtable("./DragModelsTest/Output/BagheriOutputVM_SA.txt", "Delimiter", ",");
Table_BB_Proj= readtable("./DragModelsTest/Output/BagheriOutputVM_Proj.txt", "Delimiter", ",");

MaxT_BBSA = max(Table_BB_SA.Time);
MaxT_BBProj = max(Table_BB_Proj.Time);

Table_Dio_SA= readtable("./DragModelsTest/Output/DioguardiOutputVM_SA.txt", "Delimiter", ",");
Table_Dio_Proj= readtable("./DragModelsTest/Output/DioguardiOutputVM_SA.txt", "Delimiter", ",");

MaxT_DioSA = max(Table_Dio_SA.Time);
MaxT_DioProj = max(Table_Dio_Proj.Time);

Table_Zhang_SA= readtable("./DragModelsTest/Output/ZhangOutputVM_SA.txt", "Delimiter", ",");
Table_Zhang_Proj= readtable("./DragModelsTest/Output/ZhangOutputVM_Proj.txt", "Delimiter", ",");

MaxT_ZhangSA = max(Table_Stokes_SA.Time);
MaxT_ZhangProj = max(Table_Stokes_Proj.Time);

MaxTime = zeros(4, 4);
MaxTime(1, 1) = MaxT_StokesSA;
MaxTime(1, 2) = 0.0002;
MaxTime(1, 3) = MaxT_StokesProj;
MaxTime(1, 4) = 0.0005;
MaxTime(2, 1) = MaxT_BBSA;
MaxTime(2, 2) = 0.0002;
MaxTime(2, 3) = MaxT_BBProj;
MaxTime(2, 4) = 0.0002;
MaxTime(3, 1) = MaxT_DioSA;
MaxTime(3, 2) = 0.0002;
MaxTime(3, 3) = MaxT_DioProj;
MaxTime(3, 4) = 0.0002;
MaxTime(4, 1) = MaxT_ZhangSA;
MaxTime(4, 2) = 0.0002;
MaxTime(4, 3) = MaxT_ZhangProj;
MaxTime(4, 4) = 0.0002;

Table_Implicit_Time = array2table(MaxTime, "VariableNames", ...
    {'MaxTimeSA', 'TimestepSA', 'MaxTimeProj', 'TimestepProj'});

Names = {'Stokes' 'Bagheri' 'Dioguardi' 'Zhang'};
Names = reshape(Names, [4, 1]);
Names_T = cell2table(Names);

Table_Implicit_Time = [Names_T Table_Implicit_Time];
Table_Implicit_Time.Properties.VariableNames(1) = {'Model'};

writetable(Table_Implicit_Time, './DragModelsTest/Output/ImplicitTimes.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Implicit_Time, './DragModelsTest/Output/ImplicitTimes.xls', 'WriteRowNames', true);
