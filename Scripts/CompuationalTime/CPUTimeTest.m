%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: CPUTimeTest.m
% Date created: 15.02.22
% Date last mostified: 15.02.22
% Purpose: To provide evidence for the thesis that the implicit models are
%          more computationally expensive than the implicit models (not
%          included in paper, only in thesis).
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Function to be used: cputime

% Description: Returns the total CPU time used by Matlab since it was
% started. The returned CPU time is expressed in seconds.

% Each call to cputime returns the total CPU time used by Matlab up to the
% point when the function is called. To measure the CPU time used to run
% your code, place two calls to cputime before and after the code and then
% calculate the difference between the returned values.

% This differs from the timeit or tic and toc functions which measure the
% wall-clock time needed to run the programme.

% It is recommended that you use timeit or tic and toc to measure the
% performance of your code. These functions return wall-clock time. Unlike
% tic and toc, the timeit function calls your code multiple times and
% therefore considers first-time costs. The cputime function measures the
% total CPU time and sums across all threads. 

% It is recommended that you put the code you are trying to time into a
% functoin instead of timing it at the command line or within the script.

% Unless you are trying to measure the first time cost, you should run your
% code multiple times i.e., use the timeit function.

%% Read in data:
% Tests were carried out for each model to measure the time taken for the
% calculation to be completed for all particles in Van Melkebeke's dataset
% (140 total) 1000 times. This was replicated three times for each model,
% giving a total of 3000 measurements.

clear

% Bagheri

BB_1 = readtable("DragModelsTest\Output\20230215\BB_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
BB_1.Properties.VariableNames = {'BB_SA_1', 'BB_Proj_1'};
BB_2 = readtable("DragModelsTest\Output\20230215\BB_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
BB_2.Properties.VariableNames = {'BB_SA_2', 'BB_Proj_2'};
BB_3 = readtable("DragModelsTest\Output\20230215\BB_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
BB_3.Properties.VariableNames = {'BB_SA_3', 'BB_Proj_3'};

BB_All = [BB_1, BB_2, BB_3]

% Dietrich

Dietrich_1 = readtable("DragModelsTest\Output\20230215\Dietrich_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dietrich_1.Properties.VariableNames = {'Dietrich_1'};
Dietrich_2 = readtable("DragModelsTest\Output\20230215\Dietrich_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dietrich_2.Properties.VariableNames = {'Dietrich_2'};
Dietrich_3 = readtable("DragModelsTest\Output\20230215\Dietrich_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dietrich_3.Properties.VariableNames = {'Dietrich_3'};

Dietrich_All = [Dietrich_1, Dietrich_2, Dietrich_3]

% Dioguardi

Dio_1 = readtable("DragModelsTest\Output\20230215\Dio_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dio_1.Properties.VariableNames = {'Dio_SA_1', 'Dio_Proj_1'};
Dio_2 = readtable("DragModelsTest\Output\20230215\Dio_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dio_2.Properties.VariableNames = {'Dio_SA_2', 'Dio_Proj_2'};
Dio_3 = readtable("DragModelsTest\Output\20230215\Dio_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Dio_3.Properties.VariableNames = {'Dio_SA_3', 'Dio_Proj_3'};

Dio_All = [Dio_1, Dio_2, Dio_3]

% Francalanci

Frn_1 = readtable("DragModelsTest\Output\20230215\Frn_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Frn_1.Properties.VariableNames = {'Frn_1'};
Frn_2 = readtable("DragModelsTest\Output\20230215\Frn_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Frn_2.Properties.VariableNames = {'Frn_2'};
Frn_3 = readtable("DragModelsTest\Output\20230215\Frn_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Frn_3.Properties.VariableNames = {'Frn_3'};

Frn_All = [Frn_1, Frn_2, Frn_3]

% Stokes

Stokes_1 = readtable("DragModelsTest\Output\20230215\Stokes_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Stokes_1.Properties.VariableNames = {'Stokes_SA_1', 'Stokes_Proj_1'};
Stokes_2 = readtable("DragModelsTest\Output\20230215\Stokes_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Stokes_2.Properties.VariableNames = {'Stokes_SA_2', 'Stokes_Proj_2'};
Stokes_3 = readtable("DragModelsTest\Output\20230215\Stokes_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Stokes_3.Properties.VariableNames = {'Stokes_SA_3', 'Stokes_Proj_3'};

Stokes_All = [Stokes_1, Stokes_2, Stokes_3]

% Yu

Yu_1 = readtable("DragModelsTest\Output\20230215\Yu_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Yu_1.Properties.VariableNames = {'Yu_1'};
Yu_2 = readtable("DragModelsTest\Output\20230215\Yu_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Yu_2.Properties.VariableNames = {'Yu_2'};
Yu_3 = readtable("DragModelsTest\Output\20230215\Yu_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
Yu_3.Properties.VariableNames = {'Yu_3'};

Yu_All = [Yu_1, Yu_2, Yu_3]

% Zhang

ZC_1 = readtable("DragModelsTest\Output\20230215\ZC_Timing1.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
ZC_1.Properties.VariableNames = {'ZC_SA_1', 'ZC_Proj_1'};
ZC_2 = readtable("DragModelsTest\Output\20230215\ZC_Timing2.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
ZC_2.Properties.VariableNames = {'ZC_SA_2', 'ZC_Proj_2'};
ZC_3 = readtable("DragModelsTest\Output\20230215\ZC_Timing3.txt", Delimiter=',', ReadRowNames=true, ReadVariableNames=true);
ZC_3.Properties.VariableNames = {'ZC_SA_3', 'ZC_Proj_3'};

ZC_All = [ZC_1, ZC_2, ZC_3]

%% Combine all into one table for plotting

CPU_Table_All = [BB_All, Dietrich_All, Dio_All, Frn_All, Stokes_All, Yu_All, ZC_All];

CPU_Table_All = rows2vars(CPU_Table_All);
CPU_Table_All = sortrows(CPU_Table_All, "OriginalVariableNames", "ascend");

%% Plot as 33 individual bars
TickLabels = CPU_Table_All.OriginalVariableNames

b=bar(CPU_Table_All.AverageCPUPerTest(1:33))
title(sprintf("Average CPU Time taken to calculate the settling velocity of all particles \n in Van Melkebeke's dataset"))
ylabel('CPU Time (seconds)')
set(gca, 'XTick', 1:33)
set(gca, 'XTickLabel', TickLabels)
b.FaceColor = 'flat';
for i=1:6
    b.CData([i], :) = [0.4, 0.7608, 0.6471];
end
for i = 7:9
    b.CData([i], :) = [0.9882, 0.5529, 0.6843];
end
for i=10:15
    b.CData([i], :) = [0.5529, 0.6275, 0.7961];
end
for i=16:18
    b.CData([i], :) = [0.9059, 0.5412,7647];
end
for i=19:24
    b.CData([i], :) = [0.6510, 0.8471, 0.3294];
end
for i=25:27
    b.CData([i], :) = [1, 0.8510, 0.1843];
end
for i=28:33
    b.CData([i], :) = [0.8980, 0.7686, 0.5804];
end

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230215/CPU33.jpg', 'Resolution', 300)
%% Plot as 11 groups of 3 bars

grouped = zeros(11, 3);
for i=1:11
    for j=1:3
        grouped(i, j) = CPU_Table_All.AverageCPUPerTest((3*i)-(3-j));
    end
end

b=bar(grouped, 'FaceColor', 'flat')
title(sprintf("Average CPU Time taken to calculate the settling velocity of all particles \n in Van Melkebeke's dataset"))
ylabel('CPU Time (seconds)')
set(gca, 'XTick', 1:33)
set(gca, 'XTickLabel', {'BB Proj', 'BB SA', 'Dietrich', 'Dio Proj', 'Dio SA', 'Frn', 'Stokes Proj', 'Stokes SA', 'Yu', 'ZC Proj', 'ZC SA'})

Colours(1, :) = [0.4, 0.7608, 0.6471];
Colours(2, :) = [0.4, 0.7608, 0.6471];
Colours(3, :) = [0.9882, 0.5529, 0.6843];
Colours(4, :) = [0.5529, 0.6275, 0.7961];
Colours(5, :) = [0.5529, 0.6275, 0.7961];
Colours(6, :) = [0.9059, 0.5412,7647];
Colours(7, :) = [0.6510, 0.8471, 0.3294];
Colours(8, :) = [0.6510, 0.8471, 0.3294];
Colours(9, :) = [1, 0.8510, 0.1843];
Colours(10, :) = [0.8980, 0.7686, 0.5804];
Colours(11, :) = [0.8980, 0.7686, 0.5804];
b(1).CData = Colours
b(2).CData = Colours;
b(3).CData = Colours;

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230215/CPU11by3.jpg', 'Resolution', 300)

%% Plot bars as 3 groups of 11

grouped = zeros(3, 11)
for i=1:3
    for j=1:11
        if i==1
            grouped(i, j) = CPU_Table_All.AverageCPUPerTest((i*3)+((3*j)-5));
        elseif i==2
            grouped(i, j) = CPU_Table_All.AverageCPUPerTest((i*3)+((3*j)-7));
        elseif i==3
            grouped(i, j) = CPU_Table_All.AverageCPUPerTest((i*3)+((3*j)-9));
        end
    end
end

Colours1 = zeros(3, 3);
Colours2 = zeros(3, 3);
Colours3 = zeros(3, 3);
Colours4 = zeros(3, 3);
Colours5 = zeros(3, 3);
Colours6 = zeros(3, 3);
Colours7 = zeros(3, 3);
Colours8 = zeros(3, 3);
Colours9 = zeros(3, 3);
Colours10 = zeros(3, 3);
Colours11 = zeros(3, 3);

for i=1:3
    Colours1(i, :) = [0.4, 0.7608, 0.6471];
    Colours2(i, :) = [0.4, 0.7608, 0.6471];
    Colours3(i, :) = [0.9882, 0.5529, 0.6843];
    Colours4(i, :) = [0.5529, 0.6275, 0.7961];
    Colours5(i, :) = [0.5529, 0.6275, 0.7961];
    Colours6(i, :) = [0.9059, 0.5412,7647];
    Colours7(i, :) = [0.6510, 0.8471, 0.3294];
    Colours8(i, :) = [0.6510, 0.8471, 0.3294];
    Colours9(i, :) = [1, 0.8510, 0.1843];
    Colours10(i, :) = [0.8980, 0.7686, 0.5804];
    Colours11(i, :) = [0.8980, 0.7686, 0.5804];
end

b=bar(grouped);
title(sprintf("Average CPU Time taken to calculate the settling velocity of all particles \n in Van Melkebeke's dataset"))
ylabel('CPU Time (seconds)')
set(gca, 'XTick', 1:33)
set(gca, 'XTickLabel', {'1^{st} Run', '2^{nd} Run', '3^{rd} Run'})
for k = 1:11
    b(k).FaceColor = 'flat';
end

b(1).CData = Colours1;
b(2).CData = Colours2;
b(3).CData = Colours3;
b(4).CData = Colours4;
b(5).CData = Colours5;
b(6).CData = Colours6;
b(7).CData = Colours7;
b(8).CData = Colours8;
b(9).CData = Colours9;
b(10).CData = Colours10;
b(11).CData = Colours11;

legend('BB Proj', 'BB SA', 'Dietrich', 'Dio Proj', 'Dio SA', 'Frn', 'Stokes Proj', 'Stokes SA', 'Yu', 'ZC Proj', 'ZC SA')

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230215/CPU3by11.jpg', 'Resolution', 300)

%% Plot bar graph comparing CPU of BB_Proj, Dio_Proj and Yu

BBDIOYU(1, :)= CPU_Table_All(1, :)
BBDIOYU(2, :)= CPU_Table_All(2, :)
BBDIOYU(3, :)= CPU_Table_All(3, :)
BBDIOYU(4, :)= CPU_Table_All(10, :)
BBDIOYU(5, :)= CPU_Table_All(11, :)
BBDIOYU(6, :)= CPU_Table_All(12, :)
BBDIOYU(7, :)= CPU_Table_All(25, :)
BBDIOYU(8, :)= CPU_Table_All(26, :)
BBDIOYU(9, :)= CPU_Table_All(27, :)

groupedBDY = zeros(3, 3);
for i=1:3
    for j=1:3
        if i==1
            groupedBDY(i, j) = CPU_Table_All.AverageCPUPerTest(i);
        elseif i==2
            groupedBDY(i, j) = CPU_Table_All.AverageCPUPerTest((i*2)+((j-1)));
        elseif i==3
            groupedBDY(i, j) = CPU_Table_All.AverageCPUPerTest((i*2)+j);
        end
    end
end

b=bar(groupedBDY, 'FaceColor', 'flat')
title(sprintf("Average CPU Time taken to calculate the settling velocity of all particles \n in Van Melkebeke's dataset"))
ylabel('CPU Time (seconds)')
set(gca, 'XTick', 1:33)
set(gca, 'XTickLabel', {'Bagheri and Bonadonna', 'Dioguardi et al.', 'Yu et al'})


ColoursBDY(1, :) = [0.6510, 0.8078, 0.8902];
ColoursBDY(2, :) = [0.12116, 0.4706, 0.7059];
ColoursBDY(3, :) = [0.6980, 0.8745, 0.5412];
b(1).CData = ColoursBDY;
b(2).CData = ColoursBDY;
b(3).CData = ColoursBDY;

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230215/CPU_BBDioYu.jpg', 'Resolution', 300)