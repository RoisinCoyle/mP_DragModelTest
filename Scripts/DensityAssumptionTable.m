%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DensityAssumptionTable
% Date created: 24.05.22
% Date last mostified: 25.05.22
% Purpose: Generate table for density assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data
clear
VM_Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

%% Extract 3 random particles from each particle type

rng('default')

frag = randi([1 80], 1, 3);
fibre = randi([81 100], 1, 3);
film = randi([101 140], 1, 3);

frag_1 = VM_Dataset(frag(1), :);
frag_2 = VM_Dataset(frag(2), :);
frag_3 = VM_Dataset(frag(3), :);

fibre_1 = VM_Dataset(fibre(1), :);
fibre_2 = VM_Dataset(fibre(2), :);
fibre_3 = VM_Dataset(fibre(3), :);

film_1 = VM_Dataset(film(1), :);
film_2 = VM_Dataset(film(2), :);
film_3 = VM_Dataset(film(3), :);

%% Construct new table

test_density = [1019.8, 1023.6, 1027.4, 1035.0, 1042.6, 1050.3];
test_visdyn = [0.948e-3, 0.959e-3, 0.970e-3, 0.993e-3, 1.016e-3, 1.041e-3];
test_viskin = zeros(1, 6);
for i=1:6
    test_viskin(i) = test_visdyn(i) / test_density(i);
end

frag_all_1 = [frag_1; frag_1; frag_1; frag_1; frag_1; frag_1];
frag_all_1.FluidDensity(:) = test_density;
frag_all_1.DynamicViscosity(:) = test_visdyn;
frag_all_1.KinematicViscosity(:) = test_viskin;
frag_all_1 = removevars(frag_all_1, ["CdMeasured", "Re", "Wmeasured"]);

frag_all_2 = [frag_2; frag_2; frag_2; frag_2; frag_2; frag_2];
frag_all_2.FluidDensity(:) = test_density;
frag_all_2.DynamicViscosity(:) = test_visdyn;
frag_all_2.KinematicViscosity(:) = test_viskin;
frag_all_2 = removevars(frag_all_2, ["CdMeasured", "Re", "Wmeasured"]);

frag_all_3 = [frag_3; frag_3; frag_3; frag_3; frag_3; frag_3];
frag_all_3.FluidDensity(:) = test_density;
frag_all_3.DynamicViscosity(:) = test_visdyn;
frag_all_3.KinematicViscosity(:) = test_viskin;
frag_all_3 = removevars(frag_all_3, ["CdMeasured", "Re", "Wmeasured"]);

fibre_all_1 = [fibre_1; fibre_1; fibre_1; fibre_1; fibre_1; fibre_1];
fibre_all_1.FluidDensity(:) = test_density;
fibre_all_1.DynamicViscosity(:) = test_visdyn;
fibre_all_1.KinematicViscosity(:) = test_viskin;
fibre_all_1 = removevars(fibre_all_1, ["CdMeasured", "Re", "Wmeasured"]);

fibre_all_2 = [fibre_2; fibre_2; fibre_2; fibre_2; fibre_2; fibre_2];
fibre_all_2.FluidDensity(:) = test_density;
fibre_all_2.DynamicViscosity(:) = test_visdyn;
fibre_all_2.KinematicViscosity(:) = test_viskin;
fibre_all_2 = removevars(fibre_all_2, ["CdMeasured", "Re", "Wmeasured"]);

fibre_all_3 = [fibre_3; fibre_3; fibre_3; fibre_3; fibre_3; fibre_3];
fibre_all_3.FluidDensity(:) = test_density;
fibre_all_3.DynamicViscosity(:) = test_visdyn;
fibre_all_3.KinematicViscosity(:) = test_viskin;
fibre_all_3 = removevars(fibre_all_3, ["CdMeasured", "Re", "Wmeasured"]);

film_all_1 = [film_1; film_1; film_1; film_1; film_1; film_1];
film_all_1.FluidDensity(:) = test_density;
film_all_1.DynamicViscosity(:) = test_visdyn;
film_all_1.KinematicViscosity(:) = test_viskin;
film_all_1 = removevars(film_all_1, ["CdMeasured", "Re", "Wmeasured"]);

film_all_2 = [film_2; film_2; film_2; film_2; film_2; film_2];
film_all_2.FluidDensity(:) = test_density;
film_all_2.DynamicViscosity(:) = test_visdyn;
film_all_2.KinematicViscosity(:) = test_viskin;
film_all_2 = removevars(film_all_2, ["CdMeasured", "Re", "Wmeasured"]);

film_all_3 = [film_3; film_3; film_3; film_3; film_3; film_3];
film_all_3.FluidDensity(:) = test_density;
film_all_3.DynamicViscosity(:) = test_visdyn;
film_all_3.KinematicViscosity(:) = test_viskin;
film_all_3 = removevars(film_all_3, ["CdMeasured", "Re", "Wmeasured"]);

DensityTestTable = [frag_all_1; frag_all_2; frag_all_3; ...
                    fibre_all_1; fibre_all_2; fibre_all_3; ... 
                    film_all_1; film_all_2; film_all_3];

%% Output table
writetable(DensityTestTable, './SettlingVelocity calc/DensityTestTable.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(DensityTestTable, './SettlingVelocity calc/DensityTestTable.xls', 'WriteRowNames', true);

writetable(DensityTestTable, './DragModelsTest/Output/20220517/DensityTestTable.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(DensityTestTable, './DragModelsTest/Output/20220517/DensityTestTable.xls', 'WriteRowNames', true);

