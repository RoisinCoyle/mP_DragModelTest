%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DensityAssumptionTable
% Date created: 22.07.22
% Date last mostified: 22.07.22
% Purpose: Generate table for initial velocity assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data
clear
VM_Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

%% Extract 2 random particles from each particle type

rng('default')

frag = randi([1 80], 1, 2);
fibre = randi([81 100], 1, 2);
film = randi([101 140], 1, 2);

frag_1 = VM_Dataset(frag(1), :);
frag_2 = VM_Dataset(frag(2), :);

fibre_1 = VM_Dataset(fibre(1), :);
fibre_2 = VM_Dataset(fibre(2), :);

film_1 = VM_Dataset(film(1), :);
film_2 = VM_Dataset(film(2), :);

%% Construct new table

initial_velocity = [0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001];

frag_all_1 = [frag_1; frag_1; frag_1; frag_1; frag_1; frag_1];
frag_all_1.Initial_w(:) = initial_velocity;
frag_all_1 = removevars(frag_all_1, ["CdMeasured", "Re", "Wmeasured"]);

frag_all_2 = [frag_2; frag_2; frag_2; frag_2; frag_2; frag_2];
frag_all_2.Initial_w(:) = initial_velocity;
frag_all_2 = removevars(frag_all_2, ["CdMeasured", "Re", "Wmeasured"]);

fibre_all_1 = [fibre_1; fibre_1; fibre_1; fibre_1; fibre_1; fibre_1];
fibre_all_1.Initial_w(:) = initial_velocity;
fibre_all_1 = removevars(fibre_all_1, ["CdMeasured", "Re", "Wmeasured"]);

fibre_all_2 = [fibre_2; fibre_2; fibre_2; fibre_2; fibre_2; fibre_2];
fibre_all_2.Initial_w(:) = initial_velocity;
fibre_all_2 = removevars(fibre_all_2, ["CdMeasured", "Re", "Wmeasured"]);

film_all_1 = [film_1; film_1; film_1; film_1; film_1; film_1];
film_all_1.Initial_w(:) = initial_velocity;
film_all_1 = removevars(film_all_1, ["CdMeasured", "Re", "Wmeasured"]);

film_all_2 = [film_2; film_2; film_2; film_2; film_2; film_2];
film_all_2.Initial_w(:) = initial_velocity;
film_all_2 = removevars(film_all_2, ["CdMeasured", "Re", "Wmeasured"]);

VelocityTestTable = [frag_all_1; frag_all_2; ...
                    fibre_all_1; fibre_all_2; ... 
                    film_all_1; film_all_2];

%% If the mP density is less than the fluid density, change the mP density
% to approximately the density of polyester to ensure that the particle is
% going to sink.

for i=1:36
    if (VelocityTestTable.ParticleDensity(i) < VelocityTestTable.FluidDensity(i))
        VelocityTestTable.ParticleDensity(i) = 1400;
    end
end


%% Output table
writetable(VelocityTestTable, './SettlingVelocity calc/VelocityTestTableNew.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(VelocityTestTable, './SettlingVelocity calc/VelocityTestTableNew.xls', 'WriteRowNames', true);

writetable(VelocityTestTable, './DragModelsTest/Output/20220621/VelocityTestTableNew.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(VelocityTestTable, './DragModelsTest/Output/20220621/VelocityTestTableNew.xls', 'WriteRowNames', true);

