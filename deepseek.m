% Radar Cross Section (RCS) Benchmark for Aircraft Geometries
% Searches for STL files in a subfolder and compares RCS performance
% Uses Method of Moments (MoM) and other full-wave solvers

close all;
clear;
clc;

%% Parameters Setup
freq = 10e9; % Frequency (10 GHz)
lambda = physconst('lightspeed') / freq;
azimuth = 0:1:180; % Azimuth angles (deg)
elevation = 0; % Elevation angle (deg)

% Solver settings
solverTypes = {'MoM', 'FMM', 'MLFMM'}; % Full-wave solvers to compare
polarizations = {'HH', 'VV'}; % Polarizations to analyze

%% Locate STL files in subfolder
subfolder = 'geometries'; % Name of your subfolder
stlFiles = dir(fullfile(subfolder, '*.stl')); % Find all STL files

if isempty(stlFiles)
    error('No STL files found in the %s subfolder', subfolder);
end

numAircraft = min(length(stlFiles), 3); % Use up to 3 aircraft
fprintf('Found %d aircraft models in %s folder\n', numAircraft, subfolder);

%% Load and Prepare Aircraft Geometries
platforms = cell(1, numAircraft);
aircraftNames = cell(1, numAircraft);

for i = 1:numAircraft
    p = platform;
    p.FileName = fullfile(subfolder, stlFiles(i).name);
    p.Units = 'm';
    
    % Mesh refinement (adjust based on your geometry complexity)
    fprintf('Meshing %s...\n', stlFiles(i).name);
    try
        % Attempt fine mesh
        mesh(p, 'MaxEdgeLength', lambda/5); % Fine mesh for accurate MoM
    catch ME
        warning('Meshing failed for %s: %s', stlFiles(i).name, ME.message);
        try
            % Attempt coarser mesh
            fprintf('Retrying with coarser mesh...\n');
            mesh(p, 'MaxEdgeLength', lambda); % Coarser mesh
        catch ME2
            warning('Coarser meshing also failed for %s: %s', stlFiles(i).name, ME2.message);
            continue; % Skip this geometry if meshing fails
        end
    end
    
    platforms{i} = p;
    aircraftNames{i} = stlFiles(i).name;
    
    % Display the geometry
    figure;
    show(p);
    title(['Geometry: ' strrep(stlFiles(i).name, '_', '\_')]);
end

%% RCS Calculation and Comparison
results = struct();

for solverIdx = 1:length(solverTypes)
    solver = solverTypes{solverIdx};
    
    for polIdx = 1:length(polarizations)
        pol = polarizations{polIdx};
        
        fprintf('\nCalculating with %s solver, %s polarization...\n', solver, pol);
        
        % Calculate RCS for each aircraft
        for aircraftIdx = 1:numAircraft
            p = platforms{aircraftIdx};
            fprintf('Processing %s...\n', aircraftNames{aircraftIdx});
            
            try
                % Try with GPU acceleration if available
                fprintf('Attempting RCS calculation with GPU...\n');
                sigma = rcs(p, freq, azimuth, elevation, ...
                          'Solver', solver, ...
                          'Polarization', pol, ...
                          'EnableGPU', true);
            catch ME
                fprintf('GPU failed, falling back to CPU: %s\n', ME.message);
                try
                    % Fallback to CPU
                    fprintf('Attempting RCS calculation with CPU...\n');
                    sigma = rcs(p, freq, azimuth, elevation, ...
                              'Solver', solver, ...
                              'Polarization', pol, ...
                              'EnableGPU', false);
                catch ME2
                    warning('RCS calculation failed for %s: %s', ...
                           aircraftNames{aircraftIdx}, ME2.message);
                    sigma = NaN(size(azimuth));
                end
            end
            
            % Store results
            results(aircraftIdx).(solver).(pol) = sigma;
            results(aircraftIdx).name = aircraftNames{aircraftIdx};
        end
    end
end

%% Performance Analysis
% Calculate average RCS for each aircraft (lower is better for stealth)
performance = zeros(numAircraft, length(solverTypes), length(polarizations));

for aircraftIdx = 1:numAircraft
    for solverIdx = 1:length(solverTypes)
        for polIdx = 1:length(polarizations)
            solver = solverTypes{solverIdx};
            pol = polarizations{polIdx};
            
            if isfield(results(aircraftIdx), solver) && ...
               isfield(results(aircraftIdx).(solver), pol)
                rcsData = results(aircraftIdx).(solver).(pol);

                % Validate RCS data
                if isempty(rcsData) || ~isnumeric(rcsData)
                    warning('Invalid RCS data for aircraft %d, solver %s, polarization %s', ...
                            aircraftIdx, solver, pol);
                    performance(aircraftIdx, solverIdx, polIdx) = NaN;
                    continue;
                end

                % Convert to linear scale, average, then back to dB
                linearRCS = 10.^(rcsData/10);
                meanLinearRCS = mean(linearRCS, 'omitnan');
                performance(aircraftIdx, solverIdx, polIdx) = 10*log10(meanLinearRCS);
            else
                performance(aircraftIdx, solverIdx, polIdx) = NaN;
            end
        end
    end
end

% Handle missing performance data
if isempty(performance) || all(isnan(performance(:)))
    warning('Performance data is empty or contains only NaN values.');
    return; % Exit gracefully
end

% Calculate overall performance
overallPerformance = squeeze(mean(performance, [2,3], 'omitnan'));

% Find the best performing aircraft
[bestRCS, bestIdx] = min(overallPerformance);

fprintf('\nBest performing aircraft: %s (Average RCS: %.2f dBsm)\n', ...
       strrep(aircraftNames{bestIdx}, '.stl', ''), bestRCS);

%% Visualization
% Plot RCS patterns for each aircraft
for aircraftIdx = 1:numAircraft
    figure;
    hold on;
    titleStr = sprintf('RCS Comparison for %s\nFrequency: %.2f GHz', ...
                      strrep(aircraftNames{aircraftIdx}, '_', '\_'), freq/1e9);
    
    numLegendEntries = length(solverTypes) * length(polarizations);
    legendEntries = cell(1, numLegendEntries);
    lineStyles = {'-', '--', ':', '-.'};
    colors = lines(length(solverTypes)*length(polarizations));
    colorIdx = 1;
    entryIdx = 1; % Initialize index
    for solverIdx = 1:length(solverTypes)
        for polIdx = 1:length(polarizations)
            solver = solverTypes{solverIdx};
            pol = polarizations{polIdx};

            if ~isempty(results(aircraftIdx).(solver).(pol)) && ...
               length(azimuth) == length(results(aircraftIdx).(solver).(pol))
                plot(azimuth, results(aircraftIdx).(solver).(pol), ...
                    'LineStyle', lineStyles{solverIdx}, ...
                    'Color', colors(colorIdx,:), ...
                    'LineWidth', 2);
                legendEntries{end+1} = sprintf('%s, %s', solver, pol);
                                colorIdx = colorIdx + 1;
            else
                warning('Skipping plot for aircraft %d, solver %s, polarization %s due to invalid data.', ...
                        aircraftIdx, solver, pol);
            end
        end
    end
    
    title(titleStr);
    xlabel('Azimuth Angle (deg)');
    ylabel('RCS (dBsm)');
    legend(legendEntries, 'Location', 'best');
    grid on;
    hold off;
end

%% Display Results
fprintf('\n=== Average RCS Performance Comparison (dBsm) ===\n');
fprintf('Aircraft\t\tSolver\t\tHH-pol\tVV-pol\n');

for aircraftIdx = 1:numAircraft
    [~, nameOnly, ~] = fileparts(aircraftNames{aircraftIdx});
    for solverIdx = 1:length(solverTypes)
        solver = solverTypes{solverIdx};
        
        if aircraftIdx > 1 && solverIdx == 1
            fprintf('%s\t', nameOnly);
        elseif solverIdx == 1
            fprintf('%s\t', nameOnly);
        else
            fprintf('\t');
        end
        
        fprintf('%s\t', solver);
        
        for polIdx = 1:length(polarizations)
            val = performance(aircraftIdx, solverIdx, polIdx);
            if isnan(val)
                fprintf('N/A\t');
            else
                fprintf('%.2f\t', val);
            end
        end
        fprintf('\n');
    end
end

%% Determine Best Performing Aircraft
% Calculate overall average across all solvers and polarizations
overallPerformance = squeeze(mean(performance, [2,3], 'omitnan'));
[bestRCS, bestIdx] = min(overallPerformance);

fprintf('\nBest performing aircraft: %s (Average RCS: %.2f dBsm)\n', ...
       strrep(aircraftNames{bestIdx}, '.stl', ''), bestRCS);

%% Create Summary Figure
figure;
hold on;
barData = squeeze(mean(performance, 3, 'omitnan'))'; % Average across polarizations
bar(barData);
set(gca, 'XTickLabel', solverTypes);
ylabel('Average RCS (dBsm)');
title('Aircraft RCS Performance Comparison');
legend(strrep(aircraftNames, '.stl', ''), 'Location', 'best');
grid on;
hold off;
