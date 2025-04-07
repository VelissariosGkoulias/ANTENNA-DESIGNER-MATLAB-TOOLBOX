% MATLAB Script for Antenna Designer Toolbox - STL Import and RCS Calculation

% This script demonstrates how to import custom antenna geometries from STL
% files using the Antenna Designer Toolbox and calculate their Radar Cross
% Section (RCS).

clear;
close all;
clc;

% --- 1. Define Simulation Parameters ---
frequency = 1e9;       % Operating frequency (e.g., 1 GHz)
incidentAnglePhi = 0;  % Incident angle in phi (azimuth)
incidentAngleTheta = 0; % Incident angle in theta (elevation)

polarizations = {'HH', 'VV', 'HV', 'VH'};
for p = 1:length(polarizations)
    polarization = polarizations{p};
    disp(['--- Starting Polarization: ', polarization, ' ---']);

    % --- 2. Specify STL File Paths ---
    stlFiles = {
        'geometries/aircraft_baseline.stl'; % Replace with your baseline aircraft STL file
        'geometries/aircraft_modified1.stl'; % Replace with your first modified geometry
        'geometries/aircraft_modified2.stl'; % Replace with your second modified geometry
        % Add more STL file paths as needed
    };

    % --- 3. Loop Through STL Files and Calculate RCS ---
    numGeometries = length(stlFiles);
    rcsValues = cell(numGeometries, 1);
    geometryObjects = cell(numGeometries, 1);

    for i = 1:numGeometries
        currentFile = stlFiles{i};
        fprintf('\n--- Processing Geometry: %s ---\n', currentFile);

        % Check if the file exists
        if ~isfile(currentFile)
            error('File not found: %s', currentFile);
        end

        try
            % Import geometry from STL file
            geometry = stlread(currentFile);

            % Create a custom antenna geometry object
            customAntenna = customAntennaGeometry(); % Add parentheses if it is a class
            customAntenna.Faces = geometry.ConnectivityList;
            customAntenna.Vertices = geometry.Points;

            if isempty(customAntenna.Faces) || isempty(customAntenna.Vertices)
                error('Invalid geometry for file: %s', currentFile);
            end

            % Example: Refine the mesh of the custom antenna geometry
            customAntenna = mesh(customAntenna, 'MaxEdgeLength', 0.01); % Adjust edge length as needed

            % Store the geometry object (optional, for visualization later)
            geometryObjects{i} = customAntenna;

            % Create a radar problem object using the Method of Moments (MoM) solver
            % This object is used to calculate the RCS of the target geometry
            radar = radarProblem('Solver', 'MoM', 'Frequency', frequency);
            if isempty(radar)
                error('Failed to create radar problem object.');
            end

            % Add the custom antenna as a target
            target = customAntenna;
            add(radar, target);

            % Display radar problem object
            disp(radar);

            % Display custom antenna object
            disp(customAntenna);

            try
                % Validate geometry
                if isempty(customAntenna.Faces) || isempty(customAntenna.Vertices)
                    error('Invalid geometry for file: %s', currentFile);
                end

                % Validate radar object
                radar = radarProblem('Solver', 'MoM', 'Frequency', frequency);
                if isempty(radar)
                    error('Failed to create radar problem object.');
                end

                % Validate input parameters
                fprintf('IncidentAngles: Phi=%.1f, Theta=%.1f\n', incidentAnglePhi, incidentAngleTheta);
                fprintf('Polarization: %s\n', polarization);

                % Calculate the monostatic RCS
                [sigma, patterninfo] = rcs(radar, frequency, ...
                    'IncidentAngles', [incidentAnglePhi incidentAngleTheta], ...
                    'Polarization', polarization);

                % Store the RCS value (in dBsm)
                if (~isempty(sigma) && isnumeric(sigma) && sigma > 0)
                    rcs_dbsm = 10*log10(abs(sigma));
                    rcsValues{i} = rcs_dbsm;
                else
                    warning('Invalid RCS value for file: %s', currentFile);
                    rcsValues{i} = NaN;
                end

                fprintf('  RCS at (Phi=%.1f deg, Theta=%.1f deg), Polarization=%s: %.2f dBsm\n', ...
                    incidentAnglePhi, incidentAngleTheta, polarization, rcs_dbsm);

            catch ME
                warning('Error processing file: %s\nError Message: %s', currentFile, ME.message);
                rcsValues{i} = NaN;
            end

        catch ME
            warning('Error processing file: %s\nError Message: %s\nStack Trace:\n', ...
                currentFile, ME.message);
            disp(ME.stack); % Display the stack trace for debugging
            for k = 1:length(ME.stack)
                fprintf('In %s at line %d\n', ME.stack(k).file, ME.stack(k).line);
            end
            rcsValues{i} = NaN;
        end
    end

    % --- 4. Visualize Geometries (Optional) ---
    if numGeometries > 0
        % Dynamically calculate rows and columns for subplots
        numRows = ceil(sqrt(numGeometries));
        numCols = ceil(numGeometries / numRows);

        figure;
        for i = 1:numGeometries
            subplot(numRows, numCols, i);
            if (~isempty(geometryObjects{i}))
                show(geometryObjects{i});
                title(strrep(stlFiles{i}, '_', 'TITLE')); % Replace underscore for title
                xlabel('X (m)');
                ylabel('Y (m)');
                zlabel('Z (m)');
                axis equal;
                view([135 30]); % Adjust view angle as needed
            else
                title(['Error loading: ' strrep(stlFiles{i}, '_', 'TITLE')]);
            end
        end
    end

    % --- 5. Compare RCS Values ---
    fprintf('\n--- RCS Comparison Summary ---\n');
    for i = 1:numGeometries
        fprintf('  %s: ', stlFiles{i});
        if (~isnan(rcsValues{i}))
            fprintf('%.2f dBsm\n', rcsValues{i});
        else
            fprintf('Error in RCS calculation\n');
        end
    end
    
    disp(['--- Finished Polarization: ', polarization, ' ---']);
    disp('Script execution complete.');
    
end
