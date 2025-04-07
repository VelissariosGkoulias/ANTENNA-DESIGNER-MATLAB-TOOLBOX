%% MATLAB Script for RCS Calculation of STL Geometries using Antenna Toolbox

clear;
close all;
clc;

%% --- 1. Define Simulation Parameters ---
frequency = 1e9;         % 1 GHz
incidentAnglePhi = 0;    % Incident angle in phi (azimuth)
incidentAngleTheta = 0;  % Incident angle in theta (elevation)
polarizations = {'HH', 'VV', 'HV', 'VH'};

%% --- 2. Specify STL File Paths ---
stlFiles = {
    'geometries/aircraft_baseline.stl';
    'geometries/aircraft_modified1.stl';
    'geometries/aircraft_modified2.stl'
};

%% --- 3. Initialize Results Storage ---
numGeometries = length(stlFiles);
rcsResults = struct('FileName', stlFiles, ...
                    'HH', NaN(numGeometries,1), ...
                    'VV', NaN(numGeometries,1), ...
                    'HV', NaN(numGeometries,1), ...
                    'VH', NaN(numGeometries,1));

geometryObjects = cell(numGeometries, 1); % To store the custom antenna objects

%% --- 4. Antenna Toolbox RCS Calculation ---
for i = 1:numGeometries
    currentFile = stlFiles{i};
    fprintf('\nProcessing: %s\n', currentFile);

    try
        % Load STL file
        geom = stlread(currentFile);
        vertices = geom.Points;
        faces = geom.ConnectivityList;

        % Validate the STL file
        if isempty(faces) || isempty(vertices)
            error('Invalid STL file: %s. Faces or vertices are empty.', currentFile);
        end

        % Create a triangulation object
        tri = triangulation(faces, vertices);

        % Create a custom antenna geometry object using the triangulation
        try
            customAntenna = customAntennaGeometry(tri);
        catch ME
            warning('Error creating custom antenna geometry for %s: %s', currentFile, ME.message);
            continue;
        end

        % Refine the mesh (optional, adjust as needed)
        customAntenna.Mesh.MaxEdgeLength = 0.1; % Set MaxEdgeLength property
        customAntenna = mesh(customAntenna);     % Call mesh without parameters

        % Store the geometry object for potential visualization
        geometryObjects{i} = customAntenna;

        % Create a radar problem object using the Method of Moments (MoM) solver
        radar = radarProblem('Solver', 'MoM', 'Frequency', frequency);

        % Add the custom antenna as a target
        target = customAntenna;
        add(radar, target);

        % Calculate the monostatic RCS for each polarization
        for p_idx = 1:length(polarizations)
            current_polarization = polarizations{p_idx};
            try
                sigma = rcs(radar, frequency, ...
                          'IncidentAngles', [incidentAnglePhi, incidentAngleTheta], ...
                          'Polarization', current_polarization);

                % Store the RCS value (in dBsm)
                rcs_dbsm = 10*log10(abs(sigma));
                rcsResults(i).(current_polarization) = rcs_dbsm;
            catch ME
                warning('Error calculating RCS for %s (Polarization: %s): %s', ...
                        currentFile, current_polarization, ME.message);
                rcsResults(i).(current_polarization) = NaN;
            end
        end

        % Visualize geometry (optional)
        if ~isempty(customAntenna.Faces) && ~isempty(customAntenna.Vertices)
            figure(i);
            show(customAntenna);
            title(strrep(currentFile, '_', '\_'));
            axis equal; view(3);
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        else
            warning('Invalid geometry for visualization: %s', currentFile);
        end

    catch ME
        warning('Error processing %s: %s', currentFile, ME.message);
        for k = 1:length(ME.stack)
            fprintf('In %s at line %d\n', ME.stack(k).file, ME.stack(k).line);
        end
    end
end

%% --- 5. Display and Save Results ---
fprintf('\n=== Final RCS Results (Antenna Toolbox) ===\n');
disp(struct2table(rcsResults));

% Plot comparison
figure(length(stlFiles)+1);
polarizationColors = lines(4);
hold on;
for p_idx = 1:length(polarizations)
    plot([rcsResults.(polarizations{p_idx})], 'o-', ...
         'Color', polarizationColors(p_idx,:), ...
         'LineWidth', 2, ...
         'DisplayName', polarizations{p_idx});
end
hold off;
grid on;
title('RCS Comparison (Antenna Toolbox)');
xlabel('Geometry'); ylabel('RCS (dBsm)');
legend('Location','best');
xticks(1:numGeometries);
xticklabels(strrep(stlFiles, 'geometries/', ''));

save('rcs_results_at.mat', 'rcsResults');

disp('Script execution complete.');
