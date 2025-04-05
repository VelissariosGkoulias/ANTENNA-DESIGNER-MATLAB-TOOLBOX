% MATLAB Script for RCS Calculation of STL Geometries
% Enhanced version with physical optics approximation

clear;
close all;
clc;

%% --- 1. Define Simulation Parameters ---
frequency = 1e9;               % 1 GHz
c = physconst('LightSpeed');   % Speed of light
wavelength = c/frequency;      % Wavelength
k = 2*pi/wavelength;           % Wavenumber

% Incident wave direction (azimuth, elevation)
incidentDir = [0; 0];          % Broadside incidence
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

%% --- 4. Physical Optics RCS Calculation ---
for i = 1:numGeometries
    currentFile = stlFiles{i};
    fprintf('\nProcessing: %s\n', currentFile);
    
    try
        % Load STL file
        geom = stlread(currentFile);
        vertices = geom.Points;  % No conversion needed (already in meters)
        faces = geom.ConnectivityList;
        
        % Calculate face normals and areas
        v1 = vertices(faces(:,2),:) - vertices(faces(:,1),:);
        v2 = vertices(faces(:,3),:) - vertices(faces(:,1),:);
        faceNormals = cross(v1, v2, 2);
        faceAreas = 0.5*sqrt(sum(faceNormals.^2, 2));
        faceNormals = faceNormals ./ vecnorm(faceNormals, 2, 2);  % Normalize normals
        faceCenters = (vertices(faces(:,1),:) + vertices(faces(:,2),:) + vertices(faces(:,3),:))/3;
        
        % Incident wave vector
        ki = -[cosd(incidentDir(1))*cosd(incidentDir(2)); 
               sind(incidentDir(1))*cosd(incidentDir(2));
               sind(incidentDir(2))];
        
        % Physical Optics approximation
        [rcsHH, rcsVV, rcsHV, rcsVH] = calculateFaceRCS(faceNormals, faceAreas, faceCenters, ki, k);
        rcsResults(i).HH = 10*log10(rcsHH);
        rcsResults(i).VV = 10*log10(rcsVV);
        rcsResults(i).HV = 10*log10(rcsHV);
        rcsResults(i).VH = 10*log10(rcsVH);
        
        % Visualize geometry with normals
        figure(i);
        trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
               'FaceColor', 'cyan', 'EdgeColor', 'none');
        hold on;
        quiver3(faceCenters(:,1), faceCenters(:,2), faceCenters(:,3), ...
                faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 'r');
        title(strrep(currentFile, '_', '\_'));
        axis equal; view(3);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
    catch ME
        warning('Error processing %s: %s', currentFile, ME.message);
    end
end

%% --- 5. Display and Save Results ---
fprintf('\n=== Final RCS Results ===\n');
disp(struct2table(rcsResults));

% Plot comparison
figure(length(stlFiles)+1);
polarizationColors = lines(4);
hold on;
for p = 1:length(polarizations)
    plot([rcsResults.(polarizations{p})], 'o-', ...
        'Color', polarizationColors(p,:), ...
        'LineWidth', 2, ...
        'DisplayName', polarizations{p});
end
hold off;
grid on;
title('RCS Comparison');
xlabel('Geometry'); ylabel('RCS (dBsm)');
legend('Location','best');
xticks(1:numGeometries);
xticklabels(strrep(stlFiles, 'geometries/', ''));

save('rcs_results.mat', 'rcsResults');

%% --- Helper Function: Physical Optics RCS Calculation ---
function [rcsHH, rcsVV, rcsHV, rcsVH] = calculateFaceRCS(normals, areas, centers, ki, k)
    rcsHH = 0; rcsVV = 0; rcsHV = 0; rcsVH = 0;

    % Define normalized polarization vectors
    polX = [1 0 0];                     % HH
    polY = [0 1 0];                     % VV
    pol45 = [1 1 0]/norm([1 1 0]);      % HV
    pol135 = [0 1 1]/norm([0 1 1]);     % VH

    for f = 1:size(normals,1)
        n = normals(f,:);
        A = areas(f);
        
        if dot(n, ki) > 0
            continue; % Backface (not illuminated)
        end

        rcsHH = rcsHH + A^2 * abs(dot(n, polX))^2;
        rcsVV = rcsVV + A^2 * abs(dot(n, polY))^2;
        rcsHV = rcsHV + A^2 * abs(dot(n, pol45))^2;
        rcsVH = rcsVH + A^2 * abs(dot(n, pol135))^2;
    end

    lambda = 2*pi/k;
    scale = (4*pi)/(lambda^2);
    rcsHH = scale * rcsHH;
    rcsVV = scale * rcsVV;
    rcsHV = scale * rcsHV;
    rcsVH = scale * rcsVH;
end
