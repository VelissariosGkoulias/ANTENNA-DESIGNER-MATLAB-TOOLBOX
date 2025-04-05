% Compare built-in cube, cylinder, sphere

cube = antennaRectangle('Length', 1, 'Width', 1);
cyl = antennaCylinder('Radius', 0.5, 'Height', 1);
sph = antennaSphere('Radius', 0.5);

objs = {cube, cyl, sph};
labels = {'Cube', 'Cylinder', 'Sphere'};
az = 0:1:360;
freq = 1e9;

figure;
hold on;
for i = 1:length(objs)
    rcs = radarCrossSection(objs{i}, freq, 'Azimuth', az, 'Elevation', 0);
    plot(az, 10*log10(rcs), 'LineWidth', 2);
end
hold off;
xlabel('Azimuth (°)');
ylabel('RCS (dBsm)');
title('Built-in Shape RCS Comparison');
legend(labels);
grid on;
