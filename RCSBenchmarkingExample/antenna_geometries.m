%% RCS Benchmarking Example
% Follow Mathworks ecxample and adjust it accordingly
% testtesttest

close all;
clear;
clc;

lambda = 3.25e-2;
f1 = physconst("lightspeed")/lambda;


%% Cube (aircraft_baseline) RCS Parameters Setup
p = platform;
p.FileName = "aircraft_baseline.stl";
p.Units = "mm";
figure(1)
show(p)
title("Aircraft Baseline Model")

% Analyze and Compare with Analytical result
az = 0;
el = 0.05:1:90;
sigma = rcs(p,f1,az,el,Polarization="HH"); 
asigma1 = rectPlateRCS(1,1,f1,az,90-el);
figure(2)
plot(el,sigma,el,asigma1)
grid on
xlabel("Elevation angle (deg.)")
ylabel("RCS - dBsm")
title("Cube (aircraft_baseline) - Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")


%% Cylinder (aircraft_modified1) RCS Parameters Setup
R = 1;
pc = platform;
pc.FileName = "aircraft_modified1.stl";
pc.Units = "cm";
figure(3)
show(pc)
title("Aircraft Modified 1 Model")

% Analyze and Compare with analytical result
az = 0;
el = 0.05:1:90;
sigmaV = rcs(pc,f1,az,el,Polarization="HH"); 
asigma1 = circPlateRCS(R,f1,90-el);
figure(4)
plot(el,sigmaV,el,asigma1)
grid on
xlabel("Elevation angle (deg.)")
ylabel("RCS - dBsm")
title("Cylinder (aircraft_modified1) - Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")


%% Sphere (aircraft_modified2) Setup
p = platform;
p.FileName = "aircraft_modified2.stl";
p.Units = "cm";
figure(5)
show(p)

% Analysis parameters
f2=7e9;
m=mesh(p,MaxEdgeLength=0.35);
az=0:1:180;
el=0;

%% Check mesh number off triangles for sphere (128 at the moment with 0.35 MaxEdgeLength)
% Assuming 'm' is the mesh object
triang = m.Triangulation;  % Get the triangulation object

% Now you can access the faces and vertices
numTriangles = size(triang, 1);  % Number of triangles corresponds to rows in the triangulation

% Display the result
disp(['Number of triangles: ', num2str(numTriangles)]);


% RCS Calculation with HH-Polarization
sigmahh_po = rcs(p,f2,az,el,Solver="PO",EnableGPU=false,Polarization="HH");     
sigmahh_mom = rcs(p,f2,az,el,Solver="MoM",Polarization="HH");

% RCS Calculation with VV-Polarization
sigmavv_po = rcs(p,f2,az,el,Solver="PO",EnableGPU=false,Polarization="VV");    
sigmavv_mom = rcs(p,f2,az,el,'Solver','MoM','Polarization','VV');

% Plot the results
figure(6)
plot(az,sigmahh_mom,az,sigmahh_po,az,sigmavv_mom,az,sigmavv_po,LineWidth=2)
ax = gca;
ax.YLim = [-70,-15]; 
title("RCS Comparison, MoM vs. PO")
xlabel("Azimuth, deg.")
ylabel("Magnitude, dBsm")
grid on
legend("HH-pol, MoM","HH-pol, PO", "VV-pol, MoM","VV-pol, PO","Location","best")
