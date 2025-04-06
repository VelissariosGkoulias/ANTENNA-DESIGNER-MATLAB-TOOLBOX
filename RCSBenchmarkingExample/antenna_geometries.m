close all;
clear;
clc;

% Cube (aircraft_baseline) RCS Parameters Setup
lambda = 3.25e-2;
f1 = physconst("lightspeed")/lambda;
L = 10.16e-2;
W = 10.16e-2;
p = platform;
% Load and rescale
stl = stlread('aircraft_baseline.stl');
vertices = stl.Points / 100;  % Convert mm → dm
faces = stl.ConnectivityList;
p.Geometry = triangulation(faces, vertices); % Set geometry using triangulation
p.Units = "m"; % Set units to meters
figure(3)
show(p)

% Analyze and Compare with Analytical result
az = 0;
el = 0.05:1:90;
sigma = rcs(p,f1,az,el,Polarization="HH"); 
asigma1 = rectPlateRCS(L,W,f1,az,90-el);
figure(1)
plot(el,sigma,el,asigma1)
grid on
xlabel("Elevation angle (deg.)")
ylabel("RCS - dBsm")
title("Square Plate - Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")

% Cylinder (aircraft_modified1) RCS Parameters Setup
R = 10.16e-2;
pc = platform;
pc.FileName = "aircraft_modified1.stl";
pc.Units = "m";

% Analyze and Compare with analytical result
az = 0;
el = 0.05:1:90;
sigmaV = rcs(pc,f1,az,el,Polarization="HH"); 
asigma1 = circPlateRCS(R,f1,90-el);
figure(2)
plot(el,sigmaV,el,asigma1)
grid on
xlabel("Elevation angle (deg.)")
ylabel("RCS - dBsm")
title("Circular Plate - Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")

% Sphere (aircraft_modified2) Setup
p = platform;
p.FileName = "aircraft_modified2.stl";
p.Units = "m";
figure(3)
show(p)

% Analysis parameters
f2 = 7e9;
m=mesh(p,MaxEdgeLength=0.0035);
az = 0:1:180;
el =0;

% RCS Calculation with HH-Polarization
sigmahh_po = rcs(p,f2,az,el,Solver="PO",EnableGPU=false,Polarization="HH");     
sigmahh_mom = rcs(p,f2,az,el,Solver="MoM",Polarization="HH");

% RCS Calculation with VV-Polarization
sigmavv_po = rcs(p,f2,az,el,Solver="PO",EnableGPU=false,Polarization="VV");    
sigmavv_mom = rcs(p,f2,az,el,'Solver','MoM','Polarization','VV');
