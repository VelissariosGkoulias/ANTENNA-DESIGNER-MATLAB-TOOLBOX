% Square Plate RCS Parameters Setup
lambda = 3.25e-2;
f1 = physconst("lightspeed")/lambda;
L = 10.16e-2;
W = 10.16e-2;
p = platform;
p.FileName = "square_plate.stl";
p.Units = "m";

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
title("Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")

% Circular Plate RCS Parameters Setup
R = 10.16e-2;
pc = platform;
pc.FileName = "circular_plate.stl";
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
title("Analytical vs Numerical PO")
legend("PO-Numerical","PO-Analytical","Location","best")

% NASA Almond Setup
p = platform;
p.FileName = "NASA-Almond.stl";
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

% Plot the results
figure(4)
plot(az,sigmahh_mom,az,sigmahh_po,az,sigmavv_mom,az,sigmavv_po,LineWidth=2)
ax = gca;
ax.YLim = [-70,-15]; 
title("RCS Comparison, MoM vs. PO")
xlabel("Azimuth, deg.")
ylabel("Magnitude, dBsm")
grid on
legend("HH-pol, MoM","HH-pol, PO", "VV-pol, MoM","VV-pol, PO","Location","best")

% Dielectric NASA Almond
p = platform;
p.FileName = 'NASA_Diel_Almond_24KTri_400KTets.mat'; %f = 2.58e9
p.Units='m';
p.UseFileAsMesh = 1;
figure(5)
show(p)

%% Error from this line down
% RCS Calculation for HH-polarization
f3 = 2.58e9;
sigma_almondHH = -inf*ones(1,numel(az));
h1 = figure;
plot(az,sigma_almondHH,'LineWidth',2)
xlabel('Azimuth angle (deg.)')
ylabel('RCS (dBsm)')
title('HH-pol RCS at el =0 deg.')
grid on
for i = 1:numel(az)
    sigma_almondHH(i) = rcs(p,f3,az(i),0,'Solver','FMM','Polarization','HH');
    figure(h1)
    plot(az,sigma_almondHH,'LineWidth',2)
end

% RCS Calculation with VV-Polarization
az = [0:2:130,130.5:0.5:150,152:2:180];
sigma_almondVV = -inf*ones(1,numel(az));
h2 = figure;
plot(az,sigma_almondVV,'LineWidth',2)
xlabel('Azimuth angle (deg.)')
ylabel('RCS (dBsm)')
title('VV-pol RCS at el =0 deg.')
grid on
tstart = tic
for i = 1:numel(az)
    sigma_almondVV(i) = rcs(p,f3,az(i),0,'Solver','FMM','Polarization','VV');    
    figure(h2)
    plot(az,sigma_almondVV,'LineWidth',2)
end
