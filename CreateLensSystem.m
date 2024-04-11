%% Clean up
clearvars –global
close all 
clc
clear all 
clear classes 

%% Input Values 
% Effective focal length of the system [mm]
EFL = 50; 

% K - Multiplier of EFL - positioning of the object to the left (must be 
% negative) of the first vertex of the first surface (Object is K*EFL [mm]
% to the left from fist vertex of furst lens)
K = -2;  

% Half of the effective aperture
H1=10;

% Refractive indices of chosen glasses and air for light lengths 
% (d) - basic, (e) (F') (C') - for which chromatic aberration is corrected,
% Lens system is in the air.

n1d=1.00027717; % Air refractive index for light length (d) = 0.5876 µm
n1ff=1.00027954; % Air refractive index for light length (F') = 0.4800 µm
n1cc=1.00027640; % Air refractive index for light length (C') = 0.6438 µm
n1e=1.00027791; % Air refractive index for light length (e) = 0.5461 µm

% NBAF10 - SCHOTT Catalogue - First lens material 
n2d=1.67; % NBAF10  refractive index for light length (d) = 0.5876 µm
n2ff=1.6808; % NBAF10 refractive index for light length (F') = 0.4800 µm
n2cc=1.6665; % NBAF10 refractive index for light length (C') = 0.6438 µm
n2e=1.6734; % NBAF10 refractive index for light length (e) = 0.5461 µm

% NSF10 - SCHOTT Catalogue - Second lens material
n3d=1.7283;% NSF10 refractive index for light length (d) = 0.5876 µm
n3ff=1.7480; % NSF10 refractive index for light length (F') = 0.4800 µm
n3cc=1.7221; % NSF10 refractive index for light length (C') = 0.6438 µm
n3e=1.7343; % NSF10 refractive index for light length (e) = 0.5461 µm

% Paraxial height - chosen as 1
h1=1;

% Maximum length of the system from first to last(forth) vertex 
Lmax=25; 

%% Creating Lens System object
% Lens system object mostly has symbolic equations depending on curvature 
% radii of surfaces and thickness of lenses resp. distance inbetween the 
% lenses f(Ri, ti) as its parameters 
lensSystem = LensSystemDesign(EFL,K,H1,n1d, ...
                n2d,n3d,n1ff,n1cc, ...
                n1e,n2ff,n2cc,n2e, ...
                n3ff,n3cc,n3e,h1, ...
                Lmax);
%% Optimization with GA algorithm

% Limiting conditions 
% x is vector of parameters in alphabetical order. Here example
% [curRad1 curRad2 curRad3 curRad4 prevDist2 prevDist3 prevDist4]
% curRadX - curvature radius of optical surface X
% prevDistX - previous optic surface vertex distance from optic surface X

% A*x<=b. Limits chosen by user.
% Example: A=[0 0 0 0 1 1 1], b=[lensSystem.Lmax] -> 
% -> prevDist2 + prevDist3 + prevDist4 <= Lmax
A=[0 0 0 0 1 1 1; 0 0 0 0 1 0 0; 0 0 0 0 -1 0 0;0 0 0 0 0 0 1; ...
    0 0 0 0 0 0 -1];
b=[lensSystem.Lmax 2*lensSystem.H1/7 -2*lensSystem.H1/11 ...
    2*lensSystem.H1/7 -2*lensSystem.H1/11];

% Lower boundaries and upper boundaries 
% Example: lb=[15 -450 15 15 2 2 5]; ub=[450 -15 450 450 8 8 20]; -> 
% -> //15 mm <= curRad1 <= 450 mm// (Vertex is to the left, concave)
% -> //-450 mm <= curRad2 <= -15 mm// (Vertex is to the right, convex) ... etc.
% See standard optical convention of positive and negative curvature
lb=[20 -Inf 20 -Inf 1.5 3 1.5];
ub=[Inf -20 Inf -20 3.5 20 3.5];

% Options of genetic algorithm function - see Matlab documentation

% Constraint tolerance - lower value = longer calculation
gaOptions.ConstraintTolerance = 1e-5;

% Population Size - starting pull of points, 
% higher value = longer calculation
gaOptions.PopulationSize = 500;

% Vector with found solutions for vector of parameters in alphabetical 
% order. Here example 
% [curRad1 curRad2 curRad3 curRad4 prevDist2 prevDist3 prevDist4]
systemValues = lensSystem.GAfunction(A, b, lb, ub, gaOptions);

%% Calculation of system parameters
% Calculation of parameters of the system, depending on found parameters
% curvature radii and thickness of lenses, and distance inbetween them
calculatedLensSystem = lensSystem.inserValues(systemValues);

%% Visualisation
%Graphical scheme of found lens system
calculatedLensSystem.visualiseSystem();
Answer = [calculatedLensSystem.optSurf1.objectPosition calculatedLensSystem.optSurf1.curveRadius ...
    calculatedLensSystem.optSurf2.curveRadius     calculatedLensSystem.optSurf3.curveRadius ...
    calculatedLensSystem.optSurf4.curveRadius     calculatedLensSystem.optSurf2.previousSurfaceDistance ...
    calculatedLensSystem.optSurf3.previousSurfaceDistance    calculatedLensSystem.optSurf4.previousSurfaceDistance ...
    calculatedLensSystem.optSurf4.imagePosition    calculatedLensSystem.systemCalculatedEFL ...
    calculatedLensSystem.magnification     calculatedLensSystem.systemSphericalAberration ...
    calculatedLensSystem.systemChromaticAberration   lb(1)...
   lb(2)   lb(3)   lb(4)   lb(5)   lb(6)   lb(7) gaOptions.ConstraintTolerance   gaOptions.PopulationSize;
   0 0 0 0 0 0 0 0 0 0 0 0 0 ub(1) ub(2) ub(3) ub(4) ub(5) ub(6) ub(7) 0 0];
