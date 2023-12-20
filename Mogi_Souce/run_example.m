%% This script illustrates the use of MogiLayers. First, get_scaleN_mogi.m
% is run to examine the influence of number of Hankel transform terms on
% the solution. Then the MogiLayers solution with layered elastid properties 
% is compared with the uniform halfspace solution

%% run get_scaleN_mogi to test number of terms in Hankel transform
%INPUTS:
%   m = 1x4 volume source geometry (length(km); length(km); length(km); length^3 (m^3))
%        (x-coord, y-coord, depth(+), volume change)
% rmax = farthest distance of observation point from source (km)
% scaleN = a positive number (set to 1 if not included as an argument)
%          scales the number of terms in the Hankel transforms
%          i.e., if there are nominally N terms, then there will be
%          round(scaleN*N) terms. The default number of terms was chosen to give
%          accurate results for observations within 300 km of the fault. 
m = [0 0 10 10];
rmax = 100; %km
scaleN=10;
[Uprop,Uexact]=get_scaleN_mogi(m,rmax,scaleN);


%% compare MogiLayers and Mogi (homogeneous) displacements 


%coordinates of points on surface
x=linspace(-40,40,48);
y=linspace(-40,40,48);
[X,Y]=meshgrid(x,y);
xloc=[X(:)';Y(:)'];

%Inputs:
%   m = 1x4 volume source geometry (length; length; length; length^3)
%        (x-coord, y-coord, depth(+), volume change)
%   xloc = 2xs matrix of observation coordinates (length)
%   h = 1xN vector of depths to bottom of N layers
%   mu = 1x(N+1) vector of shear moduli -- (last entry is shear modulus of halfspace)
%	 lam = 1x(N+1) vector of lame constant lamda -- (last entry is lamda of halfspace)


m = [0 0 10 100];
scaleN=10;
h=1:1:15;
mu = linspace(0.1,5,15);
lam = mu;
U = MogiLayers(m,xloc,h,mu,lam,scaleN);

%homogeneous halfspace
[Uhom,D,S]=Mogi(m',[xloc;zeros(size(xloc(1,:)))],0.25,1);

figure
scale = 10;
quiver(xloc(1,:),xloc(2,:),scale*U(1,:),scale*U(2,:),0)
hold on
quiver(xloc(1,:),xloc(2,:),scale*Uhom(1,:),scale*Uhom(2,:),0)
title('horizontal displacement vectors')
legend('layered','homogeneous')

figure
subplot(121)
scatter(xloc(1,:),xloc(2,:),50,U(3,:),'fill')
colorbar
colormap(jet)
axis equal
caxis([0 0.4])
title('vertical displacement -- layered')

subplot(122)
scatter(xloc(1,:),xloc(2,:),50,Uhom(3,:),'fill')
colorbar
colormap(jet)
axis equal
caxis([0 0.4])
title('vertical displacement -- homogeneous')
