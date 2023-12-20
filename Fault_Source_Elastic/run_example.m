
%% This script illustrates the use of LayeredGreens. First, get_scaleN.m
% is run to examine the influence of number of Hankel transform terms on
% the solution. Then the LayeredGreens solution with layered elastic properties 
% is compared with the uniform halfspace solution


%% Run get_scaleN to test the number of terms need in the Hankel transform

%%INPUTS:
%fault model:
% m(1) = length (km) (do not give distances in meters)
% m(2) = width (km)
% m(3) = depth to down-dip edge (km)
% m(4) = dip (degrees)
% rmax = farthest distance of observation point from fault patch (km)
% scaleN = a positive number (set to 1 if not included as an argument)
%          scales the number of terms in the Hankel transforms
%          i.e., if there are nominally N terms, then there will be
%          round(scaleN*N) terms. The default number of terms was chosen to give
%          accurate results for observations within 300 km of the fault. 

m = [10 10 20 40];
rmax = 500;
scaleN=10;
[Uprop,Uokada] = get_scaleN(m,rmax,scaleN);


%% compare displacements computd with LayeredGreens and the solution for a uniform 
% elastic halfspace

%coordinates of points on surface
x=linspace(-40,40,24);
y=linspace(-40,40,24);
[X,Y]=meshgrid(x,y);
xloc=[X(:)';Y(:)'];


%fault geometry  -- type 'help LayeredGreens' for more info
%m=[length,width,depth of bottom edge,dip,strike,east*,north*,strike-slip,dip-slip,opening]
% *position of midpoint of bottom edge
m=[10 10 20 30 0 0 0 0 5 0];

% d = 1x(M-1) vector of depths to bottom of layer -- M is number of layers
%     (including halfspace); note: bottom layer is halfspace, thus no depth is specified
% mu,lam = 1xM vector of normalized Lame constants (normalized with value of mu in top layer)
%          mu=lam for poisson ratio = 0.25, M is number of layers (including halfspace)
%          It is IMPORTANT that the Lame constants used are relative to (normalized by)
%          the value in the top layer
d = 2:2:20;
mu = linspace(.1,5,10);
lam = mu;

U=LayeredGreens(m,xloc,d,mu,lam);

[Uhom,D,S]= disloc3d(m,[xloc;zeros(size(xloc(1,:)))],1,0.25);

figure
scale = 50;
quiver(xloc(1,:),xloc(2,:),scale*U(1,:),scale*U(2,:),0)
hold on
quiver(xloc(1,:),xloc(2,:),scale*Uhom(1,:),scale*Uhom(2,:),0)
legend('layered','homogeneous')
axis equal


figure
subplot(121)
scatter(xloc(1,:),xloc(2,:),50,U(3,:),'fill')
colorbar
colormap(jet)
axis equal
caxis([-0.2 0.2])
title('vertical displacement -- layered')

subplot(122)
scatter(xloc(1,:),xloc(2,:),50,Uhom(3,:),'fill')
colorbar
colormap(jet)
axis equal
caxis([-0.2 0.2])
title('vertical displacement -- homogeneous')

