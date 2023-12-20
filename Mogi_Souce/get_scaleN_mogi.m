function [Uprop,Uexact]=get_scaleN_mogi(m,rmax,scaleN)

%[Uprop,Uexact]=get_scaleN_mogi(m,rmax,scaleN)
%
%Program to determine the number of terms needed in the Hankel transforms
%in MogiLayers.m. Calculates the displacements in a homogeneous
%halfspace (and poisson ratio = 0.25) using MogiLayers.m and compares with 
%the exact analytical
%solution using Mogi.m. Displacments are calculated at 50 observation points 
%spaced evenly along a linestretching from 0 to rmax kms from the mogi source. 
%Vary scaleN until the propagator matrix solution agrees with the exact Mogi
%solution within the desired accuracy.

%The number of terms in the transforms is round(N*scaleN) where N is the default 
%number of terms. If too few terms are used there may be aliasing of the solution,
%that is, short wavelength deformation is wrapped around into the long wavelengths, 
%creating errors in the solution at large distances from the fault. The 
%default value (scaleN=1) gives accurate results for observations within 
%100 km of the source. If your farthest observation point is greater than 100 km 
%from the source, you may need more than the default number (scaleN>1). If your
%farthest observation point is less than 100 km, you may be able to get away 
%with fewer than the default number of terms (scaleN<1). There is always a 
%trade-off between computation time and accuracy (number of terms), so determine 
%the fewest number of terms needed to obtain the desired accuracy if you care 
%about computation time. MogiLayers.m will run faster if you compile it to C 
%using mcc.
%
%%INPUTS:
%   m = 1x4 volume source geometry (length; length; length; length^3)
%        (x-coord, y-coord, depth(+), volume change)
% rmax = farthest distance of observation point from source (km)
% scaleN = a positive number (set to 1 if not included as an argument)
%          scales the number of terms in the Hankel transforms
%          i.e., if there are nominally N terms, then there will be
%          round(scaleN*N) terms. The default number of terms was chosen to give
%          accurate results for observations within 300 km of the fault. 
%%OUTPUTS:
%Uprop is a 3x50 matrix of displacements using propagator matrix solution
%      (first row east, second row north, third row up) 
%Uexact is a 3x50 matrix of displacements using exact Mogi solution
%      (first row east, second row north, third row up) 
%Kaj M. Johnson, Stanford and UC Berkeley, last updated 1/10/05

if nargin==2
    scaleN=1;
end

r=linspace(0,rmax,50);
x=r*cos(45*pi/180);
y=r*sin(45*pi/180);
xloc=[x;y];

xloc3d=[xloc;zeros(size(x))];

%parameters for homogeneous halfspace
h=[1];
mu=[1 1];
lam=[1 1];

Uprop=MogiLayers(m,xloc,h,mu,lam,scaleN);
[Uexact,D,S]=Mogi(m',xloc3d,.25,1);

figure
subplot(311); p=plot(r,Uprop(1,:)); hold on; plot(r,Uprop(1,:),'.'); o=plot(r,Uexact(1,:),'r'); plot(r,Uexact(1,:),'r.'); legend([p,o],'propagator','exact'); title('displacements in x direction')
subplot(312); p=plot(r,Uprop(2,:)); hold on; plot(r,Uprop(2,:),'.'); o=plot(r,Uexact(2,:),'r'); plot(r,Uexact(2,:),'r.'); legend([p,o],'propagator','exact'); title('displacements in y direction')
subplot(313); p=plot(r,Uprop(3,:)); hold on; plot(r,Uprop(3,:),'.'); o=plot(r,Uexact(3,:),'r'); plot(r,Uexact(3,:),'r.'); legend([p,o],'propagator','exact'); title('displacements in z direction')
xlabel('distance from bottom edge (km)')

figure
subplot(311); p=plot(r,(Uprop(1,:)-Uexact(1,:))); title('difference -- displacements in x direction')
subplot(312); p=plot(r,(Uprop(2,:)-Uexact(2,:))); title('difference -- displacements in y direction')
subplot(313); p=plot(r,(Uprop(3,:)-Uexact(3,:))); title('difference -- displacements in z direction')
xlabel('distance from bottom edge (km)')


figure
subplot(311); p=plot(r,100*(Uprop(1,:)-Uexact(1,:))./Uexact(1,:)); title('%difference -- displacements in x direction')
subplot(312); p=plot(r,100*(Uprop(2,:)-Uexact(2,:))./Uexact(2,:)); title('%difference -- displacements in y direction')
subplot(313); p=plot(r,100*(Uprop(3,:)-Uexact(3,:))./Uexact(3,:)); title('%difference -- displacements in z direction')
xlabel('distance from bottom edge (km)')

