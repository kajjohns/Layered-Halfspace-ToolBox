function [Uprop,Uokada]=get_scaleN(m,rmax,scaleN)

%[Uprop,Uokada]=get_scaleN(m,rmax,scaleN)
%
%Program to determine the number of terms needed in the Hankel transforms
%in LayeredGreens.m. Calculates the displacements in a homogeneous
%halfspace (and poisson ratio = 0.25) using LayeredGreens.m and 
%compares with the exact analytical
%solution using Okada's dislocation code, disloc3d. Displacments are 
%calculated at 50 observation points spaced evenly along a line at 45
%degrees to the strike of a fault and stretching from -rmax to rmax kms from 
%the bottom edge of the fault. The fault has 1 unit of strike-slip
%component and 1 unit of dip-slip compenent of slip (displacements 
%are computed in the same units as slip). Vary scaleN until
%the propagator matrix solution agrees with the exact Okada dislocation
%solution within the desired accuracy.

%The number of terms in
%the transforms is round(N*scaleN) where N is the default number of terms. 
%If too few terms are used there may be aliasing of the solution, that is, 
%short wavelength deformation is wrapped around into the long wavelengths, 
%creating errors in the solution at large distances from the fault. The 
%default value  (scaleN=1) gives accurate results for observations within 
%300 km of the fault (the solution is inaccurate within 2 km of a fault that breaks 
%the ground surface). If your farthest observation point is greater than 300 km 
%from the fault, you may need more than the default number (scaleN>1). If your
%farthest observation point is less than 300 km, you may be able to get away 
%with fewer than the default number of terms (scaleN<1). 
%There is always a trade-off between computation time 
%and accuracy (number of terms), so determine the fewest number of terms
%needed to obtain the desired accuracy if you care about computation time. 
%LayeredGreens.m will run faster if you compile it to C using mcc.
%
%To be sure you are using enough terms (but not too many if you are
%concerned about computation time) for your particular problem, run this
%program using the length and width and depths of the actual fault patch(es) 
%in your model. If you have a fault surface with numerous slip
%patches, it may be necessary to test the slip patches at every depth. For
%example if you have a fault with 10 30deg dipping, 5km X 5km patches along strike and 
%10 5km X 5km patches down dip, run this program 10 times for 5km X 5km, 30deg dipping
%patches at each of the down-dip depths. You may find that scaleN varies
%with patch depth.
%
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
%%OUTPUTS:
%Uprop is a 3x50 matrix of displacements using propagator matrix solution
%      (first row east, second row north, third row up) 
%Uokada is a 3x50 matrix of displacements using exact okada solution
%      (first row east, second row north, third row up) 
%  %%NOTE: displacements in prop solution are inaccurate within about 2 km of a fault that
%           breaks the ground surface.
%Kaj M. Johnson, Stanford and UC Berkeley, last updated 1/10/05

if nargin==2
    scaleN=1;
end

m(5)=0; %strike
m(6)=0; %east position of bottom edge
m(7)=0; %north postion of bottom edge
m(8)=1; %strike-slip
m(9)=1; %dip-slip
m(10)=0; %tensile

r=linspace(-rmax,rmax,50);
x=r*cos(45*pi/180);
y=r*sin(45*pi/180);
xloc=[x;y];

xloc3d=[xloc;zeros(size(x))];

%parameters for homogeneous halfspace
d=[1];
mu=[1 1];
lam=[1 1];

Uprop=LayeredGreens(m,xloc,d,mu,lam,scaleN);
[Uokada,D,S,flag]=disloc3d(m',xloc3d,1,.25);

figure
subplot(311); p=plot(r,Uprop(1,:)); hold on; plot(r,Uprop(1,:),'.'); o=plot(r,Uokada(1,:),'r'); plot(r,Uokada(1,:),'r.'); legend([p,o],'propagator','okada'); title('displacements in x direction')
subplot(312); p=plot(r,Uprop(2,:)); hold on; plot(r,Uprop(2,:),'.'); o=plot(r,Uokada(2,:),'r'); plot(r,Uokada(2,:),'r.'); legend([p,o],'propagator','okada'); title('displacements in y direction')
subplot(313); p=plot(r,Uprop(3,:)); hold on; plot(r,Uprop(3,:),'.'); o=plot(r,Uokada(3,:),'r'); plot(r,Uokada(3,:),'r.'); legend([p,o],'propagator','okada'); title('displacements in z direction')
xlabel('distance from bottom edge (km)')

figure
subplot(311); p=plot(r,Uprop(1,:)-Uokada(1,:)); title('difference -- displacements in x direction')
subplot(312); p=plot(r,Uprop(2,:)-Uokada(2,:)); title('difference -- displacements in y direction')
subplot(313); p=plot(r,Uprop(3,:)-Uokada(3,:)); title('difference -- displacements in z direction')
xlabel('distance from bottom edge (km)')


figure
subplot(311); p=plot(r,100*(Uprop(1,:)-Uokada(1,:))./Uokada(1,:)); title('% difference -- displacements in x direction')
subplot(312); p=plot(r,100*(Uprop(2,:)-Uokada(2,:))./Uokada(1,:)); title('% difference -- displacements in y direction')
subplot(313); p=plot(r,100*(Uprop(3,:)-Uokada(3,:))./Uokada(1,:)); title('% difference -- displacements in z direction')
xlabel('distance from bottom edge (km)')