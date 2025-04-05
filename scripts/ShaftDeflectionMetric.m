% Stolen from Dr. S

function [x, y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionMetric(F,Floc,d,dloc,Rloc,L)
%This program finds deflection y and slope, theta for a steel stepped shaft
%with point loads.  Program handles loads/deflections in 2-D. 
%Units in lb and in.
%  Required user inputs are:
%   F - list of applied forces in one plane
%   Floc - location of the forces from left end
%   d - vector of shaft diameters
%   dloc - location of beginning of shaft diameter
%   Rloc - location of bearings (limit 2)
%   L - total length
%  Function outputs are x:
%   x - vector of locations
%   y - vector of deflection (same direction as F) to match x
%   dydx - vector of slopes in radians to match x
%   M & MdEI - vector of bending moment M, or M/EI, useful for plots
%   R vector of reaction forces at the bearings
%   diam %the diameters of the shaft as a function of x
%   EI %product of modulus and moment of inertia as a function of x

dx=.0001*L; %in  Increment for numerical integration

E=200e9; %psi  Modulus of Elascticity for steel

%%Use statics to solve for reaction forces
%Sum moments about first support
R(2)=-F*(Floc-Rloc(1))'/(Rloc(2)-Rloc(1));
%Sum forces in the y direction
R(1)=-(sum(F)+R(2));

%%Generate a vector of points along the length.
x=0:dx:L;  
x=[x Floc, dloc, Rloc, L];  %Make sure key points are included.
%Omit duplicate points but have a copy of pts where the diameter changes.
x=[unique(x) dloc(2:end)]; 
x=sort(x); %Put vector in ascending order.

%%Generate V(x), M(x),M(x)/EI(x)
V=zeros(1,length(x));
M=zeros(1,length(x));

diam=ones(1,length(x));
dloc=[dloc L];
for i=1:length(dloc)-1
    j=find(x==dloc(i),1,'last');
    k=find(x==dloc(i+1),1,'first');
    diam(j:k)=d(i);
end

for i=1:length(x)
   for j=1:2  %Account for reaction forces
        if x(i)>=Rloc(j)
            V(i)=V(i)+R(j);
            M(i)=M(i)+R(j)*(x(i)-Rloc(j));
        end
   end
   for j=1:length(F)  %Account for applied forces
       if x(i)>=Floc(j)
            V(i)=V(i)+F(j);
            M(i)=M(i)+F(j)*(x(i)-Floc(j));
       end
   end
end

EI=E*(pi/4)*(diam/2).^4;
MdEI=M./EI;

y(1)=0;
dydx(1)=0;
%Integrate twice to get slope and deflection, excluding constants.
dydx=cumtrapz(x,MdEI);
y=cumtrapz(x,dydx);

%%Get the constants of integration by solving BCs
%BC1:  0=y(Rloc(1))+c0*Rloc(1)+c1  , zero deflection at bearing 1
%BC2:  0=y(Rloc(2))+c0*Rloc(2)+c1 ,  zero deflection at bearing 2
%Matrix Form:  [y(Rloc(1));y(Rloc(2))]=[Rloc(1) 1;Rloc(2) 1][c0; c1]
for j=1:2  %Get uncorrected y at supports, put in vector b
    Rindex(j)=find(x>=Rloc(j),1);  
    b(j)=y(Rindex(j));
end
C=[Rloc;1 1]';
constants=C\-b';  %Solve matrix form of simultaneous equations.

%Adjust dydx and y based on these constants
dydx=dydx+constants(1);
y=y+constants(1)*x+constants(2);

end

