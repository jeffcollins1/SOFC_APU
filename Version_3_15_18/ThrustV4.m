
function [Wprop,ThrustN] = ThrustV4(Thrustlb,Height,Velocity)
i = 1:100
ThrustN = Thrustlb*4.44822; %Convert lb thrust to Newtons
VelocityVec = zeros(size(i));
ApropVecc = zeros(size(i)); 
ApropVect = zeros(size(i));
DpropVec = zeros(size(i)); 
WpropVec1 = zeros(size(i))'; 
WpropVec2 = zeros(size(i))';
PO3Vec = zeros(size(i))'; 
Tamb = -0.0065*Height + 14.987 + 273.1; %Ambient Temperature Kelvin
gam = 1.4;
Ru = 8.314; %Gas constant, J/molK
mmair = 0.0288; %molar mass of air in kg
R = Ru/mmair;
for i = 1:100
    
mAir(i) = ThrustN/(i); %mass flow of air as a function of thrust required and velocity
PO1c = intake(Height,Velocity); 

PO2c = PO1c;

V2c(i) = Velocity + i; 

PO3c(i) = intake(Height,V2c(i)); 
sO2(i) = refpropm('S','T',Tamb,'P',PO2c,'OXYGEN','NITROGEN',[0.21,0.79]);
TO3(i) = refpropm('T','P',PO3c(i),'S',sO2(i),'OXYGEN','NITROGEN',[0.21,0.79]);
PO4 = PO3c; 
a(i) = sqrt(gam*R*TO3(i)); 
M(i) = V2c(i)/a(i); 
Dca = density(Height)

Dc(i) = Dca/(1 + 0.5*(gam-1)*M(i))^(-1/(gam-1))

WpropVec1(i) = mAir(i)*(PO3c(i)-PO2c)/Dc(i); 

ApropVecc(i) = mAir(i)/(Dc(i)*Velocity); 

DpropVecc(i) = sqrt(4*ApropVecc(i)/(pi));

VelocityVec(i) = V2c(i);
PO3Vecc(i) = PO3c(i); 
end
% plot(DpropVecc,WpropVec1)
find(DpropVecc<5,i);
Wprop = min(WpropVec1(i));
end



