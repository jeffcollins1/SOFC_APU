function[Pmax,Pr_to,Prclimb,Prcruise,time_to,time_climb] = craft('type')
if 'type' = '747'
%Altitude bands, feet
n=100; 
b0 = zeros(1,10); %Takeoff Roll
b1 = [0:n:500]; %Departure
b2 = [500:n:1500]; 
b3 = [1500:n:4500];
b4 = [4500:n:9500];
b5 = [9500:n:19500]; 
b6 = [19500:n:29500];
b7 = [29500:n:39500];
b8 = [39500:n:44500]; 
[x,m1]= size(b1);
[x,m2]= size(b2);
[x,m3]= size(b3);
[x,m4]= size(b4);
[x,m5]= size(b5);
[x,m6]= size(b6);
[x,m7]= size(b7);
[x,m8]= size(b8);
mtot = m1 + m2 + m3 + m4 + m5 +m6 + m7 + m8; 
%Altitude bands, meters
b0_m = b0;
b1_m = 0.3048*b1; 
b2_m = 0.3048*b2; 
b3_m = 0.3048*b3;
b4_m = 0.3048*b4;
b5_m = 0.3048*b5; 
b6_m = 0.3048*b6;
b7_m = 0.3048*b7;
b8_m = 0.3048*b8; 
%Change in average Velocity between each altitude band, knots
V0 = linspace(15,158,10); 
V1 = linspace(158,161,m1);
V2 = linspace(161,181,m2);
V3 = linspace(181,235,m3);
V4 = linspace(235,265,m4);
V5 = linspace(265,312,m5); 
V6 = linspace(312,375,m6);
V7 = linspace(375,425,m7); 
V8 = linspace(425,480,m8); 
%change in average Velocity across each altitude band, m/s
V0_m = 0.5114*V0; 
V1_m = 0.5144*V1;
V2_m = 0.5144*V2;
V3_m = 0.5144*V3;
V4_m = 0.5144*V4;
V5_m = 0.5144*V5; 
V6_m = 0.5144*V6;
V7_m = 0.5144*V7; 
V8_m = 0.5144*V8; 
%Pitch angle at selected flight phases
ROC_max = 7.5; %Max rate of climb, m/s
ROC_min = 0.762; %Min rate of climb, m/s
%aoa1 = (5*pi/180)*ones(1,m1); %lift off angle for 60% of recorded flight data
%aoa2 = 5*pi/180*ones(1,m2); % departure angle for 80% of recorded flight data
aoa0 = zeros(1,10); 
aoa1 = asin(ROC_max./V1_m); %lift off angle for 60% of recorded flight data
aoa2 = asin(ROC_max./V2_m); 
aoa3 = asin(ROC_max./V3_m); 
aoa4 = asin(ROC_max./V4_m);
aoa5 = asin(ROC_max./V5_m);
aoa6 = asin(ROC_max./V6_m);
aoa7 = asin(ROC_max./V7_m); 
aoa8 = asin(ROC_max./V8_m);
aoa = [aoa1,aoa2,aoa3,aoa4,aoa5,aoa6,aoa7,aoa8];
aoa_deg = (180/pi)*aoa; 
ROC1 = V1_m.*sin(aoa1); 
ROC2 = V2_m.*sin(aoa2); 
ROC3 = V3_m.*sin(aoa3);
ROC4 = V4_m.*sin(aoa4); 
ROC5 = V5_m.*sin(aoa5);
ROC6 = V6_m.*sin(aoa6); 
ROC7 = V7_m.*sin(aoa7); 
ROC8 = V8_m.*sin(aoa8); 
ROC = [ROC1,ROC2,ROC3,ROC4,ROC5,ROC6,ROC7,ROC8];
%Time spent in each altitude band
t1 = (max(b1_m) - min(b1_m))/mean(ROC1);
t2 = (max(b2_m) - min(b2_m))/mean(ROC2);
t3 = (max(b3_m) - min(b3_m))/mean(ROC3);
t4 = (max(b4_m) - min(b4_m))/mean(ROC4);
t5 = (max(b5_m) - min(b5_m))/mean(ROC5);
t6 = (max(b6_m) - min(b6_m))/mean(ROC6);
t7 = (max(b7_m) - min(b7_m))/mean(ROC7);
t8 = (max(b8_m) - min(b8_m))/mean(ROC8);
%Average acceleration in each altitude band
a0 = 0.4; %Average runway takeoff acceleration 
a1 = (max(V1) - min(V1))/t1; 
a2 = (max(V2) - min(V2))/t2; 
a3 = (max(V3) - min(V3))/t3; 
a4 = (max(V4) - min(V4))/t4; 
a5 = (max(V5) - min(V5))/t5; 
a6 = (max(V6) - min(V6))/t6; 
a7 = (max(V7) - min(V7))/t7; 
a8 = (a7);
%Total trajectory and altitude vectors
a = [a0,a1,a2,a3,a4,a5,a6,a7,a8];
t = [t1,t2,t3,t4,t5,t6,t7];
velocity = [V1_m,V2_m,V3_m,V4_m,V5_m,V6_m,V7_m,V8_m];
velocityfps = (1/0.3048)*velocity; 
altitude = [b1_m,b2_m,b3_m,b4_m,b5_m,b6_m,b7_m,b8_m];
%Atmospheric Data From Table A-16,Thermodynamics Cengel and Boles, 6e
alt1 = 0:200:7000;
alt = [alt1,8000,9000,10000,12000,14000];%
TC = [15,13.7,12.4,11.1,9.8,8.5,7.2,5.9,4.6,3.3,2.0,0.7,-0.59,-1.89,-3.19,-4.49,-5.79,-7.09,-8.39,-9.69,-10.98,-12.3,-13.6,-14.9,-16.2,-17.5,-18.8,-20.1,-21.4,-22.7,-24.0,-25.3,-26.6,-27.9,-29.2,-30.5,-36.9,-43.4,-49.9,-56.5,-56.5]'; %Temperature, Kelvin
TK = 273.1 + TC; %Temperature, Kelvin
P = [101.33,98.95,96.61,94.32,92.08,89.88,87.72,85.6,83.53,81.49,79.5,77.55,75.63,73.76,71.92,70.12,68.36,66.63,64.94,63.28,61.66,60.07,58.52,57.00,55.51,54.05,52.62,51.23,49.86,48.52,47.22,45.94,44.69,43.47,42.27,41.11,35.65,30.80,26.5,19.4,14.17]'; %Pressure, kPa
rho = [1.225,1.202,1.179,1.156,1.134,1.112,1.090,1.069,1.048,1.027,1.007,0.987,0.967,0.947,0.928,0.909,0.891,0.872,0.854,0.837,0.819,0.802,0.785,0.769,0.752,0.736,0.721,0.705,0.690,0.675,0.660,0.646,0.631,0.617,0.604,0.590,0.526,0.467,0.414,0.312,0.228]'; %Density, kg/m^3
ssound = [340.3,339.5,338.8,338,337.2,336.4,335.7,334.9,334.1,333.3,332.5,331.7,331,330.2,329.4,328.6,327.8,326.2,325.4,324.6,323.8,323,322.2,321.4,320.5,319.7,318.9,318.1,317.3,316.5,315.6,314.8,314,313.1,312.3,308.1,303.8,299.5,295.1,295.1,295.1]; %speed of sound
% %fit_rho = ; %quadratic curve fit for density as a function of altitude;
d0 = 1.22*ones(1,10);
d1 = (2.9e-9).*b1_m.^2 - 0.00011.*b1_m + 1.2;
d2 = (2.9e-9).*b2_m.^2 - 0.00011.*b2_m + 1.2;
d3 = (2.9e-9).*b3_m.^2 - 0.00011.*b3_m + 1.2;
d4 = (2.9e-9).*b4_m.^2 - 0.00011.*b4_m + 1.2;
d5 = (2.9e-9).*b5_m.^2 - 0.00011.*b5_m + 1.2;
d6 = (2.9e-9).*b6_m.^2 - 0.00011.*b6_m + 1.2;
d7 = (2.9e-9).*b7_m.^2 - 0.00011.*b7_m + 1.2;
d8 = (2.9e-9).*b8_m.^2 - 0.00011.*b8_m + 1.2;
d = [d1,d2,d3,d4,d5,d6,d7,d8]; 
%fitP = (3.5e-7).*Y^2 - 0.011.*Y + 100; kPa
P0 = 100*ones(1,10); 
P1 = (3.5e-7).*b1_m.^2 - 0.011.*b1_m + 100;
P2 = (3.5e-7).*b2_m.^2 - 0.011.*b2_m + 100;
P3 = (3.5e-7).*b3_m.^2 - 0.011.*b3_m + 100;
P4 = (3.5e-7).*b4_m.^2 - 0.011.*b4_m + 100;
P5 = (3.5e-7).*b5_m.^2 - 0.011.*b5_m + 100;
P6 = (3.5e-7).*b6_m.^2 - 0.011.*b6_m + 100;
P7 = (3.5e-7).*b7_m.^2 - 0.011.*b7_m + 100;
P8 = (3.5e-7).*b8_m.^2 - 0.011.*b8_m + 100;
P = [P0,P1,P2,P3,P4,P5,P6,P7,P8]; 
%if Y < 10000 fitT = -0.0065.*Y + 290 else fitTK = 56.5 ;
TK0 = 288*ones(1,10); %Standard atmosphere
TK1 = -0.0065.*b1_m + 290;
TK2 = -0.0065.*b2_m + 290;
TK3 = -0.0065.*b3_m + 290;
TK4 = -0.0065.*b4_m + 290;
TK5 = -0.0065.*b5_m + 290;
TK6 = -0.0065.*b6_m + 290;
TK7 = zeros(1,m7);
for i = 1:1:m7
    if b7_m(1,i) <=10000
        TK7(1,i) = -0.0065.*b7_m(1,i) + 290;
    else
        TK7(1,i) = 273.1 - 56.5; 
    end
    TK7(1,i) = TK7(1,i);
end
TK8 = 216*ones(1,m8);
T = [TK0,TK1,TK2,TK3,TK4,TK5,TK6,TK7,TK8];
%Fit Speed of Sound = (-9.4e-8)*x^2 - 0.0036*x + 340
[v,size1] = size(altitude)
ss0 = refproparray('A','T',TK0,'P',P0, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss1 = refproparray('A','T',TK1,'P',P1, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss2 = refproparray('A','T',TK2,'P',P2, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss3 = refproparray('A','T',TK3,'P',P3, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss4 = refproparray('A','T',TK4,'P',P4, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss5 = refproparray('A','T',TK5,'P',P5, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss6 = refproparray('A','T',TK6,'P',P6, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss7 = refproparray('A','T',TK7,'P',P7, 'OXYGEN','NITROGEN',[0.21,0.79]); 
ss8 = refproparray('A','T',TK8,'P',P8, 'OXYGEN','NITROGEN',[0.21,0.79]); 
M0 = V0_m./ss0;
M1 = V1_m./ss1;
M2 = V2_m./ss2;
M3 = V3_m./ss3;
M4 = V4_m./ss4;
M5 = V5_m./ss5;
M6 = V6_m./ss6;
M7 = V7_m./ss7;
M8 = V8_m./ss8;
%%Boeing 747 Data From "Simulation of a Jumbo Jet Transport Aircraft Volume
%%II: Modeling Data, NASA, 1970
SI = 5500; % Wing area, ft^2
S = SI*0.092903; %Wing area, m^2
m = 396000; %Total max takeoff weight, kg
Afront = 158; %Frontal projection of craft, m^2
W = 9.81*m; % Craft weight in Newtons 
WI = 880000; %Craft weight in pounds
e = 0.8; %Oswald wing efficiency factor
bI = 195.68; %Wingspan in feet
b = bI*0.3048; %Wingspan, meters
AR = (b^2)/S; %Aspect ratio
%%Power and Thrust relations based on Gerald Brown NASA paper, 2005
%Takeoff model, using equation 6.103 in Inroduction to Flight
Tmax = 4*60600*4.448; %Maximum takeoff thrust at sea level from 4 Rolls Royce RB211 - 524H, Newtons
PmaxHP = 4*60600*0.97*1.2; %Total horsepower required to generate max sea level thrust
takeoff_length = 3200; %Required takeoff length, meters
V_to = 82; %Takeoff speed, meters per second
mu_r = 0.02; %Rolling friction between tires and runway
cl_to = 1; % Limit lift coefficient during ground roll
h = 21.5*0.3048; %Wing height above ground,m
phi = ((16*h/b)^2)/(1 + (16*h/b)^2); %Ground effect correction factor
D_to = 0.5*1.22*((0.7*V_to)^2)*S*(0.02 + phi/(pi*e*AR)); 
L_to = 0.5*cl_to*S*1.22*(0.7*V_to)^2;
Fnet_to = 0.5*Tmax - D_to - mu_r*(W - L_to); %Average net force during ground roll 
a = Fnet_to/m; %Resulting average acceleration
time_to = sqrt(2*m*takeoff_length/Fnet_to);
Pmax = PmaxHP*0.7475*ones(1,mtot); %Total power required from power plant, kW
Pr_to = 0.5*PmaxHP*0.7475; 
%Thrust and velocity with altitude after liftoff
cl = W./(0.5*1.22.*S.*velocity.^2); %Coefficient of lift for the velocity range at sea level
%clI = 880000./(0.5*0.0023769.*5500.*velocityfps.^2);
cd = 0.019 + (cl.^2)./(pi*e*AR); 
%cdI = 0.019 + (clI.^2)/(pi*e*((195^2)/5500)); 
ratio = cl./cd; 
%ratioI = clI./cdI; 
q = 0.5*1.22*velocity.^2;
Valt = velocity.*sqrt(1.22./d); 
%D = q.*S.*cd; 
Tr = m*9.81./ratio; %Thrust required at sea level to overcome drag, Newtons
Ta = 0.5*Tmax*ones(1,mtot);  %
Taalt = 4*Tmax.*d./1.22; 
Pr_to = Tr.*velocity/1000; %Power required at sea level to achieve velocity, kW
Pa = Ta.*velocity/1000; %Power available at sea level, kW
%Pr2 = (W./ratio).*sqrt(2*W./(d1.*S.*cl));
%plot(velocity,Tr)
% 
%% Thrust Required for Climbing
TrC = W.*sin(aoa) + Tr; %Thrust required at max ROC
TaC = Tmax*d./1.22;
PrC = 0.001*TrC.*velocity.*sqrt(d./1.22);
time_climb = 0.3048*44500/7; 
%plot(velocity,PrC)

%% Cruise Performance
dcruise = min(d); 
Vcruise = max(velocity); 
clcruise = W/(0.5*dcruise*S*Vcruise^2);
cdcruise = 0.01 + (clcruise^2)/(pi*e*AR); 
ratiocruise = 15; 
TrCruise = W/ratiocruise; 
PrCruise = sqrt(2*(W^3)*(cdcruise^2)/(dcruise*S*clcruise^3))/1000; %Pr(1,453)*sqrt(dcruise/1.22); %Power required at cruise, kW  

time1 = 1:1:time_to; 
[x,s1] = size(time1);
time2= linspace(time_to,time_climb,453);
Pr_toVec = Pr_to*ones(1,s1); 
[x,s2] = size(time2);
[x,s2] = size(time2); 
time3 = max(time_to) + max(time_climb);
[timetotal] = [time1,time2,time3]; 
[Ptotal] = [Pr_toVec,PrC,PrCruise]; 
plot(timetotal,Ptotal)
end
end
% 
% % TrI = 880000./ratioI; %Thrust required, pounds


