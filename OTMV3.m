function [Pintake,PA1,TA1,HA1,PA2,TA2,PA5,TA5,HA5,HA2,WC1,PT1,QA2_3,mAirIn,Membrane] = OTMV3(QE1_2, molFlowO2,Height,Velocity,PO1)
PinMax = 3550; % Max pressure into OTM
PinMin = 550; % Min Pressure into OTM
x = 35;
PVec = zeros(x,1);
O2Rvec = zeros(x,1)'; 
WcompVec = zeros(x,1)';
PnetVec = zeros(x,1)';
mFlowVec = zeros(x,1)'; 
PturbVec = zeros(x,1)'; 
AirInVec = zeros(x,1)';
MembraneVec = zeros(x,1)';
QpreheatAirVec = zeros(x,1)';
Height = Height;
Velocity = Velocity; 
PA1 = intake(Height,Velocity); 
Pamb = intake(Height,0);
TA1 = -0.0065*Height + 14.987 + 273.1; %Ambient Temperature Kelvin
h1 = refpropm('H','T',TA1,'P',PA1,'OXYGEN','NITROGEN',[0.21,0.79])/1000; %Enthalpy at compressor intake
s1 = refpropm('S','T',TA1,'P',PA1,'OXYGEN','NITROGEN',[0.21,0.79]); %Entropy at compressor intake
CompEff = 0.85; 
TurbEff = 0.85;
PTurbOut = Pamb+10;
for i = 1:x
Pperm(i) = 50;
P2out(i) =  PinMin + 100*i;
TOperation = 1023; %Operating temperature of OTM/FC 'Hotbox'
hOTM(i) = refpropm('H','T',TOperation,'P',P2out(i),'OXYGEN','NITROGEN',[0.21,0.79])/1000; %Enthalpy of air at OTM kJ/kg 
O2RT(i) = 1-(1-0.21)*(Pperm(i))/(0.21*(P2out(i)-Pperm(i))); %Theoretical Percentage Recovery of O2 through OTM
O2RA(i) = .85*O2RT(i); %Actual Percentage Recovery of O2 
O2in(i) = molFlowO2*O2RA(i);
Membrane(i) = molFlowO2/(.064*.01)
mO2in(i) = O2in(i)*0.032; %Mass Flow of O2 required at OTM permeate
AirIn(i) = molFlowO2*(0.032 + 3.76*0.028); %Mass flow of air required at OTM intake
mAirOut(i) = AirIn(i) - mO2in(i); 
h2s(i) = refpropm('H','P',P2out(i),'S',s1,'OXYGEN','NITROGEN',[0.21,0.79])/1000; %Isentropic Enthalpy out of Compressor
h2a(i) = (h2s(i) - h1)/CompEff + h1; %Actual enthalpy out of compressor
TA2(i) = refpropm('T','P',P2out(i),'H',h2a(i)*1000,'OXYGEN','NITROGEN',[0.21,0.79]); %Temperature out of compressor
Qpreheatair(i) = AirIn(i)*(hOTM(i) - h2a(i)); 
Wcomp(i) = AirIn(i)*(h2a(i) - h1); 
WcompVec(i) = Wcomp(i); 
sT1(i) = refpropm('S','T',TOperation,'P',P2out(i),'OXYGEN','NITROGEN',[0.21,0.79]); %Enthalpy into Turbine
hT1s(i) = refpropm('H','P',PTurbOut,'S',sT1(i),'OXYGEN','NITROGEN',[0.21,0.79])/1000; %Isentropic enthalpy out of Turbine
hT1a(i) = hOTM(i) - (hOTM(i) - hT1s(i))*TurbEff;
TA5(i) = refpropm('T','P',PTurbOut,'h',hT1a(i),'OXYGEN','NITROGEN',[0.15,0.85])/1000
PTurb(i) = mAirOut(i)*(hOTM(i)-hT1a(i)); 
PturbVec(i) = PTurb(i); 
PnetVec(i) = (PTurb(i) - WcompVec(i)); 
PVec(i) = P2out(i);
MembraneVec(i) = Membrane(i)
O2Rvec(i) = O2RA(i); 
mFlowVec(i) = AirIn(i);
QpreheatAirVec(i) = Qpreheatair(i);
end
% a = PVec
% b = QpreheatAirVec;
% yyaxis right
% plot(a,b,'-')
% c = PnetVec;
% yyaxis left
% plot(a,c)
% title('Intake Power and Required Preheating vs. Intake Pressure')
% xlabel('Intake Pressure (kPa)'); 
% ylabel('(kW)'); 
% % yylabelleft('Net Intake Power (kW)');
%plot(PVec,QpreheatAirVec,PVec,PnetVec)
 X = find(QpreheatAirVec < QE1_2)
 if X > 0
 Y = MembraneVec(X)
 Z = find(min(Y))
 Q = QpreheatAirVec(X)
 Pintake = PnetVec(Z);
 PA2 = PVec(Z);
 WC1 = WcompVec(Z); 
 PT1 = PTurb(Z);
 QA2_3 = QpreheatAirVec(Z);
 PO1 = Pperm;
 mAirIn = mFlowVec(Z); 
 Membrane = MembraneVec(Z);
 HA1 = h1;
 HA2 = h2a;
 PA5 = PTurbOut;
 TA5 = TA5(Z);
 HA5 = hT1a(Z);
 else 
     Pintake = 0;
 PA2 = 0;
 WC1 = 0; 
 PT1 = 0;
 QA2_3 = 0;
 PO1 = Pperm;
 HA1 - h1;
 HA2 = h2a; 
 PA5 = PTurbOut;
 TA5 = 0; 
 end
% figure;plot(X,Y,X,Q)
% title('Intake Pressure vs Net Power')
% xlabel('Intake Pressure (kPa)')
% ylabel('Net Power (kW)')
%plot(PVec,QpreheatAirVec)
% title('Intake Pressure vs Required Preheating')
% xlabel('Intake Pressure (kPa)')
% ylabel('Required Preheating (kW)')


end

