%%  Test Cycle SOFC, air cathode w/fuel recirculation
%Version B sizes for max cruise efficiency, 'Design,' and then examines
%takeoff and climb performance in 'OffDesign.'
%Create options to test
%'craft' function outputs Velocity, ambient conditions, thrust, power and
%time as a function of altitude band

n1 = 10;
n2 = 10;
options.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
options.dP_fc = 10*ones(n1,n2); 
options.ASR = 0.15*ones(n1,n2); % Area specific resistance, ohm-cm^2
options.P_fc = 1000*ones(n1,n2); %Operating pressure for SOFC
options.T_fc = 1023*ones(n1,n2); %Inlet temperature for SOFC
options.dT_fc = 50; %Maximum temperature differential, Kelvin
options.asr = 0.2; % Area specific resistance, ohm-cm^2
options.T_motor = 77*ones(n1,n1); %temperture of H2 gas after cooling superconducting motors
options.spu = 0.4*ones(n1,n2); 
options.steamratio = 0.05; %Percentage of humidification at fuel inlet

options.C1_eff = 0.80; %Mechanical efficiency of compressor 1
options.airflow = 1*ones(n1,n2); %kmol/s
options.T1_eff = 0.88; %Mechanical efficiency of turbine
options.C2_eff = 0.80; %Mechanical efficiency of compressor 2
options.Blower_eff = 0.5*ones(n1,n1); %efficiency of blower
options.Blower_dP = 20*ones(n1,n1); %Pressure rise in blower in kPa
% options.T0= 298*ones(n1,n2); %Ambient Temperature at sea level Kelvin;
% options.P0 = 100*ones(n1,n2); %Ambient Pressure at sea level, kPa
n = 10;
T = 1048;%*ones(n,n);
options.motor_eff = 0.986*ones(10,10);
ASR = 0.15;%*ones(n,n);
A1.O2 = 0.21;%*ones(n,n); %Molar flow oxygen
A1.N2 = 0.79;%*ones(n,n); %Molar flow nitrogen

L = 9;
W = 9;
P = 10;
util = 0.4;%(ones(n,1)*linspace(0.1,0.5,n))';%*ones(n,n);
r = 0;%*ones(n,n); 
[Demand,Ppropmax,P_to,TrTO,options.time_to,climb.velocityclimb,time_climb,climb.P,climb.T,climb.d,M,climb.PrC,climb.altitude,Cruise.dcruise,Cruise.PCruise,Cruise.TCruise,Cruise.Vcruise,t1,t2,t3,t4,t5,t6,t7,t8,Fnet_to,Cruise.TrCruise,clcruise,cdcruise] = craft(options);
Preq = [Demand.Pr_to,Demand.PrCb1,Demand.PrCb2,Demand.PrCb3,Demand.PrCb4,Demand.PrCb5,Demand.PrCb6,Demand.PrCb7,Demand.PrCb8,Demand.PrCruise];
OffDesign.altdata(1,:) = [climb.altitude(1,1),climb.velocityclimb(1,1),climb.P(1,1),climb.T(1,1),climb.d(1,1),climb.PrC(1,1),M(1,1)];
OffDesign.altdata(2,:) = [climb.altitude(1,5),climb.velocityclimb(1,5),climb.P(1,5),climb.T(1,5),climb.d(1,5),climb.PrC(1,5),M(1,5)];
OffDesign.altdata(3,:) = [climb.altitude(1,15),climb.velocityclimb(1,15),climb.P(1,15),climb.T(1,15),climb.d(1,15),climb.PrC(1,15),M(1,15)];
OffDesign.altdata(4,:) = [climb.altitude(1,45),climb.velocityclimb(1,45),climb.P(1,45),climb.T(1,45),climb.d(1,45),climb.PrC(1,45),M(1,45)];
OffDesign.altdata(5,:) = [climb.altitude(1,95),climb.velocityclimb(1,95),climb.P(1,95),climb.T(1,95),climb.d(1,95),climb.PrC(1,95),M(1,95)];
OffDesign.altdata(6,:) = [climb.altitude(1,195),climb.velocityclimb(1,195),climb.P(1,15),climb.T(1,195),climb.d(1,195),climb.PrC(1,15),M(1,195)];
OffDesign.altdata(7,:) = [climb.altitude(1,295),climb.velocityclimb(1,295),climb.P(1,295),climb.T(1,295),climb.d(1,295),climb.PrC(1,295),M(1,295)];
OffDesign.altdata(8,:) = [climb.altitude(1,395),climb.velocityclimb(1,395),climb.P(1,395),climb.T(1,395),climb.d(1,395),climb.PrC(1,395),M(1,395)];
OffDesign.altdata(9,:) = [climb.altitude(1,450),climb.velocityclimb(1,450),climb.P(1,450),climb.T(1,450),climb.d(1,450),climb.PrC(1,450),M(1,450)];
OffDesign.altdata(10,:) =[44000,Cruise.Vcruise,Cruise.PCruise,Cruise.TCruise,Cruise.dcruise,Demand.PrCruise,0.82];
options.Pprop_to = Demand.Pr_to;
options.PrCruise = Demand.PrCruise; 
options.Tamb = Cruise.TCruise*ones(n1,n2);
options.Pamb = Cruise.PCruise*ones(n1,n2); 
options.T0= Cruise.TCruise*ones(n1,n2); %Ambient Temperature Kelvin;
options.P0 = Cruise.PCruise*ones(n1,n2); %Ambient Pressure as a function of altitude
Cruise.Mach = M(1,460); %Mach number at cruise speed
A1.T = options.T0(1,1)*(1 + 0.5*0.4*0.82^2);%*ones(n,n); %Stagnation temperature, Kelvin
A1.P = options.P0(1,1)*(1 + 0.5*0.4*0.82^2)^(1.4/0.4);%*ones(n,n); %Stagnation pressure into compressor inlet
options.PratiomaxComp = options.P_fc./A1.P;
options.PratiominComp = options.P_fc./100;
FCArray.Power = zeros(n,n);
FCArray.Qgen = zeros(n,n);
FCArray.Efficiency = zeros(n,n);
FCArray.H2out = zeros(n,n);
FCArray.H2Oout = zeros(n,n); 
FCArray.H2in = zeros(n,n);
FCArray.H2Oin = zeros(n,n);
FCArray.iDenArray = zeros(n,n);
FCArray.CellsArray = zeros(n,n);
FCArray.FCVoltage = zeros(n,n);
FCArray.Current = zeros(n,n)
FCArray.O2out = zeros(n,n); 
FCArray.N2out = zeros(n,n); 
FCArray.H2used = zeros(n,n);
FCArray.E1H2 = zeros(n,n); 
FCArray.E1H2O = zeros(n,n);
FCArray.O2util = zeros(n,n);
eff = 0.8; 
%Vary O2 utilization and SOFC stack size to genreate unique current
%densities for cruise optimization
for i = 1:10
    for j = 1:1:10
        O2util =  i*0.035;
        Cells = 100000 + 10000*j;
        [FC,E1] = FuelCell_H2(T,ASR,A1,L,W,n,Cells,P,O2util,util,r,options);
        FCArray.Power(i,j) = FC.Power;
        FCArray.O2used(i,j) = O2util*0.21*options.airflow(1,1); 
        FCArray.O2util(i,j) = O2util;
        FCArray.Qgen(i,j) = FC.Qgen;
        FCArray.Efficiency(i,j) = FC.Efficiency;
        FCArray.H2in(i,j) = FC.FuelFlow;
        FCArray.H2out(i,j) = FC.Flow.H2;
        FCArray.H2Oout(i,j) = FC.Flow.H2O;
        FCArray.iDenArray(i,j) = FC.iDen;
        FCArray.CellsArray(i,j) = Cells;
        FCArray.FCVoltage(i,j) = FC.Voltage;
        FCArray.O2out(i,j) = FC.O2out; 
        FCArray.N2out(i,j) = FC.N2out; 
        FCArray.H2used(i,j) = FC.H2used;
        FCArray.E1H2(i,j) = E1.H2; 
        FCArray.E1H2O(i,j) = E1.H2O; 
        FCArray.H2Oin(i,j) = FC.H2Oin; 
        FCArray.hrxnmol(i,j) = FC.hrxnmol; 
        FCArray.Cells(i,j) = Cells;
        FCArray.Current(i,j) = FC.Current;
        FCArray.airin = options.airflow; 
        FCArray.P0 = A1.P*ones(10,10);
        FCArray.T0 = A1.T*ones(10,10); 
        FCArray.Pamb = options.P0;
        FCArray.Tamb = options.T0; 
    end
end
       
[param,Design,HL,intake] = run_cycle_2b(options,FCArray,E1,A1);
mission = [];



%Takeoff and Climb Performance
Takeoff.air_in_max = Design.air_in;%*(1.22/Cruise.dcruise); %Max Molar flow of air at the same volume flow rate at cruise 
Takeoff.air_in_min = 0.5*Takeoff.air_in_max; %Min airflow at takeoff assuming compressor runs at 50% power for ground condition
vflowairnominal = Design.air_in*28.84./Cruise.dcruise; %Volumeflow at cruise, m^3/s
Takeoff.Cells = Design.Cells;
Takeoff.PrTO = P_to;
Takeoff.Mach = M(1,450); 
Takeoff.T0 = 298*ones(10,10);
Takeoff.P0 = 100*ones(10,10);
Profile = zeros(10,15);
Performance = zeros(10,15); 
OffDesign.Performance = zeros(10,15); 
OffDesign.CompMass = zeros(1,10);
OffDesign.TurbMass = zeros(1,10);
OffDesign.HXmass = zeros(1,10); 
Final.HTSM_mass = zeros(1,10);  
Final.TotalMass = zeros(1,10);
Final.Compmass = zeros(1,10);
Final.Turbmass=zeros(1,10);
Final.SOFCmass = zeros(1,10);
Final.HXmass=zeros(1,10); 
Final.Pden = zeros(1,10); 
Final.Penalty = zeros(1,10); 
Final.Payload = zeros(1,10);
Final.TSFCto= zeros(1,10);
Final.TSFCcr = zeros(1,10);
Final.SOFCmass = zeros(1,10); 
Final.TotalMass = zeros(1,10);
PerformanceTable = zeros(100,15);
Final.TotalFuel = zeros(1,10);
for q = 1:10
for y = 1:10
for i = 1:1:n
    for j = 1:1:10
        if y<5
        O2util =  i*0.02;
        else
            O2util = i*0.035;
        end
        Cells = Design.Cells*(0.65 + q/20); 
        air_in = (0.5 + 0.05*j)*Design.air_in*(OffDesign.altdata(y,5)/Cruise.dcruise); 
        A1.O2 = 0.21*air_in;
        A1.N2 = 0.79*air_in;
        A1.T = OffDesign.altdata(y,4)*(1+0.5*0.4*OffDesign.altdata(y,7)^2);
        A1.P = OffDesign.altdata(y,3)*(1+0.5*0.4*OffDesign.altdata(y,7)^2)^(1.4/0.4); 
        [FC,E1] = FuelCell_H2(T,ASR,A1,L,W,n,Cells,P,O2util,util,r,options);
        FCArray.Power(i,j) = FC.Power;
        FCArray.O2used(i,j) = O2util*0.21*air_in; 
        FCArray.O2util(i,j) = O2util;
        FCArray.Qgen(i,j) = FC.Qgen;
        FCArray.Efficiency(i,j) = FC.Efficiency;
        FCArray.H2in(i,j) = FC.FuelFlow;
        FCArray.H2out(i,j) = FC.Flow.H2;
        FCArray.H2Oout(i,j) = FC.Flow.H2O;
        FCArray.iDenArray(i,j) = FC.iDen;
        %FCArray.CellsArray(i,j) = Design.Cells;
        FCArray.FCVoltage(i,j) = FC.Voltage;
        FCArray.O2out(i,j) = FC.O2out; 
        FCArray.N2out(i,j) = FC.N2out; 
        FCArray.H2used(i,j) = FC.H2used;
        FCArray.E1H2(i,j) = E1.H2; 
        FCArray.E1H2O(i,j) = E1.H2O; 
        FCArray.H2Oin(i,j) = FC.H2Oin; 
        FCArray.hrxnmol(i,j) = FC.hrxnmol; 
        %FCArray.Cells(i,j) = Design.Cells;
        FCArray.Current(i,j) = FC.Current;
        FCArray.airin(i,j) = air_in;
        FCArray.Tamb(i,j) = A1.T;
        FCArray.Pamb(i,j) = A1.P;
        
    end
end
FCArray.P0 = OffDesign.altdata(y,3)*ones(10,10); 
P = Preq(1,y); 
[Weight,paramod] = run_cycle_2_odb(options,FCArray,E1,A1,Cells);
[OffDesign,Performance,Comp,Turb,HX] = condition_2odb(P,paramod,FCArray,y,OffDesign,Weight);
OffDesign.Performance(y,:) = Performance(y,:);
OffDesign.CompMass(1,y) = Comp; 
OffDesign.TurbMass(1,y) = Turb;
OffDesign.HXMass(1,y) = HX;
OffDesign.SOFCmass(1,q) = Weight.sofc(1,1); 
end
l = q*10;
m = l-9; 
PerformanceTable(m:l,1:15) = OffDesign.Performance;
OffDesign.Intake = OffDesign.CompMass + OffDesign.TurbMass;
OffDesign.IntakeMass = max(max(OffDesign.Intake));
[v,w] = find(OffDesign.Intake==OffDesign.IntakeMass);
OffDesign.compmass = OffDesign.CompMass(v,w);
OffDesign.turbmass = OffDesign.TurbMass(v,w); 
Final.Compmass(1,q) = OffDesign.compmass;
Final.Turbmass(1,q) = OffDesign.turbmass;
Final.SOFCmass(1,q) =OffDesign.SOFCmass(1,q); 
Final.HXmass(1,q) = max(max(OffDesign.HXMass));
Final.HTSM_mass(1,q) = 5273; 
Final.TotalMass(1,q) = (OffDesign.compmass + OffDesign.turbmass + OffDesign.SOFCmass(1,q) + Final.HXmass(1,q) + Final.HTSM_mass(1,q))*1.1; %Total mass based on maximum requirements for each component and 10% for hardware
Final.Pden(1,q) = P_to/Final.TotalMass(1,q); 
Final.Penalty(1,q) = Final.TotalMass(1,q) - 0.7*4*(9670/2.2); %True mass penalty after all components have been sized for max load
%OffDesign.specs = {'Net Power',selection(x,1);'System Efficiency',selection(x,2);'FC Efficiency',selection(x,11);'Current Density',selection(x,7); 'Voltage',selection(x,13);'FC Power',selection(x,14);'Turbine Work',selection(x,10);'Compressor Work',selection(x,9);'Recirculation Work',selection(x,15); 'O2 utilization',selection(x,12);'H2 used',selection(x,4);'Air in',selection(x,3);'Heat Balance',selection(x,5);};
OffDesign.specs = {'Net Power',OffDesign.Performance(1,1);'System Efficiency',OffDesign.Performance(1,2);'FC Efficiency',OffDesign.Performance(1,11);'Current Density',OffDesign.Performance(1,7); 'Voltage',OffDesign.Performance(1,13);'FC Power',OffDesign.Performance(1,14);'Turbine Work',OffDesign.Performance(1,10);'Compressor Work',OffDesign.Performance(1,9);'Recirculation Work',OffDesign.Performance(1,15); 'O2 Util',OffDesign.Performance(1,12) ;'H2 used',OffDesign.Performance(1,4);'Air in',OffDesign.Performance(1,3);'Heat Rejected',OffDesign.Performance(1,5);};
%Find total fuel burned from takeoff and climb
OffDesign.FB_to = OffDesign.Performance(1,4)*options.time_to; %Fuel burned during takeoff, kg
%Fuel burned in each altitude band
OffDesign.FB_b1 = 0.5*(OffDesign.Performance(1,4) + OffDesign.Performance(2,4))*t1; 
OffDesign.FB_b2 = 0.5*(OffDesign.Performance(2,4) + OffDesign.Performance(3,4))*t2; 
OffDesign.FB_b3 = 0.5*(OffDesign.Performance(3,4) + OffDesign.Performance(4,4))*t3; 
OffDesign.FB_b4 = 0.5*(OffDesign.Performance(4,4) + OffDesign.Performance(5,4))*t4; 
OffDesign.FB_b5 = 0.5*(OffDesign.Performance(5,4) + OffDesign.Performance(6,4))*t5; 
OffDesign.FB_b6 = 0.5*(OffDesign.Performance(6,4) + OffDesign.Performance(7,4))*t6; 
OffDesign.FB_b7 = 0.5*(OffDesign.Performance(7,4) + OffDesign.Performance(8,4))*t7;
OffDesign.FB_b8 = 0.5*(OffDesign.Performance(8,4) + OffDesign.Performance(9,4))*t8; 
OffDesign.FB_climb = OffDesign.FB_b1 + OffDesign.FB_b2 +OffDesign.FB_b3 + OffDesign.FB_b4 +OffDesign.FB_b5+OffDesign.FB_b6+OffDesign.FB_b7 +OffDesign.FB_b8;
OffDesign.FuelReserve = OffDesign.FB_to + OffDesign.FB_climb + OffDesign.Performance(10,6)*3600; %Reserve fuel sufficient to takeoff, climb and cruise for one hour; 
%Maximum Cruise distance
OffDesign.Penalty = (OffDesign.IntakeMass - Design.Comp_mass - Design.Turb_mass);
OffDesign.FuelStart = 101000 - Final.Penalty(1,q) - OffDesign.FuelReserve; 
OffDesign.FuelTotal = (101000 - Final.Penalty(1,q))/1.15; %Total Fuel Storage after accounting for insulated container
OffDesign.Cruise_fuel = OffDesign.FuelTotal - OffDesign.FB_to - OffDesign.FB_climb - OffDesign.FuelReserve; %Total fuel available after takeoff and climb
OffDesign.Cruise_time = (OffDesign.Cruise_fuel)/OffDesign.Performance(10,4); %Cruise time with useable fuel mass
OffDesign.Cruise_distance_m = OffDesign.Cruise_time*OffDesign.altdata(10,2); %Total Cruise distance, meters
OffDesign.Cruise_distance_nm = OffDesign.Cruise_distance_m/1852; %Total Cruise distance, nautical miles 
OffDesign.Cruise_nm_matchstandard = 7260; %standard range of 747
OffDesign.Cruise_matchtime = 7260*1852/OffDesign.altdata(10,2); 
OffDesign.Cruise_matchfuelmass = OffDesign.Cruise_matchtime*OffDesign.Performance(10,4); %Fuel required to match standard 747 cruise range
OffDesign.PayloadatRange = (101000 + 112760) - OffDesign.Cruise_matchfuelmass - Final.Penalty; %Available payload after mass penalties 
OffDesign.Rangeatpayload = (OffDesign.Cruise_fuel/OffDesign.Performance(10,4))*OffDesign.altdata(10,2); 
%OffDesign.tsfc_to =OffDesign.H2usedactual*1000000/TrTO; %g/kN*s
%OffDesign.tsfc_cruise = Design.H2usedactual*1000000/param.Thr_cruise; %g/kN*s
W0 = 396000*9.81; %Gross weight of craft
W1 = W0 - 9.81*OffDesign.Cruise_fuel; %Empty weight of craft
TSFCcruise = OffDesign.Performance(10,4)*1000000/Cruise.TrCruise; %g/kNs
TSFCtakeoff = OffDesign.Performance(1,4)*1000000/TrTO; %g/kNs
ct = 9.81*TSFCcruise/1000000; %N/N*s
Range = 2*sqrt(2/(Cruise.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))*(sqrt(W0) - sqrt(W1)); %Total Range according to eqn 6.77, introduction to flight, at standard payload
Rangenm = Range/1852;
W1b = (sqrt(W0) - (7260*1852)/(2*sqrt(2/(Cruise.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))))^2; %Total fuel consumption to match standard 747 range
matchfuelmass = (W0 - W1b)/9.81;
payloadboost = OffDesign.Cruise_fuel - matchfuelmass;
payloadatrange = 112760 + payloadboost; 
Final.Payload(1,q) = payloadatrange; 
Final.TSFCto(1,q) = TSFCtakeoff;
Final.TSFCcr(1,q) = TSFCcruise; 
Final.TotalFuel(1,q) = OffDesign.FuelTotal; 
end

Final.Table = [Final.Payload', Final.Turbmass',Final.Compmass',Final.HXmass',Final.SOFCmass', Final.HTSM_mass',Final.Penalty',Final.TotalMass',Final.TSFCto',Final.TSFCcr',Final.TotalFuel'];