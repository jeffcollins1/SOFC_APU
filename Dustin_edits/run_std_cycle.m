function param = run_std_cycle(options,mission)
%% calculate conditions for each design trade-off per kmol/s of inlet air
[m,n] = size(options.SOFC_area);
molar_flow = ones(m,n);
options.height = mission.alt(mission.design_point)*ones(m,n); %Altitude, meters
alt_tab = [0:200:7000,8000,9000,10000,12000,14000];%
atmosphere_density = [1.225,1.202,1.179,1.156,1.134,1.112,1.090,1.069,1.048,1.027,1.007,0.987,0.967,0.947,0.928,0.909,0.891,0.872,0.854,0.837,0.819,0.802,0.785,0.769,0.752,0.736,0.721,0.705,0.690,0.675,0.660,0.646,0.631,0.617,0.604,0.590,0.526,0.467,0.414,0.312,0.228]'; %Density, kg/m^3
air_den = interp1(alt_tab,atmosphere_density,options.height);
[A1,ss] = std_atmosphere(options.height,molar_flow);%Ambient conditions as a function of altitude
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
A3 = A2;
A3.T = max(options.T_fc-.5*options.dT_fc,A2.T);
[FC,A4,E1] = std_fuelcell(options,A3);
A5 = A4;
A5.T = A4.T - (A3.T-A2.T);
RC.Q = property(A3,'h','kJ') - property(A2,'h','kJ');
A5.T = find_T(A5,property(A4,'h','kJ') - RC.Q);

[A6,T1] = expander(A5,A1.P,options.T1_eff);
[HL,B1,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,[],E1);

%% Calculate nominal and mission power
P_nominal = mission.thrust(:,:,mission.design_point)*mission.mach_num(mission.design_point).*ss./options.prop_eff/1000;%nominal power in kW
mission.air_den = interp1(alt_tab,atmosphere_density,mission.alt);
[~,mission.ss] = std_atmosphere(mission.alt,1);%Ambient conditions as a function of altitude
for i = 1:1:length(mission.alt)
    mission.power(:,:,i) = mission.thrust(:,:,i).*mission.mach_num(i).*mission.ss(i)./options.prop_eff/1000;%shaft power in kW. 
end

%% scale system to meet nominal power requirements
scale = P_nominal./options.motor_eff./(FC.Power + C1.work + T1.work + B1.work);
molar_flow = scale.*molar_flow;
vol_flow = molar_flow*28.84./air_den;%Volumetric flow at the design condition
options.SOFC_area = scale.*options.SOFC_area;

%% Re-Run with scaled system parameters
[A1,~] = std_atmosphere(options.height,molar_flow);%Ambient conditions as a function of altitude
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
A3 = A2;
A3.T = max(options.T_fc-.5*options.dT_fc,A2.T);
[FC,A4,E1] = std_fuelcell(options,A3);
A5 = A4;
A5.T = A4.T - (A3.T-A2.T);
RC.Q = property(A3,'h','kJ') - property(A2,'h','kJ');
A5.T = find_T(A5,property(A4,'h','kJ') - RC.Q);

[A6,T1] = expander(A5,A1.P,options.T1_eff);
[HL,B1,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,[],E1);
weight = system_weight(options,FC,{C1;T1;B1;},[],HL);
param = NetParam(options,FC,{C1;T1;B1;},[],HL);
param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'A6',A6;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;};

%% calculate off-design power output to meet mission profile for each condition by varying permeate pressure
weight.fuel = zeros(m,n);
battery_kJ = zeros(m,n);
fuel = zeros(m*n,1);
battery = zeros(m*n,1);
P_sys_mission = zeros(m*n,length(mission.alt));
eff_mission = zeros(m*n,length(mission.alt));
param.power_mission = zeros(m,n,length(mission.alt));
param.efficiency_mission = zeros(m,n,length(mission.alt));
parallel = false;
if parallel
    parfor par_i = 1:1:m*n
        [fuel(par_i),battery(par_i),P_sys_mission(par_i,:),eff_mission(par_i,:)] = flight_profile(options,mission,vol_flow,par_i,n);
    end
else
    for i = 1:1:m*n
        i
        [fuel(i),battery(i),P_sys_mission(i,:),eff_mission(i,:)] = flight_profile(options,mission,vol_flow,i,n);
    end
end
for i = 1:1:m
    for j = 1:1:n
        battery_kJ(i,j) = battery(n*(i-1)+j);
        weight.fuel(i,j) = fuel(n*(i-1)+j);
        param.power_mission(i,j,:) = P_sys_mission(n*(i-1)+j,:);
        param.efficiency_mission(i,j,:) = eff_mission(n*(i-1)+j,:);
    end
end
weight.fuel = weight.fuel.*options.fuel_tank_mass_per_kg_fuel; %Total LH2 storage including weight of insulated container
weight.battery = battery_kJ./options.battery_specific_energy; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
weight.total = weight.sofc + weight.comp + weight.turb + weight.hx + weight.motor + weight.battery + weight.propulsor + weight.fuel; 
param.weight = weight;
param.P_den = scale*param.NetPower./(weight.sofc + weight.comp + weight.turb + weight.hx);
end%Ends function run_cycle

function [fuel,battery,P_sys_mission,eff_mission] = flight_profile(options,mission,vol_flow,par_i,n)
fuel = 0;
battery = 0;
alt_tab = [0:200:7000,8000,9000,10000,12000,14000];%
atmosphere_density = [1.225,1.202,1.179,1.156,1.134,1.112,1.090,1.069,1.048,1.027,1.007,0.987,0.967,0.947,0.928,0.909,0.891,0.872,0.854,0.837,0.819,0.802,0.785,0.769,0.752,0.736,0.721,0.705,0.690,0.675,0.660,0.646,0.631,0.617,0.604,0.590,0.526,0.467,0.414,0.312,0.228]'; %Density, kg/m^3
i = ceil(par_i/n);
j = par_i-n*(i-1);
f = fieldnames(options);
nn = length(mission.alt);
mm = 12;
for k = 1:1:length(f)
    options2.(f{k}) = ones(mm,nn)*options.(f{k})(i,j);
end
vol_flow2 = vol_flow(i,j)*(linspace(1,.1,mm)'*ones(1,nn));%reduce volume flow to 50%, then increase P_perm to reduce oxygen and power
options2.height = ones(mm,1)*mission.alt'; %Altitude, meters
air_den = interp1(alt_tab,atmosphere_density,options2.height);
molar_flow2 = vol_flow2.*air_den/28.84;%Flow rate at altitude assuming constant volumetric flow device
[A1,~] = std_atmosphere(options2.height,molar_flow2);%Ambient conditions as a function of altitude
[A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);
A3 = A2;
A3.T = max(options2.T_fc-.5*options2.dT_fc,A2.T);
[FC,A4,E1] = std_fuelcell(options2,A3);
A5 = A4;
A5.T = A4.T - (A3.T-A2.T);
RC.Q = property(A3,'h','kJ') - property(A2,'h','kJ');
A5.T = find_T(A5,property(A4,'h','kJ') - RC.Q);
[A6,T1] = expander(A5,A1.P,options2.T1_eff);

[HL,B1,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options2,FC,[],E1);
P_sys = FC.Power + C1.work + T1.work + B1.work;
P_shaft = options2.motor_eff.*P_sys;
FTE = P_sys./(FC.H2_used.*FC.hrxnmol);
P_sys_mission = zeros(1,nn);
eff_mission = zeros(1,nn);
%find permeate pressure condition that results in correct power for each flight segment
for k = 1:1:nn
    P_req = mission.power(i,j,k);%shaft power in kW.  
    if P_req>max(P_shaft(:,k))
        [P,I] = max(P_shaft(:,k));
        battery = battery + (P_req - P)*mission.duration(k)*3600;
        fuel = fuel + (FC.H2_used(I,k))*2*mission.duration(k)*3600;
        P_sys_mission(k) = P_sys(I,k);
        eff_mission(k) = FTE(I,k);
    elseif P_req<min(P_shaft(:,k))
        [h2_use,I] = min(FC.H2_used(:,k));
        fuel = fuel + P_req/P_shaft(I,k)*(FC.H2_used(I,k))*2*mission.duration(k)*3600;
        P_sys_mission(k) = P_req/min(P_sys(:,k))*P_sys(I,k);
        eff_mission(k) = FTE(I,k);
    else
        fuel = fuel + interp1(P_shaft(:,k),FC.H2_used(:,k),P_req)*2*mission.duration(k)*3600;
        P_sys_mission(k) = interp1(P_shaft(:,k),P_sys(:,k),P_req);
        eff_mission(k) = interp1(P_shaft(:,k),FTE(:,k),P_req);
    end
end
end%Ends function flight_profile