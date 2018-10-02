function param = run_cycle(options,mission)
%% calculate cruise conditions for each design trade-off
[m,n] = size(options.SOFC_area);
molar_flow = ones(m,n);
options.height = max(mission.alt)*ones(m,n); %Altitude, meters
alt_tab = [0:200:7000,8000,9000,10000,12000,14000];%
atmosphere_density = [1.225,1.202,1.179,1.156,1.134,1.112,1.090,1.069,1.048,1.027,1.007,0.987,0.967,0.947,0.928,0.909,0.891,0.872,0.854,0.837,0.819,0.802,0.785,0.769,0.752,0.736,0.721,0.705,0.690,0.675,0.660,0.646,0.631,0.617,0.604,0.590,0.526,0.467,0.414,0.312,0.228]'; %Density, kg/m^3
air_den = interp1(alt_tab,atmosphere_density,options.height);
[A1,ss] = std_atmosphere(options.height,molar_flow);%Ambient conditions as a function of altitude
[OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1);
%one-time adjustment of SOFC area per kmol air flow to ensure a feasible current density
i_den = 4000.*O5.O2.*96485.33./(options.SOFC_area*10000); %A/cm^2
i_den = max(.1,min(i_den,.4./options.asr));%Cant solve for ultra low or high current densities
options.SOFC_area = 4000.*O5.O2.*96485.33./i_den/1e4;

[FC,E1] = oxy_fuelcell(options,O5);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,OTM,E1);
param = NetParam(options,FC,OTM,HL);
param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};

%% scale system to meet cruise power requirements
T = options.TO_weight./options.Lift_2_Drag*4.448;% Thrust in N:  Thrust = Cd/Cl*Weight; cruise Cd = 0.036; Cl = .48; weight = 820,000lbs; thus 61,000lbs cruise thrust. Makes sense because max thrust at cruise is 100,000lbs.Range = 6100nautical miles, velocity = 635knts
P = T.*max(mission.mach_num).*ss./options.prop_eff/1000;%shaft power in kW.  propellor momentum theory at cruise http://164.100.133.129:81/econtent/Uploads/09-%20Ducted%20Fans%20and%20Propellers%20%5BCompatibility%20Mode%5D.pdf ,  http://web.mit.edu/16.unified/www/SPRING/systems/Lab_Notes/airpower.pdf
scale = P./options.motor_eff./(FC.Power + OTM.net_work + HL.blower_work);
molar_flow = scale.*molar_flow;
vol_flow = molar_flow*28.84./air_den;%Volumetric flow at the design condition
[weight,options] = system_weight(options,scale,FC,OTM,HL,A1);

%%calculate off-design power output to meet mission profile for each condition by varying permeate pressure

weight.fuel = zeros(m,n);
battery_kJ = zeros(m,n);
fuel = zeros(m*n,1);
battery = zeros(m*n,1);
parallel = true;
if parallel
    parfor par_i = 1:1:m*n
        [fuel(par_i),battery(par_i)] = flight_profile(options,mission,vol_flow,par_i,n);
    end
else
    for i = 1:1:m*n
        [fuel(i),battery(i)] = flight_profile(options,mission,vol_flow,i,n);
    end
end
for i = 1:1:m
    for j = 1:1:n
        battery_kJ(i,j) = battery(n*(i-1)+j);
        weight.fuel(i,j) = fuel(n*(i-1)+j);
    end
end
weight.fuel = weight.fuel*1.15; %Total LH2 storage including weight of insulated container
weight.battery = battery_kJ/1260; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
weight.total = weight.sofc + weight.otm + weight.comp + weight.turb + weight.hx + weight.motor + weight.battery + weight.propulsor + weight.fuel; 
weight.payload = 820000 - weight.total;
param.weight = weight;
param.P_den = param.NetPower./(weight.sofc + weight.otm + weight.comp + weight.turb + weight.hx);
end%Ends function run_cycle

function [fuel,battery] = flight_profile(options,mission,vol_flow,par_i,n)
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
vol_flow2 = vol_flow(i,j)*[1;.9;.8;.7;.5; .5*ones(mm-5,1);]*ones(1,nn);%reduce volume flow to 50%, then increase P_perm to reduce oxygen and power
options2.height = ones(mm,1)*mission.alt'; %Altitude, meters
air_den = interp1(alt_tab,atmosphere_density,options2.height);
molar_flow2 = vol_flow2.*air_den/28.84;%Flow rate at altitude assuming constant volumetric flow device
[A1,ss] = std_atmosphere(options2.height,molar_flow2);%Ambient conditions as a function of altitude
for k = 1:1:nn
    options2.P_perm(:,k) = [50*ones(5,1);logspace(log10(50),log10(0.99*.21*A1.P(1,k)*options.PR_comp(i,j)),mm-5)']; %Pressure of OTM oxygen stream, kPa; 
end
[OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options2,A1);
%adjust permeate pressure to be within feasible oxygen output range for SOFC area
min_O2 = 0.1*options2.SOFC_area(1,1)*10000/(96485.33*4000);
max_O2 = .4./options2.asr(1,1).*options2.SOFC_area(1,1)*10000/(96485.33*4000);%Cant solve for ultra low or high current densities
if any(any(O5.O2>max_O2)) || any(any(O5.O2<min_O2))
    R1 = max_O2./O5.O2;
    R2 = min_O2./O5.O2;
    options2.P_perm = min(max(options2.P_perm,(0.21*options2.PR_comp.*A1.P).^(1-R1).*options2.P_perm.^R1),(0.21*options2.PR_comp.*A1.P).^(1-R2).*options2.P_perm.^R2);
    for k = 1:1:nn
        options2.P_perm(:,k) = linspace(min(options2.P_perm(:,k)),max(options2.P_perm(:,k)),mm);
    end
    [OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options2,A1);
end
[FC,E1] = oxy_fuelcell(options2,O5);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options2,FC,OTM,E1);
P_shaft = options2.motor_eff.*(FC.Power + OTM.net_work + HL.blower_work);
%find permeate pressure condition that results in correct power for each flight segment
for k = 1:1:nn
    %update thrust to include climb and acceleration
    T = options.TO_weight(i,j)./options.Lift_2_Drag(i,j)*4.448;% Thrust in N
    P_req = T*mission.mach_num(k)*ss(1,k)/options.prop_eff(i,j)/1000;%shaft power in kW.  
    if P_req>max(P_shaft(:,k))
        [P,I] = max(P_shaft(:,k));
        battery = battery + (P_req - P)*mission.duration(k)*3600;
        fuel = fuel + FC.H2_used(I,k)*2*mission.duration(k)*3600;
    elseif P_req<min(P_shaft(:,k))
        fuel = fuel + P_req/min(P_shaft(:,k))*min(FC.H2_used(:,k));
    else
        fuel = fuel + interp1(P_shaft(:,k),FC.H2_used(:,k),P_req)*2*mission.duration(k)*3600;
    end
end

end%Ends function flight_profile