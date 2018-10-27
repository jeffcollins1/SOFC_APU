function [param,Air,FC] = run_std_cycle(options,mission,res_fuel)
%% calculate conditions for each design trade-off per kmol/s of inlet air
[m,n] = size(options.SOFC_area);
molar_flow = ones(m,n);
options.height = mission.alt(mission.design_point)*ones(m,n); %Altitude, meters
alt_tab = [0:200:7000,8000,9000,10000,12000,14000];%
atmosphere_density = [1.225,1.202,1.179,1.156,1.134,1.112,1.090,1.069,1.048,1.027,1.007,0.987,0.967,0.947,0.928,0.909,0.891,0.872,0.854,0.837,0.819,0.802,0.785,0.769,0.752,0.736,0.721,0.705,0.690,0.675,0.660,0.646,0.631,0.617,0.604,0.590,0.526,0.467,0.414,0.312,0.228]'; %Density, kg/m^3
air_den = interp1(alt_tab,atmosphere_density,options.height);
[Amb,ss] = std_atmosphere(options.height,molar_flow);%Ambient conditions as a function of altitude
A1 = Amb; 
A1.P = Amb.P.*((1+ 0.5*0.4.*mission.mach_num(mission.design_point).^2).^(1/0.4));
A1.T = Amb.T.*(1 + 0.5*0.4.*mission.mach_num(mission.design_point).^2);
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
A3 = A2;
%A3.T = options.T_fc - 0.5*options.dT_fc;
[FC,A4,E1,F5,A3,A5] = des_fuelcell(options,A2);
A5 = A4;
RC.Q = enthalpy(A3) - enthalpy(A2); %property(A3,'h','kJ') - property(A2,'h','kJ');
% A5.T = find_T(A5,enthalpy(A4) - RC.Q);%find_T(A5,property(A4,'h','kJ') - RC.Q); %TIT if heat is recovered from inlet, or the resulting outlet temperature from the heat pipe if considered

[A6,T1] = expander(A5,A1.P,options.T1_eff);
[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1,A5);

%% Calculate nominal and mission power
mission.air_den = interp1(alt_tab,atmosphere_density,mission.alt);
[~,mission.ss] = std_atmosphere(mission.alt,1);%Ambient conditions as a function of altitude
for i = 1:1:length(mission.alt)
    velocity = mission.mach_num(i).*mission.ss(i);
    Thrust_Coefficient = mission.thrust(:,:,i)./(options.num_engines.*.5*mission.air_den(i).*velocity.^2*pi()*options.engine_radius.^2);
    Froude_efficiency = 2./(1+(1+Thrust_Coefficient).^.5);
    net_prop_eff = options.prop_eff.*Froude_efficiency;
    mission.power(:,:,i) = mission.thrust(:,:,i)*velocity./net_prop_eff/1000 + options.electric_demand;%shaft power in kW. 
end
P_nominal = mission.power(:,:,mission.design_point);%nominal power in kW
dp_density = mission.air_den(mission.design_point);
min_density = min(mission.air_den);
%% scale system to meet nominal power requirements
scale = options.safety_factor.*P_nominal./options.motor_eff./(FC.Power + C1.work + T1.work + B1.work);
% cath_out_ideal = A4;
% cath_out_ideal_loop = cath_out_ideal; 
% cath_out_ideal_loop.T = A5.T; 
% [y,z] = size(A4.T);
% for u = 1:y %Find scale for mass flow to bring TIT below physical limit
%     for v = 1:z
%     if cath_out_ideal_loop.T(u,v) > 1048
%     cath_out_ideal_loop.N2(u,v) = 1.2*(cath_out_ideal_loop.T(u,v)/1048)*cath_out_ideal_loop.N2(u,v);
%     cath_out_ideal_loop.O2(u,v) = 1.2*(cath_out_ideal_loop.T(u,v)/1048)*cath_out_ideal_loop.O2(u,v);
%     flow_new.N2 = cath_out_ideal_loop.N2(u,v) ;
%     flow_new.O2 = cath_out_ideal_loop.O2(u,v) ;
%     H_in = FC.H_cath_out(u,v);
%     flow_new.T = 1000;
%     flow_new.T = find_T(flow_new,H_in);
%     cath_out_ideal_loop.T(u,v) = flow_new.T;
%     end
%     end
% end
[w1,w2] = max(max(mission.alt)); 
[w3,w4] = find(mission.alt ==w1);
point = mission.design_point/w3;
%if point <= 0.2
Air.scalemax= 2; %Maximum flow rate set for high altitude
Air.scalemin = 1;
% else if point > 0.2 && point <= 0.4
%         maxboost = 0.8;
%         Air.scalemax = 1+ maxboost;
%         Air.scalemin = 1- (1-maxboost);
%     else if point > 0.4 && point <= 0.6
%             maxboost = 0.6;
%             Air.scalemax = 1+maxboost;
%         Air.scalemin = (1-(1-maxboost));
%         else if point > 0.6 && point <= 0.8
%                 maxboost = 0.2;
%             Air.scalemax = 1+ maxboost;
%         Air.scalemin = 1-(1-maxboost);
%             else if point > 0.8
%                     Air.scalemax = 1;
%                     Air.scalemin = 0.5;
%                 end
%             end
%         end
%     end
% end
  Air.scale2 = dp_density/min_density;            
molar_flow = (scale).*molar_flow;
vol_flow = molar_flow*28.84./air_den;%Volumetric flow at the design condition
options.SOFC_area = scale.*options.SOFC_area;

%% Re-Run with scaled system parameters
[Amb,~] = std_atmosphere(options.height,molar_flow);%Ambient conditions as a function of altitude
A1 = Amb; 
A1.P = Amb.P.*((1+ 0.5*0.4.*mission.mach_num(mission.design_point).^2).^(1/0.4));
A1.T = Amb.T.*(1 + 0.5*0.4.*mission.mach_num(mission.design_point).^2);
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
%A3 = A2;
%A3.T = options.T_fc - 0.5*options.dT_fc;
[FC,A4,E1,F5,A3,A5] = des_fuelcell(options,A2);
A5 = A4;
Air.relation = FC.i_den./net_flow(A2);
%A5.T = A4.T - (A3.T - A2.T);
RC.Q = enthalpy(A3) - enthalpy(A2);
A5.T = find_T(A5,enthalpy(A4) - RC.Q);
[A6,T1] = expander(A5,A1.P,options.T1_eff);
[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1,A5);
HX.HP.mass = FC.Qremove.*options.heat_pipe_specific_mass; 
weight = system_weight(options,FC,{C1;T1;B1;},[],HX);
param = NetParam(options,FC,{C1;T1;B1;},[],FL);
param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5,;'A6',A6;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F2',F2;'F3',F3;'F4',F4;'F5',F5;};

%% calculate off-design power output to meet mission profile for each condition by varying permeate pressure
nn = length(mission.alt);
fuel = zeros(m*n,nn);
battery = zeros(m*n,nn);
battery_kJ = zeros(m,n,nn);
P_sys_mission = zeros(m*n,nn);
eff_mission = zeros(m*n,nn);
FCV_mission = zeros(m*n,nn);
TIT = zeros(m*n,nn);
FCiden_mission = zeros(m*n,nn);
TSFC_mission = zeros(m*n,nn);
bypass = zeros(m*n,nn); 
param.power_mission = zeros(m,n,nn);
param.efficiency_mission = zeros(m,n,nn);
param.FCV_mission = zeros(m,n,nn);
param.FCiden_mission = zeros(m,n,nn);
param.TSFC_mission = zeros(m,n,nn);
param.turbine_inlet = zeros(m,n,nn);
param.FCbypass = zeros(m,n,nn);
parallel = true;
if parallel
    parfor par_i = 1:1:m*n
        [fuel(par_i,:),battery(par_i,:),P_sys_mission(par_i,:),eff_mission(par_i,:),FCV_mission(par_i,:),FCiden_mission(par_i,:),TSFC_mission(par_i,:),TIT(par_i,:),bypass(par_i,:)] = flight_profile(options,mission,vol_flow,F5.H2,par_i,n,Air);
    end
else
    for i = 1:1:m*n
        [fuel(i,:),battery(i,:),P_sys_mission(i,:),eff_mission(i,:),FCV_mission(i,:),FCiden_mission(i,:),TSFC_mission(i,:),TIT(i,:),bypass(i,:)] = flight_profile(options,mission,vol_flow,F5.H2,i,n,Air);
    end
end
for i = 1:1:m
    for j = 1:1:n
        battery_kJ(i,j,:) = battery(n*(i-1)+j,:);
        param.fuel_by_seg(i,j,:) = fuel(n*(i-1)+j,:);
        param.power_mission(i,j,:) = P_sys_mission(n*(i-1)+j,:);
        param.efficiency_mission(i,j,:) = eff_mission(n*(i-1)+j,:);
        param.FCV_mission(i,j,:) = FCV_mission(n*(i-1)+j,:);
        param.FCiden_mission(i,j,:) = FCiden_mission(n*(i-1)+j,:);
        param.TSFC_mission(i,j,:) = TSFC_mission(n*(i-1)+j,:);
        param.battery_mass_by_segment(i,j,:) = battery_kJ(i,j,:)./options.battery_specific_energy(i,j); 
        param.turbine_inlet(i,j,:) = TIT(n*(i-1)+j,:);
        param.FCbypass(i,j,:) = bypass(n*(i-1)+j,:);
    end
end
weight.fuel_burn = sum(param.fuel_by_seg,3); 
weight.fuel_stored = weight.fuel_burn.*options.fuel_tank_mass_per_kg_fuel + res_fuel/3; %Total LH2 storage including weight of insulated container and equivalent energy reserve storage
weight.battery = sum(battery_kJ,3)./options.battery_specific_energy; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
weight.total = (weight.sofc + weight.comp + weight.turb + weight.hx + weight.motor + weight.battery + weight.propulsor + weight.fuel_stored); 
param.weight = weight;
param.P_den = param.NetPower./(weight.sofc + weight.comp + weight.turb + weight.hx);

end%Ends function run_cycle

function [fuel,battery,P_sys_mission,eff_mission,FCV_mission,FCiden_mission,TSFC_mission,turbine_inlet_temperature,bypass] = flight_profile(options,mission,vol_flow,nominal_fuel,par_i,n,Air,FC)
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
battery = zeros(1,nn);
fuel = zeros(1,nn);
dp1 = mission.design_point; 
vol_flow2 = vol_flow(i,j)*ones(mm,nn);%reduce volume flow to 50%
options2.height = ones(mm,1)*mission.alt'; %Altitude, meters
air_den = interp1(alt_tab,atmosphere_density,options2.height);
%molar_flow2 = vol_flow2.*air_den/28.84;%Flow rate at altitude assuming constant volumetric flow device
[s1,s2] = size(air_den);
max_air_den = max(max(air_den));
min_air_den = min(min(air_den));
dp_air_den = air_den(mission.design_point); 
[c1,c2] = find(air_den == max_air_den);
[c3,c4] = find(air_den==min_air_den);
descent = s1 - c4(1,1);
comp_speed1 = linspace(0.5,1,c4(1,1));
comp_speed2 = linspace(1,0.5,descent);
comp_speed = 2*ones(s1,1)*[comp_speed1,comp_speed2];%.*ones(s1,1);
vol_flow3 = vol_flow2.*(comp_speed)./(comp_speed(1,dp1));
molar_flow2 = vol_flow3.*air_den/28.84;%Flow rate at altitude assuming constant volumetric flow device
[Amb,~] = std_atmosphere(options2.height,molar_flow2);%Ambient conditions as a function of altitude
A1 = Amb; 
mach = ones(nn,1).*mission.mach_num';
A1.P = Amb.P.*((1+ 0.5*0.4.*mach.^2).^(1/0.4));
A1.T = Amb.T.*(1 + 0.5*0.4.*mach.^2);
[A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);
A3 = A2;
%% add bypass
O2_util_max = 0.45;
O2_util_min = 0.10;
max_current = O2_util_max.*A3.O2*(96485.33*4000);% 45% O2 utilization
min_current = O2_util_min.*A3.O2*(96485.33*4000);% 1% O2 utilization
max_iden = max_current./(1e4*options2.SOFC_area);
% if max(max(max_iden)) > 4
%     O2_util_max2 = (4./max_iden).*O2_util_max;
%     O2_util_min2 = (4./max_iden).*O2_util_min;
%     max_current = O2_util_max2.*A3.O2*(96485.33*4000);% 45% O2 utilization
%     min_current = O2_util_min2.*A3.O2*(96485.33*4000);% 1% O2 utilization
% end
% max_iden = max_current./(1e4*options2.SOFC_area);
%     max_fuel = max_current./(96485.33*2000)./options2.spu;
%     min_fuel = min_current./(96485.33*2000)./options2.spu;
max_fuel = min(nominal_fuel(i,j)*1.2,max_current./(96485.33*2000)./options2.spu);
min_fuel = max(nominal_fuel(i,j)*.25,min_current./(96485.33*2000)./options2.spu);
F5.T = options2.T_fc - .5*options2.dT_fc;
F5.P = A3.P;
F5.H2 = zeros(mm,nn);
for k = 1:1:nn
    F5.H2(:,k) = linspace(max_fuel(1,k),min_fuel(mm,k),mm)';
end
F5.H2O = F5.H2.*options2.steamratio./(1-options2.steamratio);
molar_flow2 = A3.N2 + A3.O2;
[Amb,~] = std_atmosphere(options2.height,molar_flow2)  ;  
 A1 = Amb; 
mach = ones(nn,1).*mission.mach_num';
A1.P = Amb.P.*((1+ 0.5*0.4.*mach.^2).^(1/0.4));
A1.T = Amb.T.*(1 + 0.5*0.4.*mach.^2);
[A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);

[FC,A4,E1] = std_fuelcell(options2,A2,F5);
%molar_flow3 = mean(mean(max_iden))./mean(mean(Air.relation));
%m = mean(mean(FC.bypass));

%  [Amb,~] = std_atmosphere(options2.height,molar_flow2)    
%  A1.N2 = Amb.N2.*(1-0.5*FC.bypass); 
%  A1.O2 = Amb.O2.*(1-0.5*FC.bypass); 
% mach = ones(nn,1).*mission.mach_num';
% A1.P = Amb.P.*((1+ 0.5*0.4.*mach.^2).^(1/0.4));
% A1.T = Amb.T.*(1 + 0.5*0.4.*mach.^2);
% [A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);
% [FC,A4,E1] = std_fuelcell(options2,A2,F5);
 A4.T = A4.T;


[A5,T1] = expander(A4,A1.P,options2.T1_eff);

[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options2,E1,F5,A1);
P_sys = FC.Power + C1.work + T1.work + B1.work;
P_shaft = options2.motor_eff.*P_sys;
FTE = P_sys./(FC.H2_used.*FC.hrxnmol);
P_sys_mission = zeros(1,nn);
eff_mission = zeros(1,nn);
FCV_mission = zeros(1,nn);
FCiden_mission = zeros(1,nn);
turbine_inlet_temperature = zeros(1,nn);
bypass = zeros(1,nn);
%find permeate pressure condition that results in correct power for each flight segment
for k = 1:1:nn
    P_req = mission.power(i,j,k);%shaft power in kW.  
    if P_req>max(P_shaft(:,k))
        [P,I] = max(P_shaft(:,k));
        battery(k) = (P_req - P)*mission.duration(k)*3600;
        fuel(k) = FC.H2_used(I,k)*2*mission.duration(k)*3600;
        P_sys_mission(k) = P_sys(I,k);
        eff_mission(k) = FTE(I,k);
        FCV_mission(k) = FC.V(I,k);
        FCiden_mission(k) = FC.i_den(I,k);
        bypass(k) = FC.bypass(I,k);
        turbine_inlet_temperature(k) = A4.T(I,k);
    elseif P_req<min(P_shaft(:,k))
        [~,I] = min(P_shaft(:,k));
        fuel(k) = P_req/P_shaft(I,k)*(FC.H2_used(I,k))*2*mission.duration(k)*3600;
        P_sys_mission(k) = P_req/min(P_sys(:,k))*P_sys(I,k);
        eff_mission(k) = FTE(I,k);
        FCV_mission(k) = FC.V(I,k);
        FCiden_mission(k) = FC.i_den(I,k);
         bypass(k) = FC.bypass(I,k); 
        turbine_inlet_temperature(k) = A4.T(I,k);
       
    else
        fuel(k) = interp1(P_shaft(:,k),FC.H2_used(:,k),P_req)*2*mission.duration(k)*3600;
        P_sys_mission(k) = interp1(P_shaft(:,k),P_sys(:,k),P_req);
        eff_mission(k) = interp1(P_shaft(:,k),FTE(:,k),P_req);
        FCV_mission(k) = interp1(P_shaft(:,k),FC.V(:,k),P_req);
        FCiden_mission(k) =interp1(P_shaft(:,k),FC.i_den(:,k),P_req);
        bypass(k) = interp1(P_shaft(:,k),FC.bypass(:,k),P_req); 
        turbine_inlet_temperature(k) = interp1(P_shaft(:,k),A4.T(:,k),P_req); 
    end
end
TSFC_mission = fuel./(squeeze(mission.thrust(i,j,:)).*mission.duration)'; % SFC in kg/N*hour; 
end  %Ends function flight_profile