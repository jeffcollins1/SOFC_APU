function [param,FC] = run_std_cycle(options,mission,res_fuel)
[m,n] = size(options.SOFC_area);
%% determine air inlet from flight conditions
options.height = mission.alt(mission.design_point)*ones(m,n); %Altitude, meters
ambient = std_atmosphere(options.height);%Ambient conditions as a function of altitude
A1.P = ambient.P.*((1+ 0.5*0.4.*mission.mach_num(mission.design_point).^2).^(1/0.4));
A1.T = ambient.T.*(1 + 0.5*0.4.*mission.mach_num(mission.design_point).^2);
A1.O2 = .21*ones(m,n);
A1.N2 = .79*ones(m,n);

%% calculate conditions for each design trade-off per kmol/s of inlet air
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
[FC,A4,E1,F5,A3] = des_fuelcell(options,A2);
[A5,T1] = expander(A4,A1.P,options.T1_eff);
[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1);

%% Calculate nominal and mission power
air = std_atmosphere(mission.alt);%Ambient conditions as a function of altitude
mission.ss = air.ss;
mission.air_den = air.rho;
for i = 1:1:length(mission.alt)
    velocity = mission.mach_num(i).*mission.ss(i);
    Thrust_Coefficient = mission.thrust(:,:,i)./(options.num_engines.*.5*mission.air_den(i).*velocity.^2*pi().*options.engine_radius.^2);
    Froude_efficiency = 2./(1+(1+Thrust_Coefficient).^.5);%http://web.mit.edu/16.unified/www/SPRING/systems/Lab_Notes/airpower.pdf
    net_prop_eff = options.prop_eff.*Froude_efficiency;
    mission.power(:,:,i) = mission.thrust(:,:,i)*velocity./net_prop_eff/1000 + options.electric_demand;%shaft power in kW. 
end
P_nominal = mission.power(:,:,mission.design_point);%nominal power in kW
%% scale system to meet nominal power requirements
scale = options.safety_factor.*P_nominal./options.motor_eff./(FC.Power + C1.work + T1.work + B1.work);
 
%% Re-Run with scaled system parameters
vol_flow = scale*28.84./ambient.rho;%Volumetric flow at the design condition
options.SOFC_area = scale.*options.SOFC_area;
A1.O2 = .21*scale;
A1.N2 = .79*scale;
[A2,C1] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
[FC,A4,E1,F5,A3] = des_fuelcell(options,A2);
[A5,T1] = expander(A4,A1.P,options.T1_eff);
[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1);
HX.HP.mass = FC.Qremove.*options.heat_pipe_specific_mass; 
weight = system_weight(options,FC,{C1;T1;B1;},[],HX);
param = NetParam(options,FC,{C1;T1;B1;},[],FL);
param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5,;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F2',F2;'F3',F3;'F4',F4;'F5',F5;};

%% calculate off-design power output to meet mission profile for each condition by varying permeate pressure
nn = length(mission.alt);
fuel = zeros(m*n,nn);
battery = zeros(m*n,nn);
battery_kJ = zeros(m,n,nn);
P_sys_mission = zeros(m*n,nn);
eff_mission = zeros(m*n,nn);
FCV_mission = zeros(m*n,nn);
FCiden_mission = zeros(m*n,nn);
TSFC_mission = zeros(m*n,nn);
TIT = zeros(m*n,nn);
bypass = zeros(m*n,nn); 
T1_work_mission = zeros(m*n,nn);
C1_work_mission = zeros(m*n,nn);
Q_balFC_mission = zeros(m*n,nn);
param.power_mission = zeros(m,n,nn);
param.efficiency_mission = zeros(m,n,nn);
param.FCV_mission = zeros(m,n,nn);
param.FCiden_mission = zeros(m,n,nn);
param.TSFC_mission = zeros(m,n,nn);
param.turbine_inlet = zeros(m,n,nn);
param.T1_work = zeros(m,n,nn);
param.C1_work = zeros(m,n,nn); 
param.Q_balFC = zeros(m,n,nn);
param.FCbypass = zeros(m,n,nn);
parallel = true;
if parallel
    parfor par_i = 1:1:m*n
        [fuel(par_i,:),battery(par_i,:),P_sys_mission(par_i,:),eff_mission(par_i,:),FCV_mission(par_i,:),FCiden_mission(par_i,:),TSFC_mission(par_i,:),TIT(par_i,:),bypass(par_i,:),C1_work_mission(par_i,:),T1_work_mission(par_i,:),Q_balFC_mission(par_i,:)] = flight_profile(options,mission,vol_flow,F5.H2,par_i,n);
    end
else
    for i = 1:1:m*n
        [fuel(i,:),battery(i,:),P_sys_mission(i,:),eff_mission(i,:),FCV_mission(i,:),FCiden_mission(i,:),TSFC_mission(i,:),TIT(i,:),bypass(i,:),C1_work_mission(i,:),T1_work_mission(i,:),Q_balFC_mission(i,:)] = flight_profile(options,mission,vol_flow,F5.H2,i,n);
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
        param.T1_work(i,j,:) = T1_work_mission(n*(i-1) +j,:);
        param.C1_work(i,j,:) = C1_work_mission(n*(i-1) + j,:); 
        param.Q_balFC(i,j,:) = Q_balFC_mission(n*(i-1) + j,:);
    end
end
weight.fuel_burn = sum(param.fuel_by_seg,3); 
weight.fuel_stored = weight.fuel_burn*(1+res_fuel).*options.fuel_tank_mass_per_kg_fuel; %Total LH2 storage including weight of insulated container and proportional reserve storage
weight.battery = sum(battery_kJ,3)./options.battery_specific_energy; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
weight.total = (weight.sofc + weight.comp + weight.turb + weight.hx + weight.motor + weight.battery + weight.propulsor + weight.fuel_stored); 
param.weight = weight;
param.P_den = param.NetPower./(weight.sofc + weight.comp + weight.turb + weight.hx);
end%Ends function run_cycle

function [fuel,battery,P_sys_mission,eff_mission,FCV_mission,FCiden_mission,TSFC_mission,turbine_inlet_temperature,bypass,C1_work_mission,T1_work_mission,Q_balFC_mission] = flight_profile(options,mission,vol_flow,nominal_fuel,par_i,n)
F = 96485.33; %Faraday's Constant in Coulombs/mol
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

options2.height = ones(mm,1)*mission.alt'; %Altitude, meters
ambient = std_atmosphere(options2.height);%Ambient conditions as a function of altitude
M = ones(mm,1).*mission.mach_num';
A1.P = ambient.P.*((1+ 0.5*0.4.*M.^2).^(1/0.4));
A1.T = ambient.T.*(1 + 0.5*0.4.*M.^2);
A1.O2 = .21*ones(mm,nn);
A1.N2 = .79*ones(mm,nn);
[A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);

nominal_O2 = 0.21*vol_flow(i,j).*ambient.rho/28.84;%full engine speed mass flow rate
max_current = 0.5.*nominal_O2*(F*4000);% 50% O2 utilization at max air flow
min_current = 0.025.*nominal_O2*(96485.33*4000);% 5% O2 utilization at 50% air flow
max_fuel = min(nominal_fuel(i,j)*1.2,max_current./(F*2000)./options2.spu);
min_fuel = max(nominal_fuel(i,j)*.1,min_current./(F*2000)./options2.spu);
F5.T = options2.T_fc - .5*options2.dT_fc;
F5.P = A2.P;
F5.H2 = zeros(mm,nn);
for k = 1:1:nn
    F5.H2(:,k) = linspace(max_fuel(1,k),min_fuel(mm,k),mm)';
end
F5.H2O = F5.H2.*options2.steamratio./(1-options2.steamratio);

[FC,A3,A4,E1] = std_fuelcell(options2,A2,F5);
%%find speed turndown or bypass amount
speed = max(.5, A3.O2./nominal_O2);
A1.O2 = max(A3.O2,.5*nominal_O2);
A1.N2 = A1.O2*.79/.21;
[A2,C1] = compressor(A1,options2.PR_comp.*A1.P,options2.C1_eff);

%bypass and mixing
bypass_flow = A2;
bypass_flow.O2 = A2.O2 - A3.O2;
bypass_flow.N2 = bypass_flow.O2*.79/.21;
turb_in.P = A4.P;
turb_in.T = A4.T;
turb_in.O2 = A4.O2+bypass_flow.O2;
turb_in.N2 = A4.N2+bypass_flow.N2;
turb_in.T = find_T(turb_in,enthalpy(A4)+enthalpy(bypass_flow));

[A5,T1] = expander(turb_in,A1.P,options2.T1_eff);
[FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options2,E1,F5,A1);

FC.H2_used = FC.i_total./(2000.*F) + FC.Q_pre_combustor./FC.hrxnmol;

P_sys = FC.Power + C1.work + T1.work + B1.work;
P_shaft = options2.motor_eff.*(P_sys-options2.electric_demand);
FTE = P_sys./(FC.H2_used.*FC.hrxnmol);
P_sys_mission = zeros(1,nn);
eff_mission = zeros(1,nn);
FCV_mission = zeros(1,nn);
FCiden_mission = zeros(1,nn);
T1_work_mission = zeros(1,nn);
C1_work_mission = zeros(1,nn); 
turbine_inlet_temperature = zeros(1,nn);
Q_balFC_mission = zeros(1,nn);
bypass = zeros(1,nn);
error = zeros(1,nn); 
%find permeate pressure condition that results in correct power for each flight segment
for k = 1:1:nn
    P_req = mission.power(i,j,k) - options2.electric_demand(1,k);%shaft power in kW.  
    if P_req>max(P_shaft(:,k))
        [P,I] = max(P_shaft(:,k));
        battery(k) = (P_req - P)*mission.duration(k)*3600;
        fuel(k) = FC.H2_used(I,k)*2*mission.duration(k)*3600;
        P_sys_mission(k) = P_sys(I,k);
    elseif P_req<min(P_shaft(:,k))
        [~,I] = min(P_shaft(:,k));
        scale_Pow = (P_req/options2.motor_eff(1,k)+options2.electric_demand(1,k))/P_sys(I,k);
        fuel(k) = scale_Pow*FC.H2_used(I,k)*2*mission.duration(k)*3600;
        P_sys_mission(k) = scale_Pow*P_sys(I,k);
    else
        I = interp1(P_shaft(:,k),1:mm,P_req);
        fuel(k) = interp1(P_shaft(:,k),FC.H2_used(:,k),P_req)*2*mission.duration(k)*3600;
        P_sys_mission(k) = interp1(P_shaft(:,k),P_sys(:,k),P_req);
    end
    if I == mm
        I = I-1;
        r = 1;
    else
        r = I - floor(I);
        I = floor(I);
    end
    
    eff_mission(k) = (1-r)*FTE(I,k) + r*FTE(I+1,k);
    FCV_mission(k) = (1-r)*FC.V(I,k) + r*FC.V(I+1,k);
    FCiden_mission(k) = (1-r)*FC.i_den(I,k) + r*FC.i_den(I+1,k);
    bypass(k) = (1-r)*(bypass_flow.O2(I,k)/A2.O2(I,k)) + r*(bypass_flow.O2(I+1,k)/A2.O2(I+1,k)); 
    T1_work_mission(k) = (1-r)*T1.work(I,k) + r*T1.work(I+1,k);
    C1_work_mission(k) = (1-r)*C1.work(I,k) + r*C1.work(I+1,k);
    turbine_inlet_temperature(k) = (1-r)*turb_in.T(I,k) + r*turb_in.T(I+1,k);
    Q_balFC_mission(k) = (1-r)*FC.Qremove(I,k) + r*FC.Qremove(I+1,k);
    Ebalsystem(k) = fuel(k).*FC.hrxnmol(k).*eff_mission(k)/2 - P_sys_mission(k).*mission.duration(k)*3600; %Check fuel burn, convert back to kmol from kg, multiply by reported efficiency and LHV fuel, and compare with total system energy output over flight segment 
    Pbalmission(k) = P_req - P_sys_mission(k) - battery(k)/(3600*mission.duration(k)); %check instananeous balance between required power, system output and battery assist
%     if bypass(k)>0
%         disp(strcat('TIT:',num2str(turbine_inlet_temperature(k))))
%     end

%     if ((1-r)*speed(I,k) + r*speed(I+1,k))>1
%         disp(strcat('overspeed:',num2str(((1-r)*speed(I,k) + r*speed(I+1,k)-1)*100),'%'))
%     end
end
TSFC_mission = fuel./(squeeze(mission.thrust(i,j,:)).*mission.duration)'; % SFC in kg/N*hour; 
end  %Ends function flight_profile