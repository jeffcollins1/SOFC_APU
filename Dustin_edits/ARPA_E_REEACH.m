
segment.power = [100,70,50,35,35,30,100,70,35,30]/100;
segment.time = [1/12,1/6,1/4,1/3,5,1/2,1/12,1/6,1/2,1/2];
segment.alt = [100 500 1000 6000 11000 5000 100 500 2000 2000]; %altitude in m
segment.mach_num = [0 .1 .3 .5 .85 .6 .4 .2 .4 .2];
estimated_parasitic = [2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000]; %just copy net_turbomachinery_work and iterate

F = 96485;
fuel_kmol_per_kg = 5.7604/1000; %fuel property
fuel_energy = 44.0*1000; %kJ/kg
net_utilization = 0.85;

max_power = 26000;
bat_power = 12000;
recharge = 200;
max_fc_power = max_power - bat_power + max(estimated_parasitic);% + 2000 to cover parasitic
FC_power = min(max_fc_power,max_power*segment.power+estimated_parasitic);
FC_power(5) = FC_power(5) + recharge;

OCV = 1.1;
V_at_max_pow = 0.8;
V = OCV - FC_power/max(FC_power)*(OCV-V_at_max_pow);
Air = std_atmosphere(segment.alt);%Ambient conditions as a function of altitude


O2_use = FC_power./V/(4*F); %kmol/s O2
H2_use = 2*O2_use;
H2mol_per_fuel_mol = 12.4+1.8*12.4;
fuel_mass = H2_use/H2mol_per_fuel_mol/fuel_kmol_per_kg/net_utilization;%kg/s


C1_PR = 4;
C1_eff = 0.75;
C2_PR = 8;
C2_eff = 0.75;
T1_eff = .85;

mass_flow = 28.84/.21/.3*O2_use;%kg/s

C1_inlet.P = Air.P.*((1+ 0.5*0.4.*segment.mach_num.^2).^(1/0.4));
C1_inlet.T = Air.T;

C1_inlet.O2 = mass_flow/28.84*0.21;
C1_inlet.N2 = mass_flow/28.84*0.79;
[C1_outlet,C1_prop] = compressor(C1_inlet,C1_inlet.P*C1_PR,C1_eff);
C2_inlet = C1_outlet;
C2_inlet.T = 0*C2_inlet.T + 273; %Intercooler before scroll
C2_outlet = C2_inlet;
C2_outlet.P = C2_inlet.P*C2_PR;
%Isothermal compression
pV = mass_flow/28.84*Ru.*C2_inlet.T;
C2_prop.work = pV*log(1/C2_PR)*1/C2_eff;

FC_inlet = C2_outlet;
FC_inlet.T = 0*FC_inlet.T + 650+273;
air_preheat = enthalpy(FC_inlet) - enthalpy(C2_outlet);

%expander
expander_inlet = C2_outlet;
expander_inlet.O2 = C2_outlet.O2 - O2_use;
expander_inlet.T = 250+273;
expander_inlet.P = C2_outlet.P - 25; %25 psi drop across stack

[expander_outlet,T1_prop] = expander(expander_inlet,C1_inlet.P,T1_eff);
net_turbomachinery_work = C1_prop.work + C2_prop.work + T1_prop.work;

bat_power = max_power*segment.power - (net_turbomachinery_work + FC_power);
net_power = net_turbomachinery_work + FC_power + bat_power;
efficiency = (net_turbomachinery_work + FC_power)./(fuel_mass*fuel_energy);

bat_energy = sum(bat_power(1:4).*segment.time(1:4))%kWh
total_fuel = sum(fuel_mass.*segment.time*3600)%kg



