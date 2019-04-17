function [design_point,take_off,component_weight] = collector_std_cycle(param,mission)
dp = mission.design_point; 
design_Q_out = param.Q_balFC(:,:,dp);
max_Q_out = param.Q_balFC(:,:,1);
Q_diff  = max_Q_out - design_Q_out;
for z1 = 1:10
    for z2 = 1:10
        if Q_diff(z1,z2) > 0
            param.weight.hx(z1,z2) = param.weight.hx(z1,z2) + Q_diff(z1,z2)/1.72; %add difference of heat pipe for largest heat rejection
            param.weight.payload(z1,z2) = param.weight.payload(z1,z2) - Q_diff(z1,z2)/1.72; %Subtract heat pipe mass from payload
        end
    end
end
[v,w] = size(mission.alt); 
maxpayload = max(max(param.weight.payload));
[x,y] = find(param.weight.payload == maxpayload);
weighttable.payload = maxpayload; 
weighttable.sofc = param.weight.sofc(x,y);
%weighttable.otm = param.weight.otm(x,y);
weighttable.hx = param.weight.hx(x,y);
weighttable.comp = param.weight.comp(x,y);
weighttable.turb = param.weight.turb(x,y);
weighttable.motor = param.weight.motor(x,y);
weighttable.fuel_stored = param.weight.fuel_stored(x,y);
weighttable.fuel_burn = param.weight.fuel_burn(x,y);
weighttable.battery = param.weight.comp(x,y);
performancetable.topower = param.power_mission(x,y,1);
performancetable.dppower = param.power_mission(x,y,dp);
performancetable.dpefficiency = param.efficiency_mission(x,y,dp); 
performancetable.toefficiency = param.efficiency_mission(x,y,1); 
performancetable.FCidento = param.FCiden_mission(x,y,1);
performancetable.FCidencruise = param.FCiden_mission(x,y,dp); 
performancetable.FCVto = param.FCV_mission(x,y,1);
performancetable.FCVcr = param.FCV_mission(x,y,dp);
performancetable.P_den = param.P_den; 
%performancetable.P_in = param.P_perm(x,y,1); 
performancetable.T1_workto = param.T1_work(x,y,1);
performancetable.T1_workdp = param.T1_work(x,y,dp);
performancetable.C1_workto = param.C1_work(x,y,1);
performancetable.C1_workdp = param.C1_work(x,y,dp);
performancetable.dp = dp; 
performancetable.tsfc = zeros(v,1);
performancetable.thrust = zeros(v,1);
performancetable.TIT = zeros(v,1);
performancetable.mission_efficiency = zeros(v,1);
for i = 1:v
    performancetable.tsfc(i,1) = param.TSFC_mission(x,y,i);
    performancetable.thrust(i,1) = mission.thrust(x,y,i); 
    performancetable.TIT(i,1) = param.turbine_inlet(x,y,i); 
    performancetable.mission_efficiency(i,1) = param.efficiency_mission(x,y,i);
performancetable.I_den(i,1) = param.FCiden_mission(x,y,i);
performancetable.QbalFC_mission(i,1) = param.Q_balFC(x,y,i);
%     performancetable.FCV(i,1) = param.FCV_mission(x,y,i); 
end
take_off = [performancetable.topower,performancetable.FCVto,performancetable.toefficiency,performancetable.FCidento,performancetable.T1_workto,performancetable.C1_workto]'; 
design_point = [performancetable.dppower,performancetable.FCVcr,performancetable.dpefficiency,performancetable.FCidencruise,performancetable.T1_workdp,performancetable.C1_workdp]';
component_weight = [weighttable.payload,weighttable.sofc,weighttable.hx,weighttable.comp,weighttable.turb,weighttable.motor,weighttable.fuel_stored,weighttable.battery,weighttable.fuel_burn]';