
function [performancetable,weighttable,missionmax] = collector_std_cycle(param,mission)
dp = mission.design_point; 
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
performancetable.dp = dp; 
performancetable.tsfc = zeros(v,1);
performancetable.thrust = zeros(v,1);
for i =1:v
    performancetable.tsfc(i,1) = param.TSFC_mission(x,y,i);
    performancetable.thrust(i,1) = mission.thrust(x,y,i); 
end
performancetable.PASTEto = [performancetable.topower,performancetable.FCVto,performancetable.toefficiency,performancetable.FCidento]'; 
performancetable.PASTEdp = [performancetable.dppower,performancetable.FCVcr,performancetable.dpefficiency,performancetable.FCidencruise]';
weighttable.PASTE = [weighttable.payload,weighttable.sofc,weighttable.hx,weighttable.comp,weighttable.turb,weighttable.motor,weighttable.fuel_stored,weighttable.battery,weighttable.fuel_burn]';