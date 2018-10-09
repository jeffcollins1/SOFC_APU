function [performancetable,weighttable] = collector(param,mission)
dp = mission.design_point; 
maxpayload = max(max(param.weight.payload));
[x,y] = find(param.weight.payload == maxpayload);
weighttable.payload = maxpayload; 
weighttable.sofc = param.weight.sofc(x,y);
weighttable.otm = param.weight.otm(x,y);
weighttable.hx = param.weight.hx(x,y);
weighttable.comp = param.weight.comp(x,y);
weighttable.turb = param.weight.turb(x,y);
weighttable.motor = param.weight.motor(x,y);
weighttable.fuel = param.weight.fuel(x,y);
weighttable.battery = param.weight.comp(x,y);
performancetable.topower = param.power_mission(x,y,1);
performancetable.dppower = param.power_mission(x,y,dp);
performancetable.dpefficiency = param.efficiency_mission(x,y,dp); 
performancetable.toefficiency = param.efficiency_mission(x,y,1); 
performancetable.FCidento = param.FCiden_mission(x,y,1);
performancetable.FCidencruise = param.FCiden_mission(x,y,dp); 
performancetable.FCVto = param.FCV_mission(x,y,1);
performancetable.FCVcr = param.FCV_mission(x,y,dp);
performancetable.dp = dp; 
end
