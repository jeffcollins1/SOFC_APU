function [performancetable,weighttable] = collector(param,mission)

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
%Q_diff = zeros(10,10,dp);
% for q = 1:dp
% mission_Q_out = param.Q_balFC(:,:,q);
% Q_diff(:,:,q) = mission_Q_out - design_Q_out;
% end
% [x1,y1,z1] = find(param.Q_balFC == max_Q_out1);
% max_Q_out = max(param.Q_balFC(:,:,z1)); 
% param.weight.hx = (max_Q_out-design_Q_out)./1.72 + param.weight.hx;
% param.weight.payload = param.weight.payload - (max_Q_out-design_Q_out)/1.72; 
% param.weight.hx;
[v,w] = size(mission.alt); 
maxpayload = max(max(param.weight.payload));
[x,y] = find(param.weight.payload == maxpayload);
weighttable.payload = maxpayload; 
weighttable.sofc = param.weight.sofc(x,y);
weighttable.otm = param.weight.otm(x,y);
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
performancetable.P_in = param.P_perm(x,y,1); 
performancetable.T1_workto = param.T1_work(x,y,1);
performancetable.T1_workdp = param.T1_work(x,y,dp);
performancetable.C1_workto = param.C1_work(x,y,1);
performancetable.C1_workdp = param.C1_work(x,y,dp);
% performancetable.P_reqto = mission.power(x,y,1);
% performancetable.P_reqdp = mission.power(x,y,dp); 
performancetable.QbalFCto = param.Q_balFC(x,y,1); 
performancetable.QbalFCdp = param.Q_balFC(x,y,dp); 
performancetable.Qbalmotorto = param.Q_balmotors(x,y,1); 
performancetable.Qbalmotordp = param.Q_balmotors(x,y,dp); 
performancetable.dp = dp; 
performancetable.tsfc = zeros(v,1);
performancetable.thrust = zeros(v,1);
for i =1:v
    performancetable.tsfc(i,1) = param.TSFC_mission(x,y,i);
    performancetable.thrust(i,1) = mission.thrust(x,y,i); 
    performancetable.missionefficiency(i,1) =param.efficiency_mission(x,y,i);
    performancetable.QbalFC_mission(i,1) = param.Q_balFC(x,y,i);
    performancetable.i_den_mission(i,1) = param.FCiden_mission(x,y,i); 
end
% design_Q_out = performancetable.QbalFC_mission(dp,1);
% max_Q_out = max(max(performancetable.QbalFC_mission));
% 
% weighttable.hxTrue = (max_Q_out-design_Q_out)/1.72 + weighttable.hx;
weighttable.payloadTrue = weighttable.payload - (max_Q_out-design_Q_out)/1.72; 
performancetable.PASTEto = [performancetable.topower,performancetable.FCVto,performancetable.toefficiency,performancetable.FCidento,performancetable.T1_workto,performancetable.C1_workto,performancetable.QbalFCto,performancetable.Qbalmotorto]'; 
performancetable.PASTEdp = [performancetable.dppower,performancetable.FCVcr,performancetable.dpefficiency,performancetable.FCidencruise,performancetable.T1_workdp,performancetable.C1_workdp,performancetable.QbalFCdp,performancetable.Qbalmotordp]';
weighttable.PASTE = [weighttable.payload,weighttable.sofc,weighttable.otm,weighttable.hx,weighttable.comp,weighttable.turb,weighttable.motor,weighttable.fuel_stored,weighttable.battery,weighttable.fuel_burn]';
end
