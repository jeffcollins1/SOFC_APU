
function [Battery] = battery(options)
[Demand,Ppropmax,P_to,TrTO,time_to,climb.velocityclimb,time_climb,climb.P,climb.T,climb.d,M,climb.PrC,climb.altitude,Cruise.dcruise,Cruise.PCruise,Cruise.TCruise,Cruise.Vcruise,t1,t2,t3,t4,t5,t6,t7,t8,Fnet_to,Cruise.TrCruise,clcruise,cdcruise] = craft(options);
Battery.Pbatt1 = (Demand.Pr_to - Demand.PrCruise)*time_to; 
Battery.Pbatt2 = (Demand.PrCb1 - Demand.PrCruise)*t1;
Battery.Pbatt3 = (Demand.PrCb2 - Demand.PrCruise)*t2;
Battery.Pbatt4 = (Demand.PrCb3 - Demand.PrCruise)*t3;
Battery.Pbatt5 = (Demand.PrCb4 - Demand.PrCruise)*t4;
Battery.Pbatt6 = (Demand.PrCb5 - Demand.PrCruise)*t5;
Battery.Pbatt7 = (Demand.PrCb6 - Demand.PrCruise)*t6;
Battery.Pbatt8 = (Demand.PrCb7 - Demand.PrCruise)*t7
Battery.Pbatt9 = (Demand.PrCb8- Demand.PrCruise)*t8;
Battery.Total = Battery.Pbatt1 + Battery.Pbatt2 + Battery.Pbatt3 + Battery.Pbatt4 + Battery.Pbatt5 + Battery.Pbatt6 + Battery.Pbatt7 + Battery.Pbatt8 + Battery.Pbatt9; 
end






