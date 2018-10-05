function Eff = efficiency(FC,OTM,HL)
Ein = FC.H2_used*240000; %Lower heating value of fuel intake, kJ/kmol
Eout = FC.Power + OTM.net_work + HL.blower_work;
Eff.FTE = Eout./Ein;
Eff.FC = FC.Power./Ein;
end
