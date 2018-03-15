function [WBlower,TF4,HF4] = Blower(mFlowH2in,Pint4)
PE3 = Pint4 - 200; 
TE3 = 420;
EffBlower = 0.5;
PE4 = Pint4;
HE3 = refpropm('H','T',TE3,'P',PE3,'HYDROGEN');
sE3 = refpropm('S','T',TE3,'P',PE3,'HYDROGEN');
HE4s = refpropm('H','T',TE3,'S',sE3,'HYDROGEN');
HF4 = HE3 + (HE4s - HE3)/EffBlower;
TF4 = refpropm('T','P',PE4,'H',HE3,'HYDROGEN');
WBlower = mFlowH2in*(HF4-HE3); 