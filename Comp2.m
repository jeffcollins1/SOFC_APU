
function [QO1_2,WC2,TO2] = Comp2(TOperation,PO1,PO3,CompEff,flowO2FC)
if PO1 > 0
TO3 = 1; %Initialize TO3
TO2 = 270;  %Initialize TO2
TO1 = TOperation;
PO2 = PO1;
HO1 = refpropm('H','T',TO1,'P',PO1,'OXYGEN')/1000; %Enthalpy of permeate O2 stream
while TO3<TOperation;
    TO2= TO2 + 1
HO2 = refpropm('H','T',TO2,'P',PO1,'OXYGEN')/1000; 
sO2 = refpropm('S','T',TO2,'H',HO2*1000,'OXYGEN'); 
HO3s = refpropm('H','P',PO3,'S',sO2,'OXYGEN')/1000;
HO3a = (HO3s -HO2)/CompEff + HO2; 
TO3 = refpropm('T','P',PO3,'H',HO3a*1000,'OXYGEN')
end
TO2; 
TO3;
mflow = flowO2FC*0.032; %Air flow in kg/s
QO1_2 = mflow*(HO1 - HO2); %Heat rejected to reach desired compressor inlet temperature
WC2 = mflow*(HO3a - HO2); 
else
    QO1_2 = 0;
    WC2 = 0; 
    mflow = 0;
    TO2 = 0; 
end
end