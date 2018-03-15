function [flowO2FC,mFlowH2in,Qgen,PFC,QE1_2,QE2_3,QE4_F3,Cells,TF2,HF2,TE2,HE2,TE3,HE3,TE4,HE4,TF3] = FCV4(Preq,TOperation,Pint4,i_den)
PFC = 1000; %Initialize loop
flowO2FC =25;
while Preq > PFC
    if Preq/PFC > 1
    flowO2FC = flowO2FC + (Preq /PFC);
    else 
       flowO2FC = flowO2FC + 1;
    end

X_O2 = 1; %Assuming Fresh Oxygen supply is always at 100% concentration
T = TOperation; %Target operating temperature from OTM permeate stream
F = 96485.33; %Faraday's Constant in Coulombs/mol
Ru = 8.314; %Universal gas constant, J/mol.K  
%i_den = 1.5; %Average current density Amp/cm^2
W = 10; %Width of cell in cm
L = 10; %Length of cell in cm
n = 10; %Number of nodes at which the voltage change is calculated
ASR = 0.2; %Area specific resistance Ohm*cm^2
%i_den = 1.5
i_total = 4*flowO2FC*F; %Total current produced at 100% O2 utilization
A_total = i_total/i_den; %cm^2
Cells = A_total/(L*W); % Number of cells based on average current density and size of cells
i_Cell = i_total/Cells; %Total amount of current per cell
h = enthalpy2(T); %Enthalpy of products and reactants at operating temperature
s = entropy2(T);  %Entropy of products and reactants at operating temperature
E0 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F); %Reference Voltage
i = (ones(n,1)).*(i_Cell./(L.*W)); %Initial current distribution per cell
J_int = zeros(n,1); %Initializing the matrix for Current distribution
error = 1;
error2 = 5;
error3 = 50; 
PrFC = 10; %Operating Pressure Ratio of FC in kPa/100
%hrxn1 = -(2*refpropm('H','T',298,'P',PrFC*100,'water')*0.018 - 2*refpropm('H','T',298,'P',PrFC*100,'HYDROGEN')*0.002 -refpropm('H','T',298,'P',PrFC*100,'OXYGEN')*0.016)/1000; %Heat of reaction kJ/mol
hrxn = enthalpy2(298); %Enthalpy to calculate heat of reaction 
hrxnmol = -(2*hrxn.H2O-2*hrxn.H2-hrxn.O2); %2H2 + O2 -->  2H2O, Heat Released in kW/kmol
SteamRatio = 0.01; %Ratio of Steam in fuel intake for optimal cell operation and thermal balancing
reactH2FC = 2*flowO2FC; %Flow of H2 Needed to React 100% O2 Flow mol/s
H2Oproduced = 2*flowO2FC; %Water/Steam Produced by Reaction, mol/s
mH2Oproduced = 0.018*H2Oproduced %Water present in exhaust, kg/s
V = .5;


 flowH2FC = reactH2FC/0.1; %Initial guess for the required Fuel flow for thermal balancing, mol/s
 flowSteam = SteamRatio*flowH2FC; %Molar flow of Steam Required in Fuel Intake, mol/s
 flowTotalInlet = flowH2FC + flowSteam; %Total molar flow of fuel stream at inlet
 CpH2mol = refpropm('C','T',TOperation,'P',Pint4,'HYDROGEN')*0.002/1000; %Specific Heat of H2 as a Function of Temperature, kJ/mol*K
 CpH2Omol = refpropm('C','T',TOperation,'P',Pint4,'WATER')*0.018/1000; %Specific Heat of H2O as a Function of Temperature, kJ/mol*K

 
    
     
     for k=1:1:n
        J_int(k) = sum(i(1:k))*(W*L/n); %integral of current density as a function of length, total current thus far
        X_H2(k,1) = (flowH2FC - J_int(k)/(2*F/Cells))/(flowTotalInlet); %Concentration of H2 as a function of position and steam concentration
        X_H2O(k,1) = 1-X_H2(k,1); %Concentration of H20 product as H2 is consumed 
        E(k,1) = E0 + Ru*T/(2*F)*log(X_H2(k,1)/X_H2O(k,1)*(Pint4/100)^0.5); %Nernst Potential as a function of product and reactant concentrations
    end
    error2 = 1;
    count=0;
      if error == 1
         V =mean(E-i*ASR);
      end
      Vold = V;
    while abs(mean(error2))>(i_total*1e-6) || count < 2
        i = max(0,(E-V)/ASR); %new current distribution
        error2 = (sum(i)/n*L*W) - i_Cell; %error in total current
        V = V + .5*(error2/(L*W)*ASR); %New voltage distribution
        count = count + 1;
    end  
  
    PFC = mean(V)*i_Cell/1000*Cells %Electric Power produced by FC, kW
    EffFC = mean(V)/mean(E); %Voltage Efficiency of FC
    X_H2out = X_H2(n,1);
    X_H2Oout = X_H2O(n,1);
    X_H2in = X_H2(1,1);
    X_H2Oin = X_H2O(1,1);
     error = Vold-V;
    
    Utilization = (reactH2FC)/(flowH2FC);
    Qgen = (i_total/(4*F)*hrxnmol - V*i_total)/1000; %heat release from electrochemistry in kW
    Qout = (1-Utilization)*flowH2FC*(CpH2mol*50) + (flowSteam)*CpH2Omol*50; %Heat removed from stack by unused H2 
    error3 = (Qout-Qgen)/Qgen;
    
    %flowH2FC = flowH2FC/(Qout/Qgen);
    flowSteam = SteamRatio*flowH2FC;
    flowTotalInlet = flowH2FC + flowSteam;  
   

mFlowH2in = flowH2FC*0.002; %Total mass flow of H2 into FC, kg/s
mFlowH2recirc = (1-Utilization)*mFlowH2in; %Amount of Hydrogen Recirculated
mFlowH2new = Utilization*mFlowH2in; % Amount of new Hydrogen Added at each pass
mFlowH2Oin = SteamRatio*flowH2FC*0.018; %Total mass flow H2O into FC, kg/s
mFlowO2 = flowO2FC*0.032; %Total mass flow O2 into FC, kg/s
mFlowH2out = (1-Utilization)*mFlowH2in; %Total mass flow of H2 in Exhaust
mFlowH2Oout = (mH2Oproduced + mFlowH2Oin); %Total mass flow of H2O in Exhaust  
mFlowInFuelTotal = mFlowH2in + mFlowH2Oin; %Total fuel inlet mass flow
mFlowOutFuelTotal = mFlowH2out + mFlowH2Oout; %Total exhaust mass flow
mX_H2in = mFlowH2in/(mFlowH2in + mFlowH2Oin); %Mass Fraction of Hydrogen in Fuel
mX_H2Oin = mFlowH2Oin/(mFlowH2in + mFlowH2Oin); %Mass Fraction of H2O in fuel
mX_H2Oout = (mFlowH2Oout/(mFlowH2out + mFlowH2Oout)); %Mass Fraction H2O in exhaust
mX_H2out = (mFlowH2out/(mFlowH2out + mFlowH2Oout)); % Mass Fraction of Hydrogen in Exhaust
CpH2 = refpropm('C','T',TOperation,'P',Pint4,'HYDROGEN')/1000; %Specific Heat of H2 as a Function of Temperature, kJ/kg*K
CpH2Og = refpropm('C','T',TOperation,'P',Pint4,'WATER')/1000; %Specific Heat of H2O as a Function of Temperature, kJ/kg*K
TE1 = TOperation + (Qgen)/((mFlowH2Oout)*CpH2Og + (1-Utilization)*mFlowH2in*CpH2); % Temperature of FC Exhaust
HE1 = refpropm('H','T',TE1,'P',Pint4,'HYDROGEN','WATER',[mX_H2out,mX_H2Oout])/1000; %Enthalpy of H2 leaving FC, kJ/kg
mX_H2E4 = (mFlowH2in/(mFlowH2in + mFlowH2Oout)); %Mass Fraction of Hydrogen before condensation
mX_H2OE4 = mFlowH2Oout/(mFlowH2in + mFlowH2Oout); %Mass Fraction of Water before condensation
TE2 = 1030;
TF5 = 1023;
PF5 = Pint4;
HF5 = refpropm('H','T',TF5,'P',PF5,'WATER','HYDROGEN',[mX_H2Oin,mX_H2in])/1000; %Temperature at fuel inlet
TF3= 450;
PF3 = Pint4 - 200; %Pressure of recirculated fuel after moving through sofc stack and heat exchanger
PF4 = Pint4;
%TF4 = TOperation;
TF2 = 77;
HF2 = refpropm('H','T',TF2,'P',Pint4,'HYDROGEN')/1000;
[WBlower,TF4,HF4] = Blower(mFlowH2in,Pint4); 
HF3 = refpropm('H','T',TF3,'P',PF3,'HYDROGEN')/1000; 
HF4 = refpropm('H','T',TF4,'P',PF4,'WATER','HYDROGEN',[mX_H2Oin,mX_H2in])/1000;
QE2_3 = mFlowH2in*(HF5 - HF4);
TE3 = TE2  - (QE2_3)/((mFlowH2Oout)*CpH2Og + mFlowH2out*CpH2);
HE3 = refpropm('H','T',TE3,'P',Pint4,'HYDROGEN','WATER',[mX_H2out,mX_H2Oout])/1000;
HE4 = (mFlowOutFuelTotal/(mFlowOutFuelTotal + mFlowH2new))*HE3 + (mFlowH2new/(mFlowOutFuelTotal + mFlowH2new))*HF2;
TE4 = refpropm('T','H',HE4*1000, 'P',Pint4,'HYDROGEN','WATER',[mX_H2E4,mX_H2OE4]);%Temperature of exhaust after mixing with new Fuel
HE2 = refpropm('H','T',TE2,'P',Pint4,'HYDROGEN','WATER',[mX_H2out,mX_H2Oout])/1000; 
QE1_2 = mFlowOutFuelTotal*(HE1 - HE2); 
TE5 = 450;
HE5 = refpropm('H','T',TE5,'P',Pint4,'HYDROGEN','WATER',[mX_H2out,mX_H2Oout])/1000;
QE4_F3 = mFlowOutFuelTotal*(HE4 - HE5);

end
end