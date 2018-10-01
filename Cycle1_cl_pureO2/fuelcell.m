function [FC,Exhaust] = fuelcell(options,Cathode)
F = 96485.33; %Faraday's Constant in Coulombs/mol
Ru = 8.314; %Universal gas constant, J/mol.K  
T = options.T_fc + 0.5*options.dT_fc;
P = options.P_fc(1:10,1:10);
[n1,n2] = size(T);
%i_den = 1.5; %Average current density Amp/cm^2
W = 10; %Width of cell in cm
L = 10; %Length of cell in cm
n = 10; %Number of nodes at which the voltage change is calculated
FC.i_total = 4000.*Cathode.O2.*F; %Total current produced at 100% O2 utilization

FC.i_den = zeros(n1,n2);
for x= 1:1:n1
    for y = 1:1:n1
        FC.i_den(x,y) = FC.i_total(x,y)/(options.SOFC_area(x,y)*10000); %A/cm^2
    end
end

FC.Cells = options.SOFC_area.*10000/(81); % Number of cells based on total active surface area of 81 cm^2 per cell
FC.i_Cell = 81*FC.i_den;%FC.i_total./FC.Cells; %Total amount of current per cell
% G_H2O = refpropm('h','T',T,'P',100,'WATER')/1000*18-241826.4-refpropm('s','T',298,'P',100,'WATER')/1000*18
% G_H2 = refpropm('h','T',T,'P',100,'HYDROGEN')/1000*2-refpropm('s','T',T,'P',100,'HYDROGEN')/1000*2.*T;
% G_O2 = refpropm('h','T',T,'P',100,'OXYGEN')/1000*32-refpropm('s','T',T,'P',100,'OXYGEN')/1000*32.*T;
% E0 = -(G_H2O - G_H2 - .5*G_O2)/(2*F);%Reference Voltage
T_avg = T + .5*options.dT_fc;
nine_to_hundred = (T_avg>=900 & T_avg < 1000);
FC.G = -198141 + (T_avg-900)*(-192652+198141)/(1000-900);
FC.G(~nine_to_hundred) = -192652 + (T_avg(~nine_to_hundred)-1000)*(-187100+192652)/(1100-1000);
FC.G0 = -228608; %Enthalpy of formation for water at 298K & 100 kPa, kJ/kmol;
E0 = -FC.G/(2*F);
i = (ones(n1,n2,n)).*(FC.i_Cell./(L.*W)); %Initial current distribution per cell
J_int = zeros(n1,n2,n); %Initializing the matrix for Current distribution
h_H2O = (refproparray('h','T',T_avg,'P',P,'WATER')*18/1000 -45839 -241847);% total enthalpy in kJ/kmol
h_H2 = ((refproparray('h','T',T_avg,'P',P,'hydrogen') - refpropm('h','T',298,'P',100,'hydrogen'))*2/1000);% total enthalpy in kJ/mol
h_O2 = ((refproparray('h','T',T_avg,'P',P,'oxygen') - refpropm('h','T',298,'P',100,'oxygen'))*32/1000);% total enthalpy in kJ/kmol
% s_H2O = ((refpropm('s','T',T,'P',P,'WATER') - refpropm('s','T',298,'P',100,'WATER'))*18/1000 -44010 -241847);% total enthalpy in kJ/kmol
% s_H2 = ((refpropm('s','T',T,'P',P,'hydrogen') - refpropm('s','T',298,'P',100,'hydrogen'))*2/1000);% total enthalpy in kJ/kmol
% s_O2 = ((refpropm('s','T',T,'P',P,'oxygen') - refpropm('s','T',298,'P',100,'oxygen'))*32/1000);% total enthalpy in kJ/kmol
hrxn = -(h_H2O - h_H2 - .5*h_O2); %H2 + 0.5*O2 -->  H2O, Heat Released in kJ/kmol
FC.hrxnmol = hrxn(1:10,1)*ones(1,10); 
% FC.Xch_in = FC.hrxnmol- 298*(s_H2O + s_H2 + 0.5*s_O2); %Chemical exergy present in fuel stream, kJ/kmol (0.019*(h_H2O-298*s_H2O)- 0.2*(h_O2 - 298*s_O2)) -
% hrxnmol = 241847;
FC.H2_used = 2*Cathode.O2; %Flow of H2 Needed to React 100% O2 Flow mol/s
FC.H2_supply = FC.H2_used./options.spu;
FC.H2O_supply = FC.H2_supply.*options.steamratio./(1-options.steamratio);
FC.n_in = FC.H2_supply + FC.H2O_supply;  

Anode.T = options.T_fc;
Anode.P = options.P_fc;
Anode.H2 = FC.H2_supply;
Anode.H2O = FC.H2O_supply; 
Anode.Y_H2 = FC.H2_supply./(FC.H2_supply+FC.H2O_supply);
Anode.Y_H2O = FC.H2O_supply./(FC.H2_supply+FC.H2O_supply);
Exhaust.T = Anode.T + options.dT_fc;
Exhaust.P = Anode.P;
Exhaust.H2 = Anode.H2 - FC.H2_used;
Exhaust.H2O = Anode.H2O + 2*Cathode.O2; %Water/Steam Produced by Reaction, mol/s
Exhaust.Y_H2O = Exhaust.H2O./(Exhaust.H2 + Exhaust.H2O); %Molar composition of water in exhaust
Exhaust.Y_H2 = Exhaust.H2./(Exhaust.H2 + Exhaust.H2O); %Molar composition of hydrogen in exhaust
error_V = 1;
V = zeros(n1,n2);
while max(max(abs(error_V)))>1e-4
    X_H2 = zeros(n1,n2,n);   
    X_H2O = zeros(n1,n2,n);
    E = zeros(n1,n2,n);
    for k=1:1:n
        J_int(:,:,k) = sum(i(:,:,1:k),3)*(W*L/n); %integral of current density as a function of length, total current thus far
        X_H2(:,:,k) = (FC.H2_supply - J_int(:,:,k)./(2000*F./FC.Cells))./(FC.n_in); %Concentration of H2 as a function of position and steam concentration
        X_H2O(:,:,k) = 1-X_H2(:,:,k); %Concentration of H20 product as H2 is consumed 
        E(:,:,k) = E0 + Ru*T_avg./(2*F).*log(X_H2(:,:,k)./(X_H2O(:,:,k).*(P/100).^0.5)); %Nernst Potential as a function of product and reactant concentrations
    end
    error = 1;
    count=0;
    if error_V == 1
        for j = 1:1:n1
            for k = 1:1:n2
                V(j,k) = mean(E(j,k,:)-i(j,k,:).*options.asr(j,k))+.05;
            end
        end
        Vold = V;
    end
    while all(abs(mean(error))>(FC.i_Cell*1e-6)) %|| count < 2
        for j = 1:1:n1
            for k = 1:1:n2
                i(j,k,:) = max(0,(E(j,k,:)-V(j,k))./options.asr(j,k)); %new current distribution
            end
        end
        error = (sum(i,3)/n*L*W) - FC.i_Cell; %error in total current
        V = V + .5*(error./(L*W).*options.asr); %New voltage distribution
        count = count + 1;
    end  
    error_V = Vold - V;
    Vold = V;
end
FC.V = V;
FC.Power = FC.V.*FC.i_total./1000%.*FC.Cells; %Electric Power produced by FC, kW
FC.Qgen =  FC.hrxnmol.*FC.i_total./(2000.*F) - FC.V.*FC.i_total./1000; %Heat generated by fuel cell, kW 
FC.Qremove = FC.Qgen - (property(Exhaust,'h','kJ')-property(Anode,'h','kJ')-property(Cathode,'h','kJ')); %heat removed by some means to maintain temperature in kW
FC.O2 = Cathode.O2;
end