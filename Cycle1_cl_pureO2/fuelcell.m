function [FC,Exhaust] = fuelcell(options,Cathode)
F = 96485.33; %Faraday's Constant in Coulombs/mol
Ru = 8.314; %Universal gas constant, J/mol.K  
T = options.T_fc;
P = options.P_fc;
%i_den = 1.5; %Average current density Amp/cm^2
W = 10; %Width of cell in cm
L = 10; %Length of cell in cm
n = 10; %Number of nodes at which the voltage change is calculated
i_total = 4000*Cathode.O2*F; %Total current produced at 100% O2 utilization
FC.i_den = i_total/(options.SOFC_area*10000); %A/cm^2

FC.Cells = options.SOFC_area*10000/(L*W); % Number of cells based on average current density and size of cells
i_Cell = i_total/FC.Cells; %Total amount of current per cell
G_H2O = refpropm('h','T',T,'P',P,'WATER')/1000*18-241826.4-refpropm('s','T',T,'P',P,'WATER')/1000*18.*T;
G_H2 = refpropm('h','T',T,'P',P,'HYDROGEN')/1000*2-refpropm('s','T',T,'P',P,'HYDROGEN')/1000*2.*T;
G_O2 = refpropm('h','T',T,'P',P,'OXYGEN')/1000*32-refpropm('s','T',T,'P',P,'OXYGEN')/1000*32.*T;
E0 = -(G_H2O - G_H2 - .5*G_O2)/(2*F);%Reference Voltage

i = (ones(n,1)).*(i_Cell./(L.*W)); %Initial current distribution per cell
J_int = zeros(n,1); %Initializing the matrix for Current distribution


hrxnmol = -(2*refpropm('h','T',T,'P',P,'WATER')-2*refpropm('h','T',T,'P',P,'HYDROGEN')-refpropm('h','T',T,'P',P,'OXYGEN')); %2H2 + O2 -->  2H2O, Heat Released in kW/kmol
FC.H2_used = 2*Cathode.O2; %Flow of H2 Needed to React 100% O2 Flow mol/s
FC.H2_supply = FC.H2_used/options.spu;
FC.H2O_supply = FC.H2_supply*options.steamratio/(1-options.steamratio);
n_in = FC.H2_supply + FC.H2O_supply;  

Anode.T = T;
Anode.P = options.P_fc;
Anode.H2 = FC.H2_supply;
Anode.H2O = FC.H2O_supply; 

Exhaust.T = Anode.T + options.dT_fc;
Exhaust.P = Anode.P;
Exhaust.H2 = Anode.H2 - FC.H2_used;
Exhaust.H2O = Anode.H2O + 2*Cathode.O2; %Water/Steam Produced by Reaction, mol/s

error_V = 1;
while abs(error_V)>1e-4
    X_H2 = zeros(n,1);   
    X_H2O = zeros(n,1);
    E = zeros(n,1);
    for k=1:1:n
        J_int(k) = sum(i(1:k))*(W*L/n); %integral of current density as a function of length, total current thus far
        X_H2(k,1) = (FC.H2_supply - J_int(k)/(2000*F/FC.Cells))/(n_in); %Concentration of H2 as a function of position and steam concentration
        X_H2O(k,1) = 1-X_H2(k,1); %Concentration of H20 product as H2 is consumed 
        E(k,1) = E0 + Ru*T/(2000*F)*log(X_H2(k,1)/X_H2O(k,1)*(options.P_fc/100)^0.5); %Nernst Potential as a function of product and reactant concentrations
    end
    error = 1;
    count=0;
    if error_V == 1
        V = mean(E-i*options.asr);
        Vold = V;
    end
    while abs(mean(error))>(i_total*1e-6) %|| count < 2
        i = max(0,(E-V)/options.asr); %new current distribution
        error = (sum(i)/n*L*W) - i_Cell; %error in total current
        V = V + .5*(error/(L*W)*options.asr); %New voltage distribution
        count = count + 1;
    end  
    error_V = Vold - V;
    Vold = V;
end
FC.V = mean(V);
FC.Power = FC.V*i_Cell/1000*FC.Cells; %Electric Power produced by FC, kW

FC.Qremove = (i_total/(4000*F)*hrxnmol - FC.V*i_total)/1000 - (property(Exhaust,'h','kJ')-property(Anode,'h','kJ')); %heat removed by some means to maintain temperature in kW

end