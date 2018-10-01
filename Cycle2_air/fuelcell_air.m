function [FC,Exhaust,E5,A4] = fuelcell_air(options,intake)
F = 96485.33; %Faraday's Constant in Coulombs/mol
Ru = 8.314; %Universal gas constant, J/mol.K  
T = options.T_fc;
P = options.P_fc;
%i_den = 1.5; %Average current density Amp/cm^2
W = 9; %Width of active cell area in cm
L = 9; %Length of active cell area in cm
n = 10; %Number of nodes at which the voltage change is calculated
i_totalmax = 4000*intake.Y_O2*F; %Total current produced at 100% O2 utilization
i_totalmin = 0.25*i_totalmax; 
i_total = ones(n1,1)*linspace(i_totalmin,i_totalmax,10);
FC.i_den = i_total./(options.SOFC_area*10000); %A/cm^2

FC.Cells = options.SOFC_area*10000/(L*W); % Number of cells based on average current density and size of cells
i_Cell = i_total./FC.Cells; %Total amount of current per cell
% G_H2O = refpropm('h','T',T,'P',100,'WATER')/1000*18-241826.4-refpropm('s','T',298,'P',100,'WATER')/1000*18
% G_H2 = refpropm('h','T',T,'P',100,'HYDROGEN')/1000*2-refpropm('s','T',T,'P',100,'HYDROGEN')/1000*2.*T;
% G_O2 = refpropm('h','T',T,'P',100,'OXYGEN')/1000*32-refpropm('s','T',T,'P',100,'OXYGEN')/1000*32.*T;
% E0 = -(G_H2O - G_H2 - .5*G_O2)/(2*F);%Reference Voltage
T_avg = T + .5*options.dT_fc;

if T_avg>=900 && T_avg < 1000
    FC.G = -198141 + (T_avg-900)*(-192652+198141)/(1000-900);
elseif T_avg >=1000 && T_avg<1100
    FC.G = -192652 + (T_avg-1000)*(-187100+192652)/(1100-1000);
end
FC.G0 = -228608; %Enthalpy of formation for water at 298K & 100 kPa, kJ/kmol;
E0 = -FC.G/(2*F);
i = (ones(n,1)).*(i_Cell./(L.*W)); %Initial current distribution per cell
J_int = zeros(n,1); %Initializing the matrix for Current distribution

h_H2O = ((refpropm('h','T',T,'P',P,'WATER') - refpropm('h','T',298,'P',100,'WATER'))*18/1000 -44010 -241847);% total enthalpy in kJ/kmol
h_H2 = ((refpropm('h','T',T,'P',P,'hydrogen') - refpropm('h','T',298,'P',100,'hydrogen'))*2/1000);% total enthalpy in kJ/mol
h_O2 = ((refpropm('h','T',T,'P',P,'oxygen') - refpropm('h','T',298,'P',100,'oxygen'))*32/1000);% total enthalpy in kJ/kmol
s_H2O = ((refpropm('s','T',T,'P',P,'WATER') - refpropm('s','T',298,'P',100,'WATER'))*18/1000 -44010 -241847);% total enthalpy in kJ/kmol
s_H2 = ((refpropm('s','T',T,'P',P,'hydrogen') - refpropm('s','T',298,'P',100,'hydrogen'))*2/1000);% total enthalpy in kJ/kmol
s_O2 = ((refpropm('s','T',T,'P',P,'oxygen') - refpropm('s','T',298,'P',100,'oxygen'))*32/1000);% total enthalpy in kJ/kmol
cp_air = 0.21*property('O2','C','kJ/kmol') + 0.79*property('N2','C','kJ/kmol'); % specific heat of air
FC.hrxnmol = -(h_H2O - h_H2 - .5*h_O2); %H2 + 0.5*O2 -->  H2O, Heat Released in kJ/kmol
FC.Xch_in = FC.hrxnmol- 298*(s_H2O + s_H2 + 0.5*s_O2); %Chemical exergy present in fuel stream, kJ/kmol (0.019*(h_H2O-298*s_H2O)- 0.2*(h_O2 - 298*s_O2)) -
% hrxnmol = 241847;
FC.H2_used = 2.*i_total./(4000*F); %Range of H2 flows needed to react range of current densities, kmol/s
FC.H2_supply = FC.H2_used./options.spu;
FC.H2O_supply = FC.H2_supply*options.steamratio/(1-options.steamratio);
FC.n_in = FC.H2_supply + FC.H2O_supply;  

Anode.T = T;
Anode.P = options.P_fc;
Anode.H2 = FC.H2_supply;
Anode.H2O = FC.H2O_supply; 
Anode.Y_H2 = FC.H2_supply/(FC.H2_supply+FC.H2O_supply);
Anode.Y_H2O = FC.H2O_supply/(FC.H2_supply+FC.H2O_supply)
Exhaust.T = Anode.T + options.dT_fc;
Exhaust.P = Anode.P;
Exhaust.H2 = Anode.H2 - FC.H2_used;
Exhaust.H2O = Anode.H2O + 2*intake.Y_O2; %Water/Steam Produced by Reaction, mol/s
Exhaust.Y_H2O = Exhaust.H2O/(Exhaust.H2 + Exhaust.H2O); %Molar composition of water in exhaust
Exhaust.Y_H2 = Exhaust.H2/(Exhaust.H2 + Exhaust.H2O); %Molar composition of hydrogen in exhaust
error_V = 1;
error_Q = 1;
E5 = Exhaust.H2O - FC.H2O_supply; 
A1 = Cathode.O2/0.21;
A2 = A1;
A2.T = intake.A2T
A3.T = T; 
while abs(error_Q) > 1e-4
   
while abs(error_V)>1e-4
    X_H2 = zeros(n,1);   
    X_H2O = zeros(n,1);
    X_O2 = zeros(n,1);
    E = zeros(n,1);
    for k=1:1:n
        J_int(k) = sum(i(1:k))*(W*L/n); %integral of current density as a function of length, total current thus far
        X_H2(k,1) = (FC.H2_supply - J_int(k)/(2000*F/FC.Cells))/(FC.n_in); %Concentration of H2 as a function of position and steam concentration
        X_H2O(k,1) = 1-X_H2(k,1); %Concentration of H20 product as H2 is consumed 
        X_O2(k,1) = 1-(Cathode.Y_O2 - 0.5*J_int(k)/(2000*F/FC.Cells))/Cathode.Y_O2; %Concentration of oxygen as a function of position
        E(k,1) = E0 + Ru*T/(2*F)*log(X_H2(k,1)*X_O2(k,1)/X_H2O(k,1)*(P/100)^0.5); %Nernst Potential as a function of product and reactant concentrations
    end
    error = 1;
    count=0;
    if error_V == 1
        V = mean(E-i*options.asr)+.05;
        Vold = V;
    end
    while abs(mean(error))>(i_Cell*1e-6) %|| count < 2
        i = max(0,(E-V)/options.asr); %new current distribution
        error = (sum(i)/n*L*W) - i_Cell; %error in total current
        V = V + .5*(error/(L*W)*options.asr); %New voltage distribution
        count = count + 1;
    end  
    error_V = Vold - V;
    Vold = V;
end
FC.Cathode_outlet = options.airflow*0.79 + options.airflow*X_O2(n,1); % flow and compostion of cathode exhaust
FC.V = mean(V);
FC.Power = FC.V*i_Cell/1000*FC.Cells; %Electric Power produced by FC, kW
FC.Qgen = i_total/(2000*F)*FC.hrxnmol - FC.V*i_total/1000; %Heat generated by fuel cell, kW
A4.N2 = A3.N2; %Total mass flow
A4.O2 = itotal./(4000*F); 
A4.T = options.T_fc; 
A4.P = options.P_fc - 20*ones(10,10); 
FC.Qremove = FC.Qgen - (property(Exhaust,'h','kJ')-property(Anode,'h','kJ')-propery(intake)); %heat removed by some means to maintain temperature in kW
error_Q = FC.Qremove - (property(A4,'h','kJ')-property(A3,'h','kJ'); %heat removed by excess air flow

end