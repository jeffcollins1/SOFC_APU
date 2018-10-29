function [Outlet,Prop,A2,A3,A4,A5,A6,IC] = compressor_2(Inlet,P_out,eff,FC,C1)
[n1,n2] = size(Inlet.T); 
%Solve for required temperature into heat pipe for air stream to remove all
%excess heat from fuel cell
for z1 = 1:n1
    for z2 = 1:n2
if FC.Qremove(z1,z2) > 0
    A6.O2(z1,z2) = Inlet.O2(z1,z2);
    A6.N2(z1,z2) = Inlet.N2(z1,z2);
    A6.T(z1,z2) = 998; % Desired inlet temperature for SOFC
    targetout.O2 = A6.O2(z1,z2);
    targetout.N2 = A6.N2(z1,z2);
    targetout.T = A6.T(z1,z2);
    targetinHP.O2 = targetout.O2;
    targetinHP.N2 = targetout.N2;
    targetinHP.T = 500; %Initial guess 
    H_targetin = enthalpy(targetout) - FC.Qremove(z1,z2);
    targetinHP.T = find_T(targetinHP,H_targetin); %Required temperature into heat pipe to remove all excess heat from FC and achieve 998K
    if targetinHP.T > Inlet.T(z1,z2) + 50 % Check to see if target temperature out of second compressor stage into heat pipe is achievable using ambient cooling 
    
    A5.T(z1,z2) = targetinHP.T;
    end
    end
    
 PR1 = sqrt(P_out(z1,z2)./Inlet.P(z1,z2)); %Pressure ratio for one stage
 PR2 = PR1^2;% Pressure ratio for two stages
 
 %Solve for temperature out of stages 1 and 2
 
Ru = 8.314; %Universal Gas constant, kJ/(kmol*K)
incomp1.O2 = Inlet.O2(z1,z2);
incomp1.N2 = Inlet.N2(z1,z2);
incomp1.T = Inlet.T(z1,z2);
outcomp1.O2 = Inlet.O2(z1,z2);
outcomp1.N2 = Inlet.N2(z1,z2);
incomp2.O2 = Inlet.O2(z1,z2);
incomp2.N2 = Inlet.N2(z1,z2);
incomp2.T = 500; %initial guess
outcomp2.O2 = Inlet.O2(z1,z2);
outcomp2.N2 = Inlet.N2(z1,z2);
[~,H_in1] = enthalpy(incomp1); %Enthalpy into first compressor inlet
S_in1 = entropy(incomp1); %
cp = SpecHeat(incomp2);
gamma = cp./(cp - Ru); 
outcomp1.T= (incomp1.T*(PR1).^((gamma-1)./gamma)); 
[~,H_ideal1] = enthalpy(outcomp1); %Ideal enthlapy out of first compressor stage without intercooling
H_out1_true = H_in1 + 1./eff(z1,z2).*(H_ideal1 - H_in1); %True enthalpy out of first compressor stage without intercooling
W_c1(z1,z2) = H_out1_true - H_in1;
outcomp1.T = find_T(outcomp1,H_out1_true); %Actual temperature out of comp 1 without intercooling;

outcomp2.T= (outcomp1.T*(PR1).^((gamma-1)./gamma)); %Ideal temperature out of comp2 without intercooling
[~,H_ideal2] = enthalpy(outcomp2); %Ideal enthlapy out of second compressor stage without intercooling
H_out2_true = H_in1 + 1./eff(z1,z2).*(H_ideal2 - H_in1); %True enthalpy out of second compressor stage without intercooling
outcomp2.T = find_T(outcomp2, H_out2_true); %Actual temperature out of second compressor stage without intercooling
W_c2(z1,z2) = H_out2_true -H_out1_true; 
%Find target temperature out of comp 1 and into comp 2 
H_in2_target_difference = H_out2_true - H_targetin;  %Difference between actual outlet enthalpy from second compressor stage without intercooling and desired outlet enthalpy;
H_out1_target = H_out1_true - H_in2_target_difference;
temp_in2 = find_T(outcomp1,H_out1_target); %required temperature into compressor 2
%Check to see if required inlet temperature to comp 2 is achievable, 
if temp_in2 < Inlet.T(z1,z2) + 10
    incomp2.T = Inlet.T(z1,z2) + 10;
    outcomp2.T= (incomp2.T*(PR1).^((gamma-1)./gamma)); %Ideal temperature out of comp2 from minimum inlet temperature
[~,H_ideal_out2] = enthalpy(outcomp2); %Ideal enthlapy out of second compressor stage with maximum interooling
H_out2_true = H_in1 + 1./eff(z1,z2).*(H_ideal_out2 - H_in1); %True enthalpy out of second compressor stage with maximum intercooling
H_in2 = enthalpy(incomp2);
W_c2(z1,z2) = H_out2_true - H_in2;
outcomp2.T = find_T(outcomp2, H_out2_true); %Actual temperature out of comp 2 with max intercooling
Q_out_IC = enthalpy(incomp2) - enthalpy(outcomp1); % heat rejected by intercooling stage
Icoolin.O2(z1,z2) = Inlet.O2(z1,z2); %Initial guess for mass flow of intercooling stream
Icoolin.N2(z1,z2) = Inlet.N2(z1,z2);
Icoolin.T(z1,z2) = Inlet.T(z1,z2); 
Icoolout.O2(z1,z2) = Inlet.O2(z1,z2); %Initial guess for mass flow of intercooling stream
Icoolout.N2(z1,z2) = Inlet.N2(z1,z2);
Icoolout.T(z1,z2) = outcomp1.T - 10; %Maximum temperature achieved by intercooling air stream
Q_out_IC_max(z1,z2) = enthalpy(Icoolout) - enthlapy(Icoolin); %Maximum heat rejection betweeen stages 1 and 2
Q_out_after_comp2(z1,z2) = H_targetin - Q_out_IC_max; %Remaining heat to be rejected before compressor 2 outlet reaches heat pipe
HPin = outcomp2;
H_HPin = enthalpy(outcomp2) - Q_out_after_comp2; %
HPin.T = find_T(outcomp2,H_HPin); 
Icoolin2 = Icoolin; %Inital guess for additional cooling air between comp2 and heat pipe; 
Icoolin2.T = Icoolin.T;
Icoolout2 = Icoolin2;
Icoolout2.T = HPin.T - 10; %Maximum temperature achieved by second cooling stream
Q_Icool2 = enthalpy(Icoolout2) - enthalpy(Icoolin2); %Actual heat aborbed by initial mass flow guess
Icoolin2.N2 = Q_out_after_comp2/Q_Icool2*Icoolin.N2; %Update mass flow required 
Icoolin2.O2 = Q_out_after_comp2/Q_Icool2*Icoolin.O2;
Prop.work(z1,z2) = W_c1(z1,z2) + W_c2(z1,z2);
IC.out2(z1,z2) = Icoolout2;
IC.in2(z1,z2) = Icoolin2;
else
    prop.work(z1,z2) = C1.work(z1,z2);
A6.O2(z1,z2) = Inlet.O2(z1,z2); %Assign A6 if not assigned in earlier if statement
    A6.N2(z1,z2) = Inlet.N2(z1,z2);
    A6.T(z1,z2) = 998; % Desired inlet temperature for SOFC
    incomp2.T = temp_in2;
    H_in2 = enthalpy(incomp2);
    outcomp2.T= (incomp2.T*(PR1).^((gamma-1)./gamma)); %Ideal temperature out of comp2 from minimum inlet temperature
[~,H_ideal_out2] = enthalpy(outcomp2); %Ideal enthlapy out of second compressor stage with maximum interooling
H_out2_true = H_in1 + 1./eff(z1,z2).*(H_ideal_out2 - H_in1); %True enthalpy out of second compressor stage with required intercooling
outcomp2.T = find_T(outcomp2, H_out2_true); %Actual temperature out of comp 2 with required intercooling
Q_out_IC(z1,z2) = enthalpy(incomp2) - enthalpy(outcomp1); % heat rejected by intercooling stage
Icoolin.O2(z1,z2) = Inlet.O2(z1,z2); %Initial guess for mass flow of intercooling stream
Icoolin.N2(z1,z2) = Inlet.N2(z1,z2);
Icoolin.T(z1,z2) = Inlet.T(z1,z2);
H_Icoolin = enthalpy(Icoolin);
Icoolout.O2(z1,z2) = Inlet.O2(z1,z2); %Initial guess for mass flow of intercooling stream
Icoolout.N2(z1,z2) = Inlet.N2(z1,z2);
Icoolout.T(z1,z2) = incomp2.T - 10; %Actual temperature out at 
Q_Icoolout = enthalpy(Icoolout) - enthalpy(Icoolin);
Icoolout.O2(z1,z2) = (Q_out_IC(z1,z2)/Q_Icoolout)*Inlet.O2(z1,z2); %Adjusted mass flow of intercooling stream
Icoolout.N2(z1,z2) = (Q_out_IC/Q_Icoolout)*Inlet.N2(z1,z2);

A2.N2(z1,z2) = outcomp1.N2;
A2.O2(z1,z2) = outcomp1.O2
A2.T(z1,z2) = outcomp1.T;
A3.N2(z1,z2) = incomp2.N2;
A3.O2(z1,z2) = incomp2.O2;
A3.T(z1,z2) = incomp2.T;
A4.O2(z1,z2)= outcomp2.O2;
A4.N2(z1,z2)= outcomp2.N2;
A4.T(z1,z2) = outcomp2.T;
A5.O2(z1,z2)= targetinHP.O2;
A5.N2(z1,z2)= targetinHP.N2;
A5.T(z1,z2) = targetinHP.T; 
IC.out(z1,z2) = Icoolout;
IC.in(z1,z2) = Icoolin;
IC.out2(z1,z2) = 0;
IC.in2(z1,z2) = 0;

    end
    end
end

Prop.pressure_ratio = P_out./Inlet.P;
Prop.eff = eff;
Prop.mass_flow = mass_flow(Inlet);
end

