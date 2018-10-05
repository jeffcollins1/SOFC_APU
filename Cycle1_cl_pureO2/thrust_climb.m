function [Ppropmax_climb,Ppropmax_cruise,time_climb,thrustact,thrusterr_climb,thrusterr_cruise,Vp,plot_V_climb,plot_P_climb,plot_T_climb,plot_D] = thrust_climb(OD,param_od)
Tamb = ones(10,1)*OD.Tamb(1,1:10); %Ambient Temperature Kelvin;
Pamb  = ones(10,1)*OD.Pamb(1,1:10); %Ambient Pressure as a function of altitude
Tamb_cruise = ones(10,1)*OD.Tamb_cruise(1,1:10); %Ambient Temperature at cruise altitude Kelvin
Pamb_cruise  = ones(10,1)*OD.Pamb_cruise(1,1:10); %Ambient Pressure at cruise altitude Kelvin
thrustNmax = 4.44822*12726; 
V_o = 85; %takeoff speed, 747
V_f = 200; %final velocity at the end of climb
thrustN_climb = 4.44822*OD.thrust_climb; 
thrustN_cruise = 4.44822*OD.thrust_cruise; 
A = (0.25*pi()*6^2); % Cross sectional area of prop
R = 8.314/(28.8); % gas constant for air, kJ/kgK 
gamma1 = (refproparray('C','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79])./refproparray('O','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]));
gamma2 = (refproparray('C','T',Tamb_cruise,'P',Pamb_cruise,'OXYGEN','NITROGEN',[0.21,0.79])./refproparray('O','T',Tamb_cruise,'P',Pamb_cruise,'OXYGEN','NITROGEN',[0.21,0.79]));
gamma_climb = gamma1(1:10,1)*ones(1,10); 
gamma_cruise = gamma2(1:10,1)*ones(1,10); 
density1 =  refproparray('D','T',Tamb(1,1:10),'P',Pamb(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]); %property(A1,'D','kJ/(kmol K)'); 
density2 =  refproparray('D','T',Tamb_cruise(1,1:10),'P',Pamb_cruise(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]); %property(A1,'D','kJ/(kmol K)'); 
density_climb = (ones(10,1)*density1(1,1:10));
density_cruise = (ones(10,1)*density2(1,1:10));
speedsound1 = refproparray('A','T',Tamb(1,1:10),'P',Pamb(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
speedsound2 = refproparray('A','T',Tamb_cruise(1,1:10),'P',Pamb_cruise(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
speedsound_climb = (ones(10,1)*speedsound1(1,1:10)); 
speedsound_cruise = (ones(10,1)*speedsound2(1,1:10)); 
%climb phase
height = OD.height;
climb_angle = 15*pi/180;
climb_angle_comp = 75*pi/180;
d_horizontal = height(10,10)*cos(climb_angle); %total horizontal distance travelled based on climb angle and altitude
d_horizontalVec = (ones(10,1)*linspace(0,d_horizontal,10));
d_total = sqrt(d_horizontalVec.^2 + height.^2);
V_craft_climb = OD.velocity_climb;%Magnitude of velocity during climb
V_craft_vert = 7.5; %Max climb rate of craft, m/s
Ve = sqrt(2.*thrustN_climb./(density_climb.*A) + V_craft_climb.^2);
Vp = 0.5*Ve + 0.5*V_craft_climb; %Velcocity at propeller, average of free stream velocity and jet exit velocity
Mair = Vp./speedsound_climb; 
m = 362870*ones(10,1); %take off mass, kg
ag = 9.81*ones(10,1); % acceleration in m/s

  
Pp = Pamb.*(1 + 0.5.*(gamma_climb -1).*Mair.^2).^(gamma_climb./(gamma_climb-1)); %Isentropic stagnation pressure into propeller
Tp = Tamb.*(Pp./Pamb).^((gamma_climb-1)./gamma_climb); %Isentropic Temperature into propeller
mflow = A.*V_craft_climb.*density_climb.*(1 + 0.5.*(gamma_climb -1).*Mair.^2).^(1./(gamma_climb-1));
thrustact = mflow.*(Ve-V_craft_climb); 
Pe = Pp.*(1 + 0.5.*(gamma_climb -1).*Mair.^2).^(gamma_climb./(gamma_climb-1));
Pprop_climb = mflow.*(Pe-Pp)./density_climb; 
   

Ppropmax = max(max(Pprop_climb)); 

 thrusterr_climb = thrustN_climb - thrustact;
 min_err_climb = min(min(abs(thrusterr_climb)))
 [i,j] = find(abs(thrusterr_climb)==min_err_climb);
 Ppropmax_climb = Pprop_climb(i,j)*ones(10,1);
 plot_V_climb = V_craft_climb(1,1:10);
 plot_P_climb = Pprop_climb(1,1:10);
 plot_T_climb =thrustact(1,1:10); 
 plot_D = d_total(1,1:10); 
 time_climb = d_total(10,10)/140; %Total time for climb based on average velocity and total distance    
%Cruise 
V_craft_cruise = 0.8*speedsound_cruise;%Magnitude of velocity during climb
Ve_cruise = sqrt(2.*thrustN_cruise./(density_cruise.*A) + V_craft_cruise.^2);
Vp_cruise = 0.5*Ve_cruise + 0.5*V_craft_cruise; %Velcocity at propeller, average of free stream velocity and jet exit velocity
Mair_cruise = 0.85; 

Pp_cruise = Pamb_cruise.*(1 + 0.5.*(gamma_cruise -1).*Mair_cruise.^2).^(gamma_cruise./(gamma_cruise-1)); %Isentropic stagnation pressure into propeller
Tp_cruise = Tamb_cruise.*(Pp_cruise./Pamb_cruise).^((gamma_cruise-1)./gamma_cruise); %Isentropic Temperature into propeller
mflow_cruise = A.*V_craft_cruise.*density_cruise.*(1 + 0.5.*(gamma_cruise -1).*Mair_cruise.^2).^(1./(gamma_cruise-1));
thrustact_cruise = mflow_cruise.*(Ve_cruise-V_craft_cruise); 
Pe_cruise = Pp_cruise.*(1 + 0.5.*(gamma_cruise -1).*Mair_cruise.^2).^(gamma_cruise./(gamma_cruise-1));
Pprop_cruise = mflow_cruise.*(Pe_cruise-Pp_cruise)./density_cruise; 
thrusterr_cruise = thrustN_cruise - thrustact_cruise;
 min_err_cruise = min(min(abs(thrusterr_cruise)))
 [k,l] = find(abs(thrusterr_cruise)==min_err_cruise);
 Ppropmax_cruise = Pprop_cruise(k,l)*ones(10,1);
thrusterr_cruise = thrustN_cruise - thrustact; 

end