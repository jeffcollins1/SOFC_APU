function [Ppropmax,thrustN,time_to,Vp,plot_V,plot_P,plot_T,plot_D] = thrust_to(options)
Tamb = 298*ones(10,10); %Ambient Temperature Kelvin;
Pamb  = 100*ones(10,10); %Ambient Pressure as a function of altitude
thrustNmax =  options.thrust_to(1,10); %Maximum thrust generated per engine
thrustN = options.thrust_to; 
A = (0.25*pi()*4^2); % Cross sectional area of prop
R = 8.314/(28.8); % gas constant for air, kJ/kgK 
gamma1 = (refproparray('C','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79])./refproparray('O','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]));%property(A1,'C','kJ/(kmol K)')./property(Inlet,'O','kJ/(kmol K)');
gamma = gamma1(1:10,1)*ones(1,10); 
density1 =  refproparray('D','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]); %property(A1,'D','kJ/(kmol K)'); 
density = density1(1:10,1)*ones(1,10);
speedsound1 = refproparray('A','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]);
speedsound = speedsound1(1:10,1)*ones(1,10); 
d = (ones(10,1)*linspace(0,3200,10));
Vcraft_to = (ones(10,1)*linspace(10,80,10));
Vcraft_to_max = 85; %Takeoff speed, m/s
Ve = sqrt(2*thrustN./(density.*A) + Vcraft_to.^2); %Velocity of exit stream 
Vp = 0.5*Ve + 0.5*Vcraft_to; %Velcocity at propeller, average of free stream velocity and jet exit velocity

Mair = Vp./speedsound; 

Pp = Pamb.*(1 + 0.5.*(gamma -1).*Mair.^2).^(gamma./(gamma-1)); %Isentropic stagnation pressure into propeller
Tp = Tamb.*(Pp./Pamb).^((gamma-1)./gamma); %Isentropic Temperature into propeller
if Mair <= 0.3
    mflow = density.*A.*Vp;
else
mflow = density.*((1 + 0.5.*(gamma -1).*Mair.^2).^(1./(gamma-1))).*A.*Vp;
end
%Pp = Pamb.*(1 + 0.5.*(gamma -1).*Mair.^2).^(gamma./(gamma-1)); %Isentropic stagnation pressure into propeller
Pe = Pp.*(1 + 0.5.*(gamma -1).*Mair.^2).^(gamma./(gamma-1));%Isentropic stagnation pressure leaving properller
Pprop = mflow.*(Pe-Pp)./density; 
Ppropmax = Pprop(10,10);
% Pe = Pp.*(1 + 0.5.*(gamma -1).*Mair.^2).^(gamma./(gamma-1));
% Pprop = mflow.*(Pe-Pp)./density; 
% Ppropmax = max(max(Pprop)); 
 
 plot_V = Vcraft_to(1,1:10);
 plot_P = Pprop(1,10);
 plot_T =thrustN(1,1:10); 
 plot_D = d(1,1:10); 

% weight_to = 362870; %maximum takeoff weight, kg;
 time_to = 3340/40; %Total time for takeoff assuming constant acceleration
end
