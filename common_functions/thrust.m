function [Pprop,thrustact,thrusterr] = thrust_to(options)
Tamb = -0.0065.*options.height + 14.987 + 273.1; %Ambient Temperature Kelvin;
Pamb  = 107.*exp(-0.0001.*options.height)- 10; %Ambient Pressure as a function of altitude
thrustN = 4.44822*options.thrust; 
A = (0.25*pi()*6^2); % Cross sectional area of prop
R = 8.314/(28.8); % gas constant for air, kJ/kgK 
gamma1 = (refproparray('C','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79])./refproparray('O','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]));%property(A1,'C','kJ/(kmol K)')./property(Inlet,'O','kJ/(kmol K)');
gamma = gamma1(1:10,1)*ones(1,10); 
density1 =  refproparray('D','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]); %property(A1,'D','kJ/(kmol K)'); 
density = density1(1:10,1)*ones(1,10);
speedsound1 = refproparray('A','T',Tamb,'P',Pamb,'OXYGEN','NITROGEN',[0.21,0.79]);
speedsound = speedsound1(1:10,1)*ones(1,10); 
M = options.velocity./speedsound; 

P2 = Pamb.*(1 + 0.5.*(gamma -1).*M.^2).^(gamma./(gamma-1)); %Isentropic stagnation pressure into propeller
T2 = Tamb.*(P2./Pamb).^((gamma-1)./gamma); %Isentropic Temperature into propeller; 

if M <= 0.3
    mflow = density.*A.*options.velocity
else
mflow = density.*(1 + 0.5.*(gamma -1).*M^2).^(1./(gamma-1)).*A.*options.velocity;
end
V3 = thrustN./mflow + options.velocity; 
%V3 = (ones(10,1)*linspace(30,300,10))';
thrustact = mflow.*(V3-options.velocity);
thrusterr = thrustN - thrustact;
M3 = V3./speedsound; 
P3 = P2.*(1 + 0.5.*(gamma -1).*M3.^2).^(gamma./(gamma-1));
T3 = T2.*(1 + 0.5.*(gamma -1).*M3.^2);
%density3a = refproparray('D','T',T3,'P',P3,'OXYGEN','NITROGEN',[0.21,0.79]);
%density3 = density3a(1:10,1)*ones(1,10);
Pprop = mflow.*(P3-P2)./density; 

end
