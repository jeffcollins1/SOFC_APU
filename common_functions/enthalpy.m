function [Out1,Out2] = enthalpy(varargin) % enthalpy (h) and sensible enthalpy(h_s), if species and flow are provided it returns H (total enthalpy(kW)) and H_s Sensible enthalpy (kW).
% returns both total and sensible enthalpy in kJ/kmol, or in rate form kJ/s
% Option 1: provide a vector of temperatures and it returns total and specific enthalpy at those temperatures for all species CH4, CO, CO2, H2, H2O, N2, O2, C, NO, OH, H
% Option 2: provide a structure where __.T coresponds to temperature, ___.CH4 coresponds to the flow rate of methane ____.H2 to the flow rate of hydrogen...
% Option 2, returns the rate of energy flow in (kJ/s)

if isfield(varargin{1},'T')
    Inlet = varargin{1};
    T = Inlet.T;
else
    T = varargin{1};
end

CH4a = [85.81217,11.26467,-2.114146,0.138190,-26.42221,-153.5327;];
CH4b = [-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376;];
CH4 = (T>1300)*CH4a+(T<=1300)*CH4b;

COa = [35.15070,1.300095,-.205921,0.013550,-3.282780,-127.8375;];
COb = [25.56759,6.096130,4.054656,-2.671301,0.131021,-118.0089;];
CO = (T>1300)*COa+(T<=1300)*COb;

CO2a = [58.16639,2.720074,-0.492289,0.038844,-6.447293,-425.9186;];
CO2b = [24.99735,55.18696,-33.69137,7.948387,-0.136638,-403.6075;];
CO2 = (T>1200)*CO2a+(T<=1200)*CO2b;

H2a = [43.413560,-4.293079,1.272428,-.096876,-20.533862,-38.515158;];
H2b = [18.563083,12.257357,-2.859786,0.268238,1.977990,-1.147438;];
H2c = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797;];
H2 = (T>2500)*H2a+(T<=2500).*(T>1000)*H2b+(T<=1000)*H2c;

H2Oa = [41.96426,8.622053,-1.499780,0.098119,-11.15764,-272.1797;];
H2Ob = [30.09200,6.832514,6.793435,-2.534480,0.082139,-250.8810;];
H2O = (T>1700)*H2Oa+(T<=1700)*H2Ob;

N2a = [35.51872,1.128728,-0.196103,0.014662,-4.553760,-18.97091;];
N2b = [19.50583,19.88705,-8.598535,1.369784,0.527601,-4.935202;];
N2c = [28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914;];
N2 = (T>2000)*N2a+(T<=2000).*(T>500)*N2b+(T<=500)*N2c;

O2a = [20.91111,10.72071,-2.020498,0.146449,9.245722,5.337651;];
O2b = [30.03235,8.772972,-3.988133,0.788313,-0.741599,-11.32468;];
O2c = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471;];
O2 = (T>2000)*O2a+(T<=2000).*(T>700)*O2b+(T<=700)*O2c;

C = [21.1751,-0.812428,0.448537,-0.043256,-0.013103,710.347,183.8734];

NOa = [35.99169,0.95717,-0.148032,0.009974,-3.004088,73.10787,246.1619];
NOb = [23.83491,12.58878,-1.139011,-1.497459,0.21419,83.35783,237.1219];
NO = (T>1200)*NOa+(T<=1200)*NOb;

OHa = [28.74701,4.7144,-0.814725,0.054748,-2.747829,26.41439,214.1166];
OHb = [32.27768,-11.36291,13.60545,-3.846486,-0.001335,29.75113,225.5783];
OH = (T>1300)*OHa+(T<=1300)*OHb;

H =[20.78603,4.85E-10,-1.58E-10,1.53E-11,3.20E-11,2.12E+02,1.40E+02];

C2H6 =[6.9,172.7,-64.06,7.285,9.173];
C3H8 =[-4.04,304.8,-157.2,31.74,11.05];
C6H6 =[-50.24,568.2244,-442.503,134.5489,6.6206];

T1 = (T/1000);
T2 = ((T/1000).^2)/2; 
T3 = ((T/1000).^3)/3; 
T4 = ((T/1000).^4)/4;
T5 = -1./(T/1000); 
T6 = ones(length(T),1);

h.CH4 = sum([T1.*CH4(:,1), T2.*CH4(:,2), T3.*CH4(:,3), T4.*CH4(:,4), T5.*CH4(:,5), T6.*CH4(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.CO = sum([T1.*CO(:,1), T2.*CO(:,2), T3.*CO(:,3), T4.*CO(:,4), T5.*CO(:,5), T6.*CO(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.CO2 = sum([T1.*CO2(:,1), T2.*CO2(:,2), T3.*CO2(:,3), T4.*CO2(:,4), T5.*CO2(:,5), T6.*CO2(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.H2 = sum([T1.*H2(:,1), T2.*H2(:,2), T3.*H2(:,3), T4.*H2(:,4), T5.*H2(:,5), T6.*H2(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.H2O = sum([T1.*H2O(:,1), T2.*H2O(:,2), T3.*H2O(:,3), T4.*H2O(:,4), T5.*H2O(:,5), T6.*H2O(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.N2 = sum([T1.*N2(:,1), T2.*N2(:,2), T3.*N2(:,3), T4.*N2(:,4), T5.*N2(:,5), T6.*N2(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.O2 = sum([T1.*O2(:,1), T2.*O2(:,2), T3.*O2(:,3), T4.*O2(:,4), T5.*O2(:,5), T6.*O2(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.C = sum([T1.*C(:,1), T2.*C(:,2), T3.*C(:,3), T4.*C(:,4), T5.*C(:,5), T6.*C(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.NO = sum([T1.*NO(:,1), T2.*NO(:,2), T3.*NO(:,3), T4.*NO(:,4), T5.*NO(:,5), T6.*NO(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.OH = sum([T1.*OH(:,1), T2.*OH(:,2), T3.*OH(:,3), T4.*OH(:,4), T5.*OH(:,5), T6.*OH(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.H = sum([T1.*H(:,1), T2.*H(:,2), T3.*H(:,3), T4.*H(:,4), T5.*H(:,5), T6.*H(:,6)],2)*1000; %convert from kJ/mol to kJ/kmol
h.C2H6 = sum([T1.*C2H6(:,1), T2.*C2H6(:,2), T3.*C2H6(:,3), T4.*C2H6(:,4)],2)*1000 - 9173.70215 -83700; %convert from kJ/mol to kJ/kmol
h.C3H8 = sum([T1.*C3H8(:,1), T2.*C3H8(:,2), T3.*C3H8(:,3), T4.*C3H8(:,4)],2)*1000 - 11005.6940 - 104700; %convert from kJ/mol to kJ/kmol
h.C6H6 = sum([T1.*C6H6(:,1), T2.*C6H6(:,2), T3.*C6H6(:,3), T4.*C6H6(:,4)],2)*1000 - 6620.6417 - 82930; %convert from kJ/mol to kJ/kmol

h_s.CH4 = h.CH4+74873.1;
h_s.CO = h.CO+110527.1;
h_s.CO2 = h.CO2+393522.4;
h_s.H2 = h.H2;
h_s.H2O = h.H2O+241826.4; %water vapor
h_s.N2 = h.N2;
h_s.O2 = h.O2;
h_s.C = h.C-716669;
h_s.NO = h.NO-90291.14;
h_s.OH = h.OH-38987.06;
h_s.H = h.H-218194;
h_s.C2H6 = h.C2H6+83700;
h_s.C3H8 = h.C3H8+104700;
h_s.C6H6 = h.C6H6+82930; %benzene vapor

if ~exist('Inlet','var')
    Out1 = h;
    Out2 = h_s;
else
    Out1 = 0;
    Out2 = 0;
    speciesName = fieldnames(Inlet);
    for i = 1:1:length(speciesName)
        if ~strcmp(speciesName{i},'T')
            Out1 = Out1 + h.(speciesName{i}).*Inlet.(speciesName{i});
            Out2 = Out2 + h_s.(speciesName{i}).*Inlet.(speciesName{i});
        end
    end
end