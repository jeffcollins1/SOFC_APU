%% Geometery for planar SOFC 
L = 1e-1; %cell length
W = 1e-1; %Cell width
% c_area = L*W; %area of each SOFC cell
t_w = .0125;% vessel wall thickness
d_w = 1800; %kg/m^3 for vessel wall
t_FC = 170e-6; %anode supported fuel cell thickness
d_FC = 6090; %kg/m^3 for YSZ
t_inter = 0.2e-3; %thickness of interconnect
d_inter = 7700; %kg/m^3 for interconnect material
h_channel = 500e-6; %height of channels
% c_w_ratio = 0.75; %channel to wall ratio 
w_seal = 5e-3; %width of seal
t_seal = 200e-6; %compressed seal thickness
d_seal = 1200;% density of compressed/sintered thermiculit 870

T_SOFC = 750; %Temperature
T_wall = 85; %Temperature
k_insul = 0.11; %W/m*K thermal conductivity of insulation
d_insul = 64; %Density of insulation
ASR = 0.25;%area specific resistance
active_area = (L-2*w_seal)*(W-2*w_seal)*1e4;% active area in cm^2
voltage = 0.8;
current = active_area*(1.1-voltage)/ASR;
power = current*voltage; %W per cell
heat_gen = 1.25*current-power; %heat generation (use 1.25 for H2, 0.95 for CH4)
h_loss = .05; % Percent of heat generation lost therough vessel wall


h_repeat = h_channel+2*t_inter+2*t_seal+t_FC;
channels = floor(W/(3*h_channel));%%assume trapezoid channels:
W_inter = channels*4.2361*h_channel + (W-channels*3*h_channel);
inter_mass = t_inter*W_inter*L*d_inter;

FC_mass = L*W*t_FC*d_FC;

seal_mass = w_seal*t_seal*(2*W+2*L-4*w_seal)*d_seal;

%radial heat transfer
d_cell = W*1.1; %'effective' diameter of cell

R = (T_SOFC - T_wall)/(h_loss*heat_gen); %Resistance of insulation R = (T1-T2)/Q
diam_inner = d_cell*exp(2*pi*k_insul*h_repeat*R);
insulation_mass = h_repeat*(pi*diam_inner^2/4 - L*W)*d_insul;

diam_outer = diam_inner+t_w;
vessel_mass = h_repeat*(pi*diam_outer^2/4 - pi*diam_inner^2/4)*d_w;

total = inter_mass + FC_mass + 2*seal_mass + insulation_mass + vessel_mass;
pow_den = power/total/1000;
%% arranging 4 cells per layer
n = 4;
d_cell = W*2.5; %'effective' diameter of cell

% n = 13;
% d_cell = W*4; %'effective' diameter of cell

R = (T_SOFC - T_wall)/(h_loss*n*heat_gen); %Resistance of insulation R = (T1-T2)/Q
diam_inner = d_cell*exp(2*pi*k_insul*h_repeat*R);
insulation_mass = h_repeat*(pi*diam_inner^2/4 - n*L*W)*d_insul;

diam_outer = diam_inner+t_w;
vessel_mass = h_repeat*(pi*diam_outer^2/4 - pi*diam_inner^2/4)*d_w;

total2 = (n*inter_mass + n*FC_mass + n*2*seal_mass + insulation_mass + vessel_mass)/n;
