clear; close all; clc;

% [file,folder] = uigetfile('.xlsx','Input Spectral Datasheet');

 

if exist("Spectral_Database.mat","file") == 2
    Spectral_data = importdata("Spectral_Database.mat");
else
    [file,folder] = uigetfile('.xlsx','Input Spectral Datasheet');
    Spectral_data = {};
    Spectral = importdata(file,folder);
    Spectral_data{1} = Spectral.data.H2O;
    Spectral_data{2} = Spectral.data.CO2;
    Spectral_data{3} = Spectral.data.CO;
    save('Spectral_database.mat',"Spectral_data");
end

T = 3000;
P = 1;
X = 0.03;
L = 10;
v_start = 2004;
v_end = 2012;
v_step = .001;
SimType = 'HITEMP Sim';
Species = 3 ; %1 H20, 2 CO2, 3 CO

Screensize = get(0, 'ScreenSize'); 
W = Screensize(3); 
H = Screensize(4);
% make centered
fig = uifigure('Position',[100,100,.8*W,.8*H]);
g = uigridlayout(fig);
g.RowHeight = {'1x','1x','25x','1x','1x'};
g.ColumnWidth = {'1x','1x','1x','1x','1x','1x','1x','1x','1x'};

% Dropdown 1
d1label = uilabel(g);
d1label.Text = 'Species';
d1label.Layout.Row = 1;
d1label.Layout.Column = 1;

d1 = uidropdown(g);
d1.Items = {'H2O','CO2','CO'};
d1.Layout.Row = 2;
d1.Layout.Column = 1;

% Button 2 - Temperature
b2label = uilabel(g);
b2label.Text = "Temperature (K)";
b2label.Layout.Row = 1;
b2label.Layout.Column = 2;

b2 = uieditfield(g, 'numeric');
b2.Value = T;
b2.Layout.Row = 2;
b2.Layout.Column = 2;

% Button 3 - Pressure
b3label = uilabel(g);
b3label.Text = "Pressure (atm)";
b3label.Layout.Row = 1;
b3label.Layout.Column = 3;

b3 = uieditfield(g, 'numeric');
b3.Value = P;
b3.Layout.Row = 2;
b3.Layout.Column = 3;

% Button 4 - Mole Fraction
b4label = uilabel(g);
b4label.Text = "Mole Fraction (X)";
b4label.Layout.Row = 1;
b4label.Layout.Column = 4;

b4 = uieditfield(g, 'numeric');
b4.Value = X;
b4.Layout.Row = 2;
b4.Layout.Column = 4;

% Button 5 - Path Length
b5label = uilabel(g);
b5label.Text = "Path Length (cm)";
b5label.Layout.Row = 1;
b5label.Layout.Column = 5;

b5 = uieditfield(g, 'numeric');
b5.Value = L;
b5.Layout.Row = 2;
b5.Layout.Column = 5;

% Button 6 - Start Wavenumber
b6label = uilabel(g);
b6label.Text = "Start Wavenumber";
b6label.Layout.Row = 1;
b6label.Layout.Column = 6;

b6 = uieditfield(g, 'numeric');
b6.Value = v_start;
b6.Layout.Row = 2;
b6.Layout.Column = 6;

% Button 7 - End Wavenumber
b7label = uilabel(g);
b7label.Text = "End Wavenumber";
b7label.Layout.Row = 1;
b7label.Layout.Column = 7;

b7 = uieditfield(g, 'numeric');
b7.Value = v_end;
b7.Layout.Row = 2;
b7.Layout.Column = 7;

% Button 8 - Step Size
b8label = uilabel(g);
b8label.Text = "Step Size";
b8label.Layout.Row = 1;
b8label.Layout.Column = 8;

b8 = uieditfield(g, 'numeric');
b8.Value = v_step;
b8.Layout.Row = 2;
b8.Layout.Column = 8;

ax = uiaxes(g);
ax.Layout.Row = [3,5]; 
ax.Layout.Column = [1,8];


[alpha_total,W_array] = SpectralSim(Spectral_data,T,P,X,L,v_start,v_end,v_step,SimType,Species);

plot (ax, W_array,alpha_total);
