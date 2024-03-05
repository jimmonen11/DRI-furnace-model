close all

Fe2O3col = "#A13D2D";
Fe3O4col = "#654321";
FeOcol = 'k';
Fecol = '#808080';
Ccol = '#ebba34' ;
Gancol = '#9e6a16' ;

H2col = [0 .8 0];
H2Ocol = 'b';
COcol = '#44734c';
CO2col = '#98afb8';
CH4col = '#36cfb9';
N2col = '#87CEEB';

gasTcol = '#0072BD';
gasmcol = '#077b7d';

solidsTcol = '#D95319';
solidsmcol = '#7d072d' ;

fontsize = 16;
fontsizetit = 19;

z = linspace(0, h_furnace, n_furnace+2);

hour_aim = 35;
time = out.w_Fe2O3.time(end);
time_interp = linspace(0,time(end), 1000);

[minval,hour_id] = min(abs(time_interp-3600*hour_aim));
hour_id = length(out.w_Fe.data(:,1));

%%
w_Fe2O3 =  out.w_Fe2O3.data(hour_id,:);
w_Fe3O4 = out.w_Fe3O4.data(hour_id,:);
w_FeO =  out.w_FeO.data(hour_id,:);
w_Fe =  out.w_Fe.data(hour_id,:) ;
w_C =  out.w_C.data(hour_id,:) ;

x_H2 = out.x_H2.data(hour_id,:);
x_H2O = out.x_H2O.data(hour_id,:);
x_N2 = out.x_N2.data(hour_id,:);
x_CO = out.x_CO.data(hour_id,:);
x_CO2 = out.x_CO2.data(hour_id,:);
x_CH4 = out.x_CH4.data(hour_id,:);


T_g = out.T_g.data(hour_id,:)-273;
T_s = out.T_s.data(hour_id,:)-273;

%%

figure(1)
subplot(1,3,1)
box on
plot(w_Fe2O3, z(2:end) , 'linewidth', 6, 'color', Fe2O3col )
hold on
plot(w_Fe3O4, z(2:end), 'linewidth', 6, 'color', Fe3O4col)
plot(w_FeO, z(2:end), 'linewidth', 6, 'color', FeOcol)
plot(w_Fe, z(2:end), 'linewidth', 6, 'color', Fecol )
%plot(w_C, z(1:end), 'linewidth', 3, 'color', 'k', 'LineStyle', '--' )
xlabel('Weight Fraction')
ylabel('Furnace Height (m)')

xlim([0, 1])
ylim([0, h_furnace])

legend('Fe_2O_3','Fe_3O_4', 'FeO', 'Fe', 'Location', 'best')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

subplot(1,3,2)
box on
plot(x_H2, z(1:end-1), 'linewidth', 6, 'color', H2col )
hold on
plot(x_H2O, z(1:end-1), 'linewidth', 6, 'color', H2Ocol)
plot(x_CO, z(1:end-1), 'linestyle', ':','linewidth', 6, 'color', COcol)
plot(x_CO2, z(1:end-1), 'linestyle', ':', 'linewidth', 6, 'color', CO2col)
plot(x_CH4, z(1:end-1), 'linestyle', ':', 'linewidth', 6, 'color', CH4col)
xlabel('Mole Fraction')
ylabel('Furnace Height (m)')

xlim([0, 0.6])
ylim([0, h_furnace])

legend('H_2','H_2O', 'CO', 'CO_2', 'CH_4', 'Location', 'best')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(1,3,3)
box on
hold on
plot(T_g, z(1:end-1), 'linewidth', 6,'color', gasTcol)
plot(T_s, z(2:end), 'linewidth', 6,'color', solidsTcol)

xlabel('Temperature (^oC)')
ylabel('Furnace Height (m)')

xlim([250, 1000]);
ylim([0, h_furnace])

legend('T_g','T_s', 'Location', 'south')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)



