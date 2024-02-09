close all

Fe2O3col = "#A13D2D";
Fe3O4col = "#654321";
FeOcol = 'k';
Fecol = '#808080';
H2col = [0 0.8 0];

H2Ocol = 'b';

COcol = '#44734c';
CO2col = '#98afb8';

gascol = "#0072BD";
solidscol = "#D95319";


hour_aim = 35;
time = out.w_Fe2O3.time(end);
time_interp = linspace(0,time(end), 1000);

[minval,hour_id] = min(abs(time_interp-3600*hour_aim));

hour_id = length(out.w_Fe.data(:,1));


z = linspace(0, h_furnace, n_furnace+2);


%rho_b = out.Fe2O3_conc.data(hour_id,:) + out.Fe3O4_conc.data(hour_id,:) + out.FeO_conc.data(hour_id,:) + out.Fe_conc.data(hour_id,:); 

w_Fe2O3 = [out.w_Fe2O3.data(hour_id,:) out.w_Fe2O3in.data(1) ];
w_Fe3O4 = [ out.w_Fe3O4.data(hour_id,:) out.w_Fe3O4in.data(1)];
w_FeO = [ out.w_FeO.data(hour_id,:) out.w_FeOin.data(1)];
w_Fe = [ out.w_Fe.data(hour_id,:) out.w_Fein.data(1)];

Xred = [ out.Xred.data(hour_id,:) out.w_Fein.data(1)];


%gas_conc = out.H2_conc.data(hour_id,:) + out.H2O_conc.data(hour_id,:);
%gas_conc_in = out.H2_concin.data(1) + out.H2O_concin.data(1);

x_H2 = [out.x_H2in.data  out.x_H2.data(hour_id,:)];
x_H2O = [out.x_H2Oin.data out.x_H2O.data(hour_id,:)];
x_N2 = [out.x_N2in.data out.x_N2.data(hour_id,:)];
x_CO = [out.x_COin.data out.x_CO.data(hour_id,:)];
x_CO2 = [out.x_CO2in.data out.x_CO2.data(hour_id,:)];

T_g = [(out.T_gin.data(1)-273) out.T_g.data(hour_id,:)-273];
T_s = [out.T_s.data(hour_id,:)-273 (out.T_sin.data(1)-273)];

%Y = 1000 * [4.1962 3.5087 2.8783 2.3026 1.7775 1.2967 0.8516 0.4318]';
%Y = interp1(Y, linspace(1,8))';

% figure (1)
% 
% subplot(1,4,1)
% imagesc (flip(w_Fe2O3'))
% axis off
% 
% subplot(1,4,2)
% imagesc (flip(w_Fe3O4'))
% axis off
% 
% subplot(1,4,3)
% imagesc (flip(w_FeO'))
% axis off
% subplot(1,4,4)
% imagesc (flip(w_Fe'))
% axis off
% colorbar

figure(2)
subplot(1,3,1)
box on
plot(w_Fe2O3, z(2:end) , 'linewidth', 6, 'color', Fe2O3col )
hold on
plot(w_Fe3O4, z(2:end), 'linewidth', 6, 'color', Fe3O4col)
plot(w_FeO, z(2:end), 'linewidth', 6, 'color', FeOcol)
plot(w_Fe, z(2:end), 'linewidth', 6, 'color', Fecol )
plot(Xred, z(2:end), 'linewidth', 3, 'color', 'k', 'LineStyle', '--' )
xlabel('Weight Fraction')
ylabel('Furnace Height (m)')

xlim([0, 1])
ylim([0, h_furnace])

legend('Fe_2O_3','Fe_3O_4', 'FeO', 'Fe', 'Location', 'southwest')
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
xlabel('Mole Fraction')
ylabel('Furnace Height (m)')

xlim([0, 0.6])
ylim([0, h_furnace])

legend('H_2','H_2O', 'CO', 'CO_2', 'Location', 'south')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(1,3,3)
box on
hold on
plot(T_g, z(1:end-1), 'linewidth', 6,'color', gascol)
plot(T_s, z(2:end), 'linewidth', 6,'color', solidscol)

xlabel('Temperature (^oC)')
ylabel('Furnace Height (m)')

xlim([250, 1000]);
ylim([0, h_furnace])

legend('T_g','T_s', 'Location', 'south')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

% figure(3)
% plot(x_H2, z, 'linewidth', 3)
% hold on
% plot(x_H2O, z, 'linewidth', 3)
% 
% xlim([0, 1])
% 
% legend('H2','H2O')
% 
% 
% figure(4)
% plot(T_s, z, 'linewidth', 3)
% hold on
% plot(T_g, z, 'linewidth', 3)
% 
% %xlim([0, 1])
% 
% legend('Solid','Gas')


