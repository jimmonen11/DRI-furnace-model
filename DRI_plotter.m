close all

Fe2O3col = "#A13D2D";
Fe3O4col = "#654321";
FeOcol = 'k';
Fecol = '#808080';
H2col = [0 0.8 0];
H2Ocol = 'b';

z = linspace(0, h_furnace, n_furnace+2);

rho_b = out.Fe2O3_conc.data(end,:) + out.Fe3O4_conc.data(end,:) + out.FeO_conc.data(end,:) + out.Fe_conc.data(end,:); 

w_Fe2O3 = [out.Fe2O3_conc.data(end,:)./rho_b out.w_Fe2O3in.data(end) ];
w_Fe3O4 = [ out.Fe3O4_conc.data(end,:)./rho_b out.w_Fe3O4in.data(end)];
w_FeO = [ out.FeO_conc.data(end,:)./rho_b out.w_FeOin.data(end)];
w_Fe = [ out.Fe_conc.data(end,:)./rho_b out.w_Fein.data(end)];

gas_conc = out.H2_conc.data(end,:) + out.H2O_conc.data(end,:);
gas_conc_in = out.H2_concin.data(end) + out.H2O_concin.data(end);

x_H2 = [out.H2_concin.data(end)./gas_conc_in  out.H2_conc.data(end,:)./gas_conc];
x_H2O = [out.H2O_concin.data(end)./gas_conc_in out.H2O_conc.data(end,:)./gas_conc];

T_gin = out.T_g.data(end) -273;

T_s = out.T_s.data(end,:) -273;
T_g = [T_gin out.T_g.data(end,:) -273];

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
subplot(1,2,1)
box on
plot(w_Fe2O3, z(2:end) , 'linewidth', 6, 'color', Fe2O3col )
hold on
plot(w_Fe3O4, z(2:end), 'linewidth', 6, 'color', Fe3O4col)
plot(w_FeO, z(2:end), 'linewidth', 6, 'color', FeOcol)
plot(w_Fe, z(2:end), 'linewidth', 6, 'color', Fecol )
xlabel('Weight Fraction')
ylabel('Furnace Height (m)')

xlim([0, 1])
ylim([0, h_furnace])

legend('Fe_2O_3','Fe_3O_4', 'FeO', 'Fe', 'Location', 'southwest')
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

subplot(1,2,2)
box on
plot(x_H2, z(1:end-1), 'linewidth', 6, 'color', H2col )
hold on
plot(x_H2O, z(1:end-1), 'linewidth', 6, 'color', H2Ocol)
xlabel('Mole Fraction')
ylabel('Furnace Height (m)')

xlim([0, 1])
ylim([0, h_furnace])

legend('H_2','H_2O', 'Location', 'south')
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


