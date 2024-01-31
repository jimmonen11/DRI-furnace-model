close all

Fe2O3col = "#A13D2D";
Fe3O4col = "#654321";
FeOcol = 'k';
Fecol = '#808080';
H2col = [0 .8 0];
H2Ocol = 'b';
gascol = "#0072BD";
solidscol = "#D95319";


time = out.Fe2O3_conc.time;
time_interp = linspace(0,time(end), 1000);

z = linspace(0, h_furnace, n_furnace+2);

rho_b = out.Fe2O3_conc.data + out.Fe3O4_conc.data + out.FeO_conc.data + out.Fe_conc.data; 


w_Fe2O3 = out.Fe2O3_conc.data./rho_b;
w_Fe3O4 = out.Fe3O4_conc.data./rho_b;
w_FeO = out.FeO_conc.data./rho_b;
w_Fe = out.Fe_conc.data./rho_b;

gas_conc = out.H2_conc.data + out.H2O_conc.data;
gas_conc_in = out.H2_concin.data + out.H2O_concin.data;

x_H2 = out.H2_conc.data./gas_conc;
x_H2O = out.H2O_conc.data./gas_conc;

T_s = out.T_s.data -273;
T_g = out.T_g.data -273;

T_s_plot = zeros(length(time_interp), n_furnace);
T_g_plot = zeros(length(time_interp), n_furnace);

x_H2_plot = zeros(length(time_interp), n_furnace);
x_H2O_plot = zeros(length(time_interp), n_furnace);

w_Fe2O3_plot = zeros(length(time_interp), n_furnace) ;
w_Fe3O4_plot = zeros(length(time_interp), n_furnace) ;
w_FeO_plot = zeros(length(time_interp), n_furnace) ;
w_Fe_plot = zeros(length(time_interp), n_furnace) ;


for i =1:n_furnace
    properties = [x_H2(:,i), x_H2O(:,i), w_Fe2O3(:,i),  w_Fe3O4(:,i),  w_FeO(:,i), w_Fe(:,i), T_g(:,i), T_s(:,i), out.H2_concin.data, out.H2O_concin.data];
    vq = interp1(time, properties, time_interp);

    x_H2_plot(:,i) = vq(:,1);
    x_H2O_plot(:,i) = vq(:,2);

    w_Fe2O3_plot(:,i) = vq(:,3);
    w_Fe3O4_plot(:,i) = vq(:,4);
    w_FeO_plot(:,i) = vq(:,5);
    w_Fe_plot(:,i) = vq(:,6);

    T_g_plot(:,i) = vq(:,7);
    T_s_plot(:,i) = vq(:,8);

end

% x_H2_plot =  [ones(size(time_interp'))*out.H2_concin.data x_H2_plot];
% x_H2O_plot =  [ones(size(time_interp'))*out.H2O_concin.data x_H2O_plot];

x_H2_plot =  [vq(:, 9) x_H2_plot];
x_H2O_plot =  [vq(:,10) x_H2O_plot];

w_Fe2O3_plot =  [w_Fe2O3_plot ones(size(time_interp'))*out.w_Fe2O3in.data];
w_Fe3O4_plot =  [w_Fe3O4_plot ones(size(time_interp'))*out.w_Fe3O4in.data];
w_FeO_plot =  [w_FeO_plot ones(size(time_interp'))*out.w_FeOin.data];
w_Fe_plot =  [w_Fe_plot ones(size(time_interp'))*out.w_Fein.data];

T_g_plot = [ones(size(time_interp'))*(out.T_gin.data-273) T_g_plot];
T_s_plot = [T_s_plot ones(size(time_interp'))*(out.T_sin.data-273)];


%% Find Closest Index

[minval,hour1] = min(abs(time_interp-3600*2));
[minval,hour3] = min(abs(time_interp-3600*4));
[minval,hour10] = min(abs(time_interp-3600*11));

%%

widthsize = 3.5;
%figure('units','normalized','outerposition',[0 0 1 1])
figure('units','inches','innerposition',[-5 -5 20 200])


subplot(2,4,1)
%subplot(4,4,[1, 5])


box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('Fe_2O_3')
plot(w_Fe2O3_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col)
plot(w_Fe2O3_plot(hour1,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', '-.')
plot(w_Fe2O3_plot(hour3,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', '--')
plot(w_Fe2O3_plot(hour10,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Weight Fraction')
ylabel('Furnace Height (m)')
xlim([0, 1])
ylim([0, h_furnace])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

subplot(2,4,2)
%subplot(4,4,[2 6])

box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('Fe_3O_4')
plot(w_Fe3O4_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col)
plot(w_Fe3O4_plot(hour1,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', '-.')
plot(w_Fe3O4_plot(hour3,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', '--')
plot(w_Fe3O4_plot(hour10,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(2,4,3)
%subplot(4,4,[3 7])

box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('FeO')
plot(w_FeO_plot(1,:), z(2:end), 'linewidth', widthsize,'color', FeOcol)
plot(w_FeO_plot(hour1,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', '-.')
plot(w_FeO_plot(hour3,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', '--')
plot(w_FeO_plot(hour10,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(2,4,4)
%subplot(4,4,[4 8])

box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('Fe')
plot(w_Fe_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fecol)
plot(w_Fe_plot(hour1,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', '-.')
plot(w_Fe_plot(hour3,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', '--')
plot(w_Fe_plot(hour10,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

%% H2 H2O and Temp 

%figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,4,5)
%subplot(4,4,[9, 13])

box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('H_2')
plot(x_H2_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', H2col)
plot(x_H2_plot(hour1,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', '-.')
plot(x_H2_plot(hour3,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', '--')
plot(x_H2_plot(hour10,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Mole Fraction')
xlim([0.7, 1])
xticks([0.7:0.1:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

subplot(2,4,6)
%subplot(4,4,[10 14])
box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('H_2O')
plot(x_H2O_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol)
plot(x_H2O_plot(hour1,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', '-.')
plot(x_H2O_plot(hour3,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', '--')
plot(x_H2O_plot(hour10,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'south')

xlabel('Weight Fraction')
xlim([0, 0.3])
xticks([0:0.1:0.3]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(2,4,7)
%subplot(4,4,[11 15])
box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('Gas')
plot(T_g_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', gascol)
plot(T_g_plot(hour1,:), z(1:end-1), 'linewidth', widthsize,'color', gascol, 'linestyle', '-.')
plot(T_g_plot(hour3,:), z(1:end-1), 'linewidth', widthsize,'color', gascol, 'linestyle', '--')
plot(T_g_plot(hour10,:), z(1:end-1), 'linewidth', widthsize,'color', gascol, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'southwest')

ylabel('Furnace Height (m)')
xlabel('Temperature (^oC)')
xlim([250, 850]);
xticks([250:100:850]);


H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)


subplot(2,4,8)
%subplot(4,4,[12 16])

box on
hold on
set(gca,'FontWeight', 'bold','FontSize',18)

title('Solids')
plot(T_s_plot(1,:), z(2:end), 'linewidth', widthsize,'color', solidscol)
plot(T_s_plot(hour1,:), z(2:end), 'linewidth', widthsize,'color', solidscol, 'linestyle', '-.')
plot(T_s_plot(hour3,:), z(2:end), 'linewidth', widthsize,'color', solidscol, 'linestyle', '--')
plot(T_s_plot(hour10,:), z(2:end), 'linewidth', widthsize,'color', solidscol, 'linestyle', ':')

legend('t = 0 hrs','t = 1 hrs', 't = 3 hrs', 't = 10 hrs', 'Location', 'southwest')

xlabel('Temperature (^oC)')
xlim([250, 850]);
xticks([250:100:850]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',18)

