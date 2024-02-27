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


time = out.Fe2O3_conc.time;
time_interp = linspace(0,time(end), 1000);

z = linspace(0, h_furnace, n_furnace+2);

T_s = out.T_s.data -273;
T_g = out.T_g.data -273;


T_s_plot = zeros(length(time_interp), n_furnace+1) ;
T_g_plot = zeros(length(time_interp), n_furnace+1) ;


x_H2_plot = zeros(length(time_interp), n_furnace+1);
x_H2O_plot = zeros(length(time_interp), n_furnace+1);
x_CO2_plot = zeros(length(time_interp), n_furnace+1);
x_CO_plot = zeros(length(time_interp), n_furnace+1);
x_CH4_plot = zeros(length(time_interp), n_furnace+1);
x_N2_plot = zeros(length(time_interp), n_furnace+1);

m_g_plot = zeros(length(time_interp), n_furnace+1) ;



w_Fe2O3_plot = zeros(length(time_interp), n_furnace+1) ;
w_Fe3O4_plot = zeros(length(time_interp), n_furnace+1) ;
w_FeO_plot = zeros(length(time_interp), n_furnace) ;
w_Fe_plot = zeros(length(time_interp), n_furnace+1) ;
w_C_plot = zeros(length(time_interp), n_furnace+1) ;
w_Gan_plot = zeros(length(time_interp), n_furnace+1) ;
m_s_plot = zeros(length(time_interp), n_furnace+1) ;


for i =1:n_furnace+1
    properties = [out.x_H2.data(:,i), out.x_H2O.data(:,i), out.x_CO.data(:,i), out.x_CO2.data(:,i), out.x_CH4.data(:,i), out.x_N2.data(:,i),...
        out.w_Fe2O3.data(:,i),  out.w_Fe3O4.data(:,i),  out.w_FeO.data(:,i), out.w_Fe.data(:,i), out.w_C.data(:,i), out.w_Gan.data(:,i), ...
        T_g(:,i), T_s(:,i), out.m_g.data(:,i), out.m_s.data(:,i)];
    
    vq = interp1(time, properties, time_interp);

    x_H2_plot(:,i) = vq(:,1);
    x_H2O_plot(:,i) = vq(:,2);
    x_CO_plot(:,i) = vq(:,3);
    x_CO2_plot(:,i) = vq(:,4);
    x_CH4_plot(:,i) = vq(:,5);
    x_N2_plot(:,i) = vq(:,6);

    w_Fe2O3_plot(:,i) = vq(:,7);
    w_Fe3O4_plot(:,i) = vq(:,8);
    w_FeO_plot(:,i) = vq(:,9);
    w_Fe_plot(:,i) = vq(:,10);
    w_C_plot(:,i) = vq(:,11);
    w_Gan_plot(:,i) = vq(:,12);

    T_g_plot(:,i) = vq(:,13);
    T_s_plot(:,i) = vq(:,14);

    m_g_plot(:,i) = vq(:,15);
    m_s_plot(:,i) = vq(:,16);

end

%% Find Closest Index 

%Choose where want time 'slices' to be

time1 = 0.5; %hours
time2 = 1; %hours
time3 = 1.5; % hours

[minval,split1] = min(abs(time_interp-(3600*time1 + 1*3600)));
[minval,split2] = min(abs(time_interp-(3600*time2 + 1*3600)));
[minval,split3] = min(abs(time_interp-(3600*time3 + 1*3600)));

%% Solids

widthsize = 3.5;
figure('units','inches','innerposition',[-5 -5 20 200])


subplot(2,4,1)
box on
hold on

plot(w_Fe2O3_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col)
plot(w_Fe2O3_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', '-.')
plot(w_Fe2O3_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', '--')
plot(w_Fe2O3_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', Fe2O3col, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')


xlabel('Weight Fraction')
ylabel('Furnace Height (m)')
xlim([0, 1])
ylim([0, h_furnace])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Fe_2O_3', 'FontSize', fontsizetit)


subplot(2,4,2)
box on
hold on

title('Fe_3O_4')
plot(w_Fe3O4_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col)
plot(w_Fe3O4_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', '-.')
plot(w_Fe3O4_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', '--')
plot(w_Fe3O4_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', Fe3O4col, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Fe_3O_4', 'FontSize', fontsizetit)


subplot(2,4,3)
box on
hold on

plot(w_FeO_plot(1,:), z(2:end), 'linewidth', widthsize,'color', FeOcol)
plot(w_FeO_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', '-.')
plot(w_FeO_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', '--')
plot(w_FeO_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', FeOcol, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('FeO', 'FontSize', fontsizetit)


subplot(2,4,4)
box on
hold on

plot(w_Fe_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Fecol)
plot(w_Fe_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', '-.')
plot(w_Fe_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', '--')
plot(w_Fe_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', Fecol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Weight Fraction')
xlim([0, 1])
xticks([0:0.25:1]);
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Fe', 'FontSize',fontsizetit)


subplot(2,4,5)
box on
hold on

plot(w_C_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Ccol)
plot(w_C_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', Ccol, 'linestyle', '-.')
plot(w_C_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', Ccol, 'linestyle', '--')
plot(w_C_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', Ccol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

ylabel('Furnace Height (m)')
xlabel('Weight Fraction')
xlim([0, 0.02])
xticks([0:0.005:0.02]);
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('C', 'FontSize', fontsizetit)


subplot(2,4,6)
box on
hold on

plot(w_Gan_plot(1,:), z(2:end), 'linewidth', widthsize,'color', Gancol)
plot(w_Gan_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', Gancol, 'linestyle', '-.')
plot(w_Gan_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', Gancol, 'linestyle', '--')
plot(w_Gan_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', Gancol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Weight Fraction')
%xlim([0, 0.05])
%xticks([0:0.01:0.05]);
H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Gangue', 'FontSize', fontsizetit)


subplot(2,4,7)
box on
hold on

plot(T_s_plot(1,:), z(2:end), 'linewidth', widthsize,'color', solidsTcol)
plot(T_s_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', solidsTcol, 'linestyle', '-.')
plot(T_s_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', solidsTcol, 'linestyle', '--')
plot(T_s_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', solidsTcol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Temperature (^oC)')
xlim([250, 1000]);
%xticks([250:100:850]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Solids Temperature', 'FontSize', fontsizetit)


subplot(2,4,8)
box on
hold on

plot(m_s_plot(1,:), z(2:end), 'linewidth', widthsize,'color', solidsmcol)
plot(m_s_plot(split1,:), z(2:end), 'linewidth', widthsize,'color', solidsmcol, 'linestyle', '-.')
plot(m_s_plot(split2,:), z(2:end), 'linewidth', widthsize,'color', solidsmcol, 'linestyle', '--')
plot(m_s_plot(split3,:), z(2:end), 'linewidth', widthsize,'color', solidsmcol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Mass Flow (kg/s)')
%xlim([250, 850]);
%xticks([250:100:850]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Solids Mass Flow', 'fontsize', fontsizetit)



%%

%% H2 H2O CO CO2 

figure(2)

subplot(2,4,1)
box on
hold on

title('H_2')
plot(x_H2_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', H2col)
plot(x_H2_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', '-.')
plot(x_H2_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', '--')
plot(x_H2_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', H2col, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

ylabel('Furnace Height (m)')
xlabel('Mole Fraction')
xlim([0.3, 0.9])
%xticks([0.3:0.1:0.9]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('H_2', 'fontsize', fontsizetit)


subplot(2,4,2)
box on
hold on

plot(x_H2O_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol)
plot(x_H2O_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', '-.')
plot(x_H2O_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', '--')
plot(x_H2O_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', H2Ocol, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Mole Fraction')
%xlim([0, 0.4])
%xticks([0:0.1:0.4]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('H_2O', 'fontsize', fontsizetit)


subplot(2,4,3)
box on
hold on

plot(x_CO_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', COcol)
plot(x_CO_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', COcol, 'linestyle', '-.')
plot(x_CO_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', COcol, 'linestyle', '--')
plot(x_CO_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', COcol, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'south')

xlabel('Mole Fraction')
xlim([0, 0.4])
%xticks([0:0.1:0.4]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize', fontsize)
title('CO', 'FontSize', fontsizetit)


subplot(2,4,4)
box on
hold on

plot(x_CO2_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', CO2col)
plot(x_CO2_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', CO2col, 'linestyle', '-.')
plot(x_CO2_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', CO2col, 'linestyle', '--')
plot(x_CO2_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', CO2col, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'south')

xlabel('Mole Fraction')
xlim([0, 0.2])
%xticks([0:0.1:0.4]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('CO_2', 'FontSize',fontsizetit)

subplot(2,4,5)
box on
hold on

plot(x_CH4_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', CH4col)
plot(x_CH4_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', CH4col, 'linestyle', '-.')
plot(x_CH4_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', CH4col, 'linestyle', '--')
plot(x_CH4_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', CH4col, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

ylabel('Furnace Height (m)')
xlabel('Mole Fraction')
%xlim([0, 0.4])
%xticks([0:0.1:0.4]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('CH_4', 'FontSize',fontsizetit)



subplot(2,4,6)
box on
hold on

plot(x_N2_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', N2col)
plot(x_N2_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', N2col, 'linestyle', '-.')
plot(x_N2_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', N2col, 'linestyle', '--')
plot(x_N2_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', N2col, 'linestyle', ':')


legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Mole Fraction')
%xlim([0, 0.4])
%xticks([0:0.1:0.4]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('N_2', 'FontSize',fontsizetit)


subplot(2,4,7)
box on
hold on

plot(T_g_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', gasTcol)
plot(T_g_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', gasTcol, 'linestyle', '-.')
plot(T_g_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', gasTcol, 'linestyle', '--')
plot(T_g_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', gasTcol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Temperature (^oC)')
xlim([250, 1000]);
%xticks([250:100:850]);


H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Gas Temperature', 'FontSize',fontsizetit)


subplot(2,4,8)
box on
hold on

plot(m_g_plot(1,:), z(1:end-1), 'linewidth', widthsize,'color', gasmcol)
plot(m_g_plot(split1,:), z(1:end-1), 'linewidth', widthsize,'color', gasmcol, 'linestyle', '-.')
plot(m_g_plot(split2,:), z(1:end-1), 'linewidth', widthsize,'color', gasmcol, 'linestyle', '--')
plot(m_g_plot(split3,:), z(1:end-1), 'linewidth', widthsize,'color', gasmcol, 'linestyle', ':')

legend('t = 0 hrs', append('t = ', num2str(time1, '%.1f'), ' hrs'),...
     append('t = ', num2str(time2, '%.0f'), ' hrs'),...
     append('t = ', num2str(time3, '%.0f'), ' hrs'), 'Location', 'best')

xlabel('Mass Flow (kg/s)')
xlim([10, 60]);
%xticks([250:100:850]);

H = gca;
grid on
H.LineWidth = 3; %change to the desired value   
set(gca,'FontWeight', 'bold','FontSize',fontsize)
title('Gas Mass Flow', 'FontSize',fontsizetit)

