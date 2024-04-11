close all

Fe2O3col = "#A13D2D";
Fe3O4col = "#654321";
FeOcol = 'k';
Fecol = '#808080';
H2col = [0 .8 0];
H2Ocol = 'b';
COcol = '#44734c';
CO2col = '#98afb8';
CH4col = '#36cfb9';
gascol = "#0072BD";
solidscol = "#D95319";


time = out.Fe2O3_conc.time;
time_interp = linspace(0,time(end), 250);

z = linspace(-dz, h_furnace + dz, n_furnace+2);

T_s = out.T_s.data -273;
T_g = out.T_g.data -273;


T_s_plot = zeros(length(time_interp), n_furnace+1) ;
T_g_plot = zeros(length(time_interp), n_furnace+1) ;


x_H2_plot = zeros(length(time_interp), n_furnace+1);
x_H2O_plot = zeros(length(time_interp), n_furnace+1);
x_CO2_plot = zeros(length(time_interp), n_furnace+1);
x_CO_plot = zeros(length(time_interp), n_furnace+1);
x_CH4_plot = zeros(length(time_interp), n_furnace+1);


w_Fe2O3_plot = zeros(length(time_interp), n_furnace+1) ;
w_Fe3O4_plot = zeros(length(time_interp), n_furnace+1) ;
w_FeO_plot = zeros(length(time_interp), n_furnace) ;
w_Fe_plot = zeros(length(time_interp), n_furnace+1) ;
w_C_plot = zeros(length(time_interp), n_furnace+1) ;
w_Gan_plot = zeros(length(time_interp), n_furnace+1) ;
Xred_plot = zeros(length(time_interp), n_furnace+1) ;



for i =1:n_furnace+1
    properties = [out.x_H2.data(:,i), out.x_H2O.data(:,i), out.x_CO.data(:,i), out.x_CO2.data(:,i), out.x_CH4.data(:,i),...
        out.w_Fe2O3.data(:,i),  out.w_Fe3O4.data(:,i),  out.w_FeO.data(:,i), out.w_Fe.data(:,i), out.w_C.data(:,i), out.w_Gan.data(:,i), out.Xred.data(:,i), ...
        T_g(:,i), T_s(:,i)];
    
    vq = interp1(time, properties, time_interp);

    x_H2_plot(:,i) = vq(:,1);
    x_H2O_plot(:,i) = vq(:,2);
    x_CO_plot(:,i) = vq(:,3);
    x_CO2_plot(:,i) = vq(:,4);
    x_CH4_plot(:,i) = vq(:,5);

    w_Fe2O3_plot(:,i) = vq(:,6);
    w_Fe3O4_plot(:,i) = vq(:,7);
    w_FeO_plot(:,i) = vq(:,8);
    w_Fe_plot(:,i) = vq(:,9);
    w_C_plot(:,i) = vq(:,10);
    w_Gan_plot(:,i) = vq(:,11);
    Xred_plot(:,i) = vq(:,12);

    T_g_plot(:,i) = vq(:,13);
    T_s_plot(:,i) = vq(:,14);


end


figure('units','normalized','outerposition',[0 0 1 1])
gif('stepcase.gif')

[minval,hour1] = min(abs(time_interp-3600));

time_end = 4; %hours
[minval, hour_end] = min(abs(time_interp-(3600*time_end + 1*3600)));


for i = hour1:hour_end
    
    subplot(1,3,1)
    box on
    xlim([0, 1])
    hold on
    set(gca,'FontWeight', 'bold','FontSize',18)
    plot(w_Fe2O3_plot(i,:), z(2:end), 'linewidth', 6,'color', Fe2O3col)
    plot(w_Fe3O4_plot(i,:), z(2:end), 'linewidth', 6,'color', Fe3O4col)
    plot(w_FeO_plot(i,:), z(2:end), 'linewidth', 6, 'color',FeOcol)
    plot(w_Fe_plot(i,:), z(2:end), 'linewidth', 6,'color', Fecol)
    
    xlabel('Weight Fraction')
    ylabel('Furnace Height (m)')
    ylim([0, h_furnace])

    xticks([0:0.25:1]);
    legend('Fe_2O_3','Fe_3O_4', 'FeO', 'Fe','Location', 'east')
    H = gca;
    grid on
    H.LineWidth = 3; %change to the desired value   
    set(gca,'FontWeight', 'bold','FontSize',18)
    

    subplot(1,3,2)
    box on
    hold on
    xlim([0, 1])
    plot(x_H2_plot(i,:), z(1:end-1), 'linewidth', 6,'color', H2col)
    plot(x_H2O_plot(i,:), z(1:end-1), 'linewidth', 6,'color', H2Ocol)
    if H2only == false
        plot(x_CO_plot(i,:), z(1:end-1), 'linewidth', 6, 'color', COcol, 'LineStyle', ':')
        plot(x_CO2_plot(i,:), z(1:end-1), 'linewidth', 6,'color', CO2col, 'LineStyle', ':')
        plot(x_CH4_plot(i,:), z(1:end-1), 'linewidth', 6,'color', CH4col, 'LineStyle', ':')
        xlim([0, 0.8])

    end
   
    xlabel('Mole Fraction')
    %ylabel('Furnace Height (m)')
    
    ylim([0, h_furnace])
    xticks([0:0.25:1]);

    legend('H_2','H_2O', 'CO', 'CO_2', 'CH_4', 'Location', 'east')
    H = gca;
    grid on
    H.LineWidth = 3; %change to the desired value   
    set(gca,'FontWeight', 'bold','FontSize',18)

    subplot(1,3,3)
    box on
    hold on
    plot(T_g_plot(i,:), z(1:end-1), 'linewidth', 6,'color', gascol)
    plot(T_s_plot(i,:), z(2:end), 'linewidth', 6,'color', solidscol)
    
    xlabel('Temperature (^oC)')
    %ylabel('Furnace Height (m)')
    
    xlim([-10, 1000]);
    ylim([0, h_furnace])
    
    xticks(0:250:1000);
    xtickangle(0)
    
    legend('T_g','T_s', 'Location', 'south')
    H = gca;
    grid on
    H.LineWidth = 3; %change to the desired value   
    set(gca,'FontWeight', 'bold','FontSize',18)
    
   
    %text(0.1, 1, num2str(time_interp(i)/3600, '%.2f' ) + "hours", 'FontSize',14,'FontWeight', 'bold')
    time_interp(time_interp < 0) = 0 ;
    sgtitle(num2str(time_interp(i)/3600-1+.005, '%.2f' ) + " hours", 'FontSize',24,'FontWeight', 'bold')
    
    gif('DelayTime',1/10)

    if i < length(time_interp)-1
        clf
    end
end

close all