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
time_interp = linspace(0,time(end), 250);

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


T_s_plot = zeros(length(time_interp), n_furnace) ;
T_g_plot = zeros(length(time_interp), n_furnace) ;

x_H2_plot = zeros(length(time_interp), n_furnace);
x_H2O_plot = zeros(length(time_interp), n_furnace);

w_Fe2O3_plot = zeros(length(time_interp), n_furnace) ;
w_Fe3O4_plot = zeros(length(time_interp), n_furnace) ;
w_FeO_plot = zeros(length(time_interp), n_furnace) ;
w_Fe_plot = zeros(length(time_interp), n_furnace) ;


for i =1:n_furnace
    properties = [x_H2(:,i), x_H2O(:,i), w_Fe2O3(:,i),  w_Fe3O4(:,i),  w_FeO(:,i), w_Fe(:,i), T_g(:,i), T_s(:,i)];
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

x_H2_plot =  [ones(size(time_interp'))*out.H2_concin.data x_H2_plot];
x_H2O_plot =  [ones(size(time_interp'))*out.H2O_concin.data x_H2O_plot];

% x_H2_plot =  [vq(:, 9) x_H2_plot];
% x_H2O_plot =  [vq(:,10) x_H2O_plot];

w_Fe2O3_plot =  [w_Fe2O3_plot ones(size(time_interp'))*out.w_Fe2O3in.data];
w_Fe3O4_plot =  [w_Fe3O4_plot ones(size(time_interp'))*out.w_Fe3O4in.data];
w_FeO_plot =  [w_FeO_plot ones(size(time_interp'))*out.w_FeOin.data];
w_Fe_plot =  [w_Fe_plot ones(size(time_interp'))*out.w_Fein.data];

T_g_plot = [ones(size(time_interp'))*(out.T_gin.data-273) T_g_plot];
T_s_plot = [T_s_plot ones(size(time_interp'))*(out.T_sin.data-273)];

% figure(1)
% figure('units','normalized','outerposition',[0 0 1 1])
% %subplot(1,2,1)
% xlim([0, 1])
% plot(w_Fe2O3(1,:), z, 'linewidth', 3, 'color','r')
% hold on
% plot(w_Fe3O4(1,:), z, 'linewidth', 3,'color', 'g')
% plot(w_FeO(1,:), z, 'linewidth', 3, 'color','b')
% plot(w_Fe(1,:), z, 'linewidth', 3,'color', 'k')
% set(gca,'FontWeight', 'bold','FontSize',18)


figure('units','normalized','outerposition',[0 0 1 1])
gif('DRI_step_change75nodes.gif')

%clf
%legend('Fe2O3','Fe3O4', 'FeO', 'Fe')

[minval,hour1] = min(abs(time_interp-3600));


for i = hour1:length(time_interp)-1
    
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
    xticks([0:0.25:1]);
    legend('Fe_2O_3','Fe_3O_4', 'FeO', 'Fe', 'Location', 'south')
    H = gca;
    grid on
    H.LineWidth = 3; %change to the desired value   
    set(gca,'FontWeight', 'bold','FontSize',18)
    

    subplot(1,3,2)
    box on
    hold on
    plot(x_H2_plot(i,:), z(1:end-1), 'linewidth', 6,'color', H2col)
    plot(x_H2O_plot(i,:), z(1:end-1), 'linewidth', 6,'color', H2Ocol)
    
    xlabel('Mole Fraction')
    ylabel('Furnace Height (m)')
    
    xlim([0, 1])
    ylim([0, h_furnace])
    xticks([0:0.25:1]);

    legend('H_2','H_2O', 'Location', 'south')
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
    ylabel('Furnace Height (m)')
    
    xlim([250, 850]);
    ylim([0, h_furnace])
    
    legend('T_g','T_s', 'Location', 'south')
    H = gca;
    grid on
    H.LineWidth = 3; %change to the desired value   
    set(gca,'FontWeight', 'bold','FontSize',18)
    
   
    %text(0.1, 1, num2str(time_interp(i)/3600, '%.2f' ) + "hours", 'FontSize',14,'FontWeight', 'bold')
    sgtitle(num2str(time_interp(i)/3600-1, '%.2f' ) + " hours", 'FontSize',24,'FontWeight', 'bold')
    
    gif('DelayTime',1/6)
    if i < length(time_interp)-1
        clf
    end
end




% for i = 1:length(time)-1
%     subplot(1,2,1)
%     xlim([-0.01, 1])
%     hold on
%     plot(w_Fe2O3(i,:), z, 'linewidth', 3,'color', 'r')
%     plot(w_Fe3O4(i,:), z, 'linewidth', 3,'color', 'g')
%     plot(w_FeO(i,:), z, 'linewidth', 3, 'color','b')
%     plot(w_Fe(i,:), z, 'linewidth', 3,'color', 'k')
% 
%     legend('Fe2O3','Fe3O4', 'FeO', 'Fe')
% 
%     subplot(1,2,2)
%     xlim([300, 900])
%     hold on
%     plot(T_s(i,:), z, 'linewidth', 3,'color', 'r')
%     plot(T_g(i,:), z, 'linewidth', 3,'color', 'g')
% 
%     legend('T_s', 'T_g')
% 
%     tpause = (time(i+1) - time(i))/3600;
%     text(0.1, 1, string(time(i)/3600))
% 
%     pause(tpause)
%     if i < length(time)-1
%         clf
%     end
% end
