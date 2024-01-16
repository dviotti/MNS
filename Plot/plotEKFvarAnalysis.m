close all

rad2arcsec = 180/pi*3600;
g0 = 9.780327;
mps22mg = 1e3/g0;
radps2degph = 180/pi*3600; 
time = 0:0.01:80;

sample = 1;
% t1 = 30;
% t2 = 60;
% t3 = 90;
% tend = 150;
t1 = 20;
t2 = 40;
t3 = 60;
tend = 80;

plotType = 2;

switch plotType
    
    case 0
        
        figure
        % DR
        subplot(3,3,1)
        stairs(time,sqrt(squeeze(Pu(1,1,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_N$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,2)
        stairs(time,sqrt(squeeze(Pu(2,2,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_E$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,3)
        stairs(time,sqrt(squeeze(Pu(3,3,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_D$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % DV
        subplot(3,3,4)
        stairs(time,sqrt(squeeze(Pu(4,4,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_N$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,5)
        stairs(time,sqrt(squeeze(Pu(5,5,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_E$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,6)
        stairs(time,sqrt(squeeze(Pu(6,6,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_D$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Psi
        subplot(3,3,7)
        stairs(time,sqrt(squeeze(Pu(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_x$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,8)
        stairs(time,sqrt(squeeze(Pu(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_y$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,9)
        stairs(time,sqrt(squeeze(Pu(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_z$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
%         print -depsc2 departure_gpsmagA.eps
%         print -dpng -r400 departure_camera.png
        
        figure
        % Nabla
        subplot(3,3,1)
        stairs(time,sqrt(squeeze(Pu(10,10,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_x$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,2)
        stairs(time,sqrt(squeeze(Pu(11,11,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_y$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,3)
        stairs(time,sqrt(squeeze(Pu(12,12,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_z$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Espilon
        subplot(3,3,4)
        stairs(time,sqrt(squeeze(Pu(13,13,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_x$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,5)
        stairs(time,sqrt(squeeze(Pu(14,14,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_y$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,6)
        stairs(time,2*sqrt(squeeze(Pu(15,15,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_z$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % GPS
        subplot(3,3,7)
        stairs(time,sqrt(squeeze(Pu(16,16,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        ylim([0 10])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');

        subplot(3,3,8)
        stairs(time,sqrt(squeeze(Pu(17,17,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        ylim([0 4])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
%       print -depsc2 departure_lps.eps
%       print -dpng -r400 departure_camera.png
    
    case 1
        
        figure
        % DR
        subplot(5,3,1)
        stairs(time,sqrt(squeeze(Pu(1,1,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_N$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,2)
        stairs(time,sqrt(squeeze(Pu(2,2,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_E$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,3)
        stairs(time,sqrt(squeeze(Pu(3,3,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_D$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % DV
        subplot(5,3,4)
        stairs(time,sqrt(squeeze(Pu(4,4,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_N$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,5)
        stairs(time,sqrt(squeeze(Pu(5,5,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_E$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,6)
        stairs(time,sqrt(squeeze(Pu(6,6,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_D$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Psi
        subplot(5,3,7)
        stairs(time,sqrt(squeeze(Pu(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_x$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,8)
        stairs(time,sqrt(squeeze(Pu(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_y$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,9)
        stairs(time,sqrt(squeeze(Pu(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_z$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % figure
        % Nabla
        subplot(5,3,10)
        stairs(time,sqrt(squeeze(Pu(10,10,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_x$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,11)
        stairs(time,sqrt(squeeze(Pu(11,11,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_y$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,12)
        stairs(time,sqrt(squeeze(Pu(12,12,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_z$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Espilon
        subplot(5,3,13)
        stairs(time,sqrt(squeeze(Pu(13,13,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_x$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,14)
        stairs(time,sqrt(squeeze(Pu(14,14,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_y$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,15)
        stairs(time,2*sqrt(squeeze(Pu(15,15,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_z$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % GPS
%         subplot(3,3,7)
%         stairs(time,sqrt(squeeze(Pu(16,16,:,sample))'),'LineWidth',2)
%         xline(t1,'k--')
%         xline(t2,'k--')
%         xline(t3,'k--')
%         xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
%         ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%         grid on
%         xlim([0 tend])
%         ylim([0 10])
%         set(gca,'FontSize',14);
%         set(gca,'TickLabelInterpreter','latex');
% 
%         subplot(3,3,8)
%         stairs(time,sqrt(squeeze(Pu(17,17,:,sample))'),'LineWidth',2)
%         xline(t1,'k--')
%         xline(t2,'k--')
%         xline(t3,'k--')
%         xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
%         ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
%         grid on
%         xlim([0 tend])
%         ylim([0 4])
%         set(gca,'FontSize',14);
%         set(gca,'TickLabelInterpreter','latex');
      
        pause  
      print -depsc2 departure_all.eps
%       print -dpng -r400 departure_camera.png
        
    case 2
        
        load Puall.mat
        load Puall-g.mat
        
        figure
        % DR
        subplot(5,3,1)
        hold on
        stairs(time,sqrt(squeeze(Pu1(1,1,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(1,1,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_N$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,2)
        hold on
        stairs(time,sqrt(squeeze(Pu1(2,2,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(2,2,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_E$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,3)
        hold on
        stairs(time,sqrt(squeeze(Pu1(3,3,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(3,3,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_D$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        legend({'ALL','ALL-GPS'},'FontSize',14,'interpreter','latex')
        
        % DV
        subplot(5,3,4)
        hold on
        stairs(time,sqrt(squeeze(Pu1(4,4,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(4,4,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_N$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,5)
        hold on
        stairs(time,sqrt(squeeze(Pu1(5,5,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(5,5,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_E$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,6)
        hold on
        stairs(time,sqrt(squeeze(Pu1(6,6,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(6,6,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_D$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Psi
        subplot(5,3,7)
        hold on
        stairs(time,sqrt(squeeze(Pu1(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_x$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,8)
        hold on
        stairs(time,sqrt(squeeze(Pu1(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_y$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,9)
        hold on
        stairs(time,sqrt(squeeze(Pu1(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_z$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
%         pause
%         print -depsc2 covPugPugm_RotA_S.eps
%         print -dpng -r400 covPugPugm_RotA_S.png
        
%         figure
        % Nabla
        subplot(5,3,10)
        hold on
        stairs(time,sqrt(squeeze(Pu1(10,10,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(10,10,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_x$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,11)
        hold on
        stairs(time,sqrt(squeeze(Pu1(11,11,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(11,11,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_y$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,12)
        hold on
        stairs(time,sqrt(squeeze(Pu1(12,12,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(12,12,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_z$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
%         legend({'GPS+MAG','GPS+MAG+CAMLM'},'FontSize',14,'interpreter','latex')
        
        % Espilon
        subplot(5,3,13)
        hold on
        stairs(time,sqrt(squeeze(Pu1(13,13,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(13,13,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_x$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,14)
        hold on
        stairs(time,sqrt(squeeze(Pu1(14,14,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(14,14,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_y$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(5,3,15)
        hold on
        stairs(time,2*sqrt(squeeze(Pu1(15,15,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(15,15,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_z$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        %GPS
%         subplot(3,3,7)
%         hold  on
%         stairs(time,sqrt(squeeze(Pug(16,16,:,sample))'),'LineWidth',2)
%         stairs(time,sqrt(squeeze(Pugm(16,16,:,sample))'),'LineWidth',2)
%         xline(t1,'k--')
%         xline(t2,'k--')
%         xline(t3,'k--')
%         xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
%         ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
%         grid on
%         xlim([0 tend])
%         ylim([0 10])
%         set(gca,'FontSize',14);
%         set(gca,'TickLabelInterpreter','latex');
% 
%         subplot(3,3,8)
%         hold on
%         stairs(time,sqrt(squeeze(Pug(17,17,:,sample))'),'LineWidth',2)
%         stairs(time,sqrt(squeeze(Pugm(17,17,:,sample))'),'LineWidth',2)
%         xline(t1,'k--')
%         xline(t2,'k--')
%         xline(t3,'k--')
%         xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
%         ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
%         grid on
%         xlim([0 tend])
%         ylim([0 4])
%         set(gca,'FontSize',14);
%         set(gca,'TickLabelInterpreter','latex');
        
        pause
        print -depsc2 departure_all.eps
%         print -dpng -r400 covPugPugm_RotB_S.png

    case 3
        
%         load Puc.mat
%         load Pugmc.mat
        
        figure
        % DR
        hold on
        subplot(3,3,1)
        stairs(time,sqrt(squeeze(Pu1(1,1,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(1,1,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_N$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,2)
        hold on
        stairs(time,sqrt(squeeze(Pu1(2,2,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(2,2,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_E$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,3)
        hold on
        stairs(time,sqrt(squeeze(Pu1(3,3,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(3,3,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta R_D$ (m)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        legend({'GPS','GPS+MAG'},'FontSize',14,'interpreter','latex')
        
        % DV
        subplot(3,3,4)
        hold on
        stairs(time,sqrt(squeeze(Pu1(4,4,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(4,4,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_N$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,5)
        hold on
        stairs(time,sqrt(squeeze(Pu1(5,5,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(5,5,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_E$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,6)
        hold on
        stairs(time,sqrt(squeeze(Pu1(6,6,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(6,6,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\Delta V_D$ (m/s)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        % Psi
        subplot(3,3,7)
        hold on
        stairs(time,sqrt(squeeze(Pu1(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(7,7,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_x$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,8)
        hold on
        stairs(time,sqrt(squeeze(Pu1(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(8,8,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_y$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,9)
        hold on
        stairs(time,sqrt(squeeze(Pu1(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(9,9,:,sample))')*rad2arcsec,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\psi_z$ (arcsec)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        pause
        print -depsc2 departure_gpsmagcamA.eps
%         print -dpng -r400 covPugPugm_RotA_S.png
        
        figure
        % Nabla
        subplot(3,3,1)
        hold on
        stairs(time,sqrt(squeeze(Pu1(10,10,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(10,10,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_x$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,2)
        hold on
        stairs(time,sqrt(squeeze(Pu1(11,11,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(11,11,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_y$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,3)
        hold on
        stairs(time,sqrt(squeeze(Pu1(12,12,:,sample))')*mps22mg,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(12,12,:,sample))')*mps22mg,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\nabla_z$ ($mg$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        legend({'GPS','GPS+MAG'},'FontSize',14,'interpreter','latex')
        
        % Espilon
        subplot(3,3,4)
        hold on
        stairs(time,sqrt(squeeze(Pu1(13,13,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(13,13,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_x$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,5)
        hold on
        stairs(time,sqrt(squeeze(Pu1(14,14,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(14,14,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_y$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        subplot(3,3,6)
        hold on
        stairs(time,2*sqrt(squeeze(Pu1(15,15,:,sample))')*radps2degph,'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(15,15,:,sample))')*radps2degph,'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\epsilon_z$ ($deg/h$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        %GPS
        subplot(3,3,7)
        hold  on
        stairs(time,sqrt(squeeze(Pu1(16,16,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(16,16,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$b_{GPS}$ ($m$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        ylim([0 10])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');

        subplot(3,3,8)
        hold on
        stairs(time,sqrt(squeeze(Pu1(17,17,:,sample))'),'LineWidth',2)
        stairs(time,sqrt(squeeze(Pu2(17,17,:,sample))'),'LineWidth',2)
        xline(t1,'k--')
        xline(t2,'k--')
        xline(t3,'k--')
        xlabel('Time ($s$)','FontSize',14,'interpreter','latex')
        ylabel('$\dot{b}_{GPS}$ ($m/s$)','FontSize',14,'interpreter','latex')
        grid on
        xlim([0 tend])
        ylim([0 4])
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','latex');
        
        pause
        print -depsc2 departure_gpsmagcamA.eps
%         print -dpng -r400 covPugPugm_RotB_S.png
        
end