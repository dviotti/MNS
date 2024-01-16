close all
clear all
clc

g0 = 9.780327;            % [m/s2]

rad2deg = 180/pi;
rad2arcsec = 180/pi*3600;
mps22mg = 1e3/g0;
degph2radps=pi/180/3600;
radps2degph = 1/degph2radps;

ft2m = 0.3048;
m2ft = 1/ft2m;

TestNames = {'SIM_1031_1356_FTE2_HQTG1_GPSTC_M1',...
    'SIM_1031_1410_FTE2_HQTG1_GPSTC_HDG_M1',...
    'SIM_1031_1431_FTE2_HQTG1_GPSTC_HDG_ALT_M1',...
    'SIM_1031_1449_FTE2_HQTG1_GPSTC_HDG_ALT_LPS_M1'};

NC = numel(TestNames);

TimeFormat = 2;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end

stddev_hat_plot = cell(NC,1);

for m=1:NC
    
    eval(['load ' TestNames{m} '.mat'])  
    
    [nKF,N] = size(x_til);

    stddev_hat = NaN(nKF,N);
    for i=1:N
        stddev_hat(:,i) = sqrt(diag(P(:,:,i)));
    end
    
    nanvec = NaN(1,N);
    if nKF==17
        stddev_hat = [stddev_hat; nanvec];
    elseif nKF==16
        stddev_hat = [stddev_hat(1:end-1,:); nanvec; nanvec; stddev_hat(end,:)];
    end
    
    stddev_hat_plot{m} = stddev_hat;
    
    time_plot = time/TS;
    

%     eval(['clear ' TestNames{m} '.mat'])
    
    S = whos('-file',[TestNames{m} '.mat']);
    for k = 1:length(S)
        eval(['clear ' S(k).name])
    end
end

%%

q = 1.96; % 95% double-sided chi-squared hypothesis test
% q = 2.58; % 99% double-sided chi-squared hypothesis test 

% colors = {'#000000','#0072BD','#EDB120','#7E2F8E','#D95319','#77AC30'};
colors = {'#000000','#0072BD','#77AC30','#A2142F','#EDB120'};
labely = {'$\widetilde{\Delta R}_N$ ($m$)','$\widetilde{\Delta R}_E$ ($m$)','$\widetilde{\Delta R}_D$ ($m$)',...
    '$\widetilde{\Delta V}_N$ ($m/s$)','$\widetilde{\Delta V}_E$ ($m/s$)','$\widetilde{\Delta V}_D$ ($m/s$)',...
    '$\widetilde{\Psi}_N$ ($arcmin$)','$\widetilde{\Psi}_E$ ($arcmin$)','$\widetilde{\Psi}_D$ ($arcmin$)',...
    '$\widetilde{\nabla}_{b,x}$ ($mg$)','$\widetilde{\nabla}_{b,y}$ ($mg$)','$\widetilde{\nabla}_{b,z}$ ($mg$)',...
    '$\widetilde{\varepsilon}_{b,x}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,y}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,z}$ ($deg/h$)',...
    '$\widetilde{b}_{GPS}$ ($m$)','$\widetilde{\dot{b}}_{GPS}$ ($m/s$)','$\widetilde{b}_{ALT}$'};
convU = [1 1 1 1 1 1 rad2deg*60 rad2deg*60 rad2deg*60 mps22mg mps22mg mps22mg radps2degph radps2degph radps2degph 1 1 1];
% COV
% GPS / GPS_ALT
% limy = {[-10 10],[-10 10],[-10 10],[-1 1],[-1 1],[-0.5 0.5],[-10 10],[-10 10],[-75 75],...
%     [-5 5],[-5 5],[-5 5],[-10 10],[-10 10],[-10 10],[-10 10],[-1 1],[-0.025 0.025]};
% GPS_ALT_MAG / GPS_ALT_MAG_CAM
% limy = {[-10 10],[-10 10],[-10 10],[-1 1],[-1 1],[-0.5 0.5],[-10 10],[-10 10],[-10 10],...
%     [-2 2],[-2 2],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.02 0.02]};
% GPS_ALT_MAG_CAM_LPS
% limy = {[-2.5 2.5],[-2.5 2.5],[-2.5 2.5],[-0.25 0.25],[-0.25 0.25],[-0.25 0.25],[-5 5],[-5 5],[-5 5],...
%     [-1 1],[-1 1],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.01 0.01]};

% MEAN
% limy = {[1e-3 10],[1e-3 10],[1e-3 10],[1e-4 1],[1e-4 1],[1e-4 1],[1e-3 10],[1e-3 10],[1e-3 10],...
%     [-1 1],[-1 1],[-1 1],[-2 2],[-2 2],[-2 2],[-10 10],[-1 1],[-0.01 0.01]};

% style = {'-','--',':','-.'};
%%

% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
    hold on
    for n=1:NC
%         semilogy(time_plot,abs(Ex_plot{n}(j,:))*convU(j),'Color',colors{n})
%         semilogy(time_plot,abs(x_plot{n}(j,:))*convU(j),'LineStyle','--','Color',colors{n})
%         semilogy(time_plot,q*stddev_sta_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
        semilogy(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'yscale','log')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    if j==1 || j==2
        ylim([-inf 20])
        yticks([1e-1 1 5 10])
    end
    if j==3
        ylim([-inf 20])
        yticks([1e-1 1 10 20])
    end
    if j==4 || j==5 || j==6
        ylim([-inf 5])
        yticks([1e-1 1 5])
    end
    if j==7 || j==8
        yticks([10 100])
        ylim([10 inf])
    end
    if j==9
        ylim([10 2e3])
        yticks([10 1e2 1e3])
    end
    xlim([0 58])
    xticks([0 15 30 45 60])
    xline([8.5 15 26.5 47.5],'--k')
%     xlim([0 750])
%     xticks([0 250 500 750])
%     ylim([1e-2 10])
%     yticks([1e-4 1e-3 1e-2 1e-1 1])
%     ylim(limy{j})
end
legend('1','2','3','4','FontSize',10,'NumColumns',2)

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4covComp1.eps
% print -depsc2 FP4meanComp1.eps
% print -depsc2 FTcovComp1.eps

%%
% nabla, epsilon, and bias
figure
for j=1:9
    subplot(3,3,j)
    jj = 9 + j;
    hold on
    for n=1:NC
%         semilogy(time_plot,abs(Ex_plot{n}(jj,:))*convU(jj),'Color',colors{n})
%         semilogy(time_plot,abs(x_plot{n}(jj,:))*convU(jj),'LineStyle','--','Color',colors{n})
%         semilogy(time_plot,q*stddev_sta_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
        semilogy(time_plot,q*stddev_hat_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{jj},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'yscale','log')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    if j==3
        ylim([-inf 10])
%         yticks([1e-1 1 10])
    end
    if j==4 || j==5 || j==6
%         ylim([5 1e2])
        yticks([10 1e2])
    end
    if j==8
        ylim([-inf 1])
    end
    xlim([0 58])
    xticks([0 15 30 45 60])
    xline([8.5 15 26.5 47.5],'--k')
%     xlim([0 750])
%     xticks([0 250 500 750])
%     ylim([1e-4 10])
%     yticks([1e-4 1e-3 1e-2 1e-1 1])
%     ylim(limy{jj})
end
legend('1','2','3','4','FontSize',10,'NumColumns',2)

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4covComp2.eps
% print -depsc2 FP4meanComp2.eps
% print -depsc2 FTcovComp2.eps

%% --------------------
% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
    hold on
    for n=1:NC
%         stairs(time_plot,Ex_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,x_plot{n}(j,:)*convU(j),'LineStyle','--','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*stddev_sta_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
        stairs(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*stddev_sta_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*stddev_hat_plot{n}(j,:)*convU(j),'LineStyle','--','Color',colors{n},'LineWidth',2)
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4_GPS1.eps

% nabla, epsilon, and bias
figure
for j=1:9
    subplot(3,3,j)
    jj = 9 + j;
    hold on
    for n=1:NC
%         stairs(time_plot,Ex_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,x_plot{n}(jj,:)*convU(jj),'LineStyle','--','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*stddev_sta_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
        stairs(time_plot,q*stddev_hat_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*stddev_sta_plot{n}(jj,:)*convU(jj),'LineStyle','-.','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*stddev_hat_plot{n}(jj,:)*convU(jj),'LineStyle',':','Color',colors{n},'LineWidth',2)
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{jj},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4_GPS2.eps

%% --------------------
% ANEES_NAV ANIS
figure
for j=1:6
    subplot(6,2,2*j-1)
    hold on
    plot(time_plot,ANEES_NAV_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,9*ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANEES$','FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(6,2,2*j)
    hold on
    plot(time_plot,ANIS_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    % plot(time_plot,NIS_bar_plot{j}*ones(size(time_plot)),'k','LineWidth',2)
    % plot(time_plot,limInf_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    % plot(time_plot,limSup_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANIS$','FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([0 2])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.1 0.1 0.5 1])
% pause
% print -depsc2 FP4_ANEESANIS.eps

%% --------------------
% ANEES_TR ANEES_NAV ANIS
figure
for j=1:6
    % ANEES TR
    subplot(6,3,3*j-2)
    hold on
    plot(time_plot,ANEES_TR_plot{j},'Color','#0072BD','LineWidth',2)
    plot(time_plot,nKF*ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANEES_TR_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANEES_TR_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANEES \ TR$','FontSize',14,'interpreter','latex')
    % ylim([0 25])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    % ANEES NAV
    subplot(6,3,3*j-1)
    hold on
    plot(time_plot,ANEES_NAV_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,9*ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANEES$','FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    % ANIS
    subplot(6,2,3*j)
    hold on
    plot(time_plot,ANIS_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    % plot(time_plot,NIS_bar_plot{j}*ones(size(time_plot)),'k','LineWidth',2)
    % plot(time_plot,limInf_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    % plot(time_plot,limSup_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANIS$','FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([0 2])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.1 0.1 0.5 1])
% pause
% print -depsc2 FP4_GPS1.eps

%% --------------------
% % DR, DV, and Psi
% figure
% for j=1:9
%     subplot(3,3,j)
%     hold on
%     for n=1:NC
%         stairs(time_plot,Ex_plot{n}(j,:),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,x_plot{n}(j,:),'LineStyle','--','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*sqrt(Covx_plot{n}(j,:)),'LineStyle','-.','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*sqrt(P_plot{n}(j,:)),'LineStyle',':','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*sqrt(Covx_plot{n}(j,:)),'LineStyle','-.','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*sqrt(P_plot{n}(j,:)),'LineStyle',':','Color',colors{n},'LineWidth',2)
%     end
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel(labely{j},'FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end
% 
% % nabla, epsilon, and bias
% figure
% for j=1:9
%     subplot(3,3,j)
%     jj = 9 + j;
%     hold on
%     for n=1:NC
%         stairs(time_plot,Ex_plot{n}(jj,:),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,x_plot{n}(jj,:),'LineStyle','--','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*sqrt(Covx_plot{n}(jj,:)),'LineStyle','-.','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*sqrt(P_plot{n}(jj,:)),'LineStyle',':','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*sqrt(Covx_plot{n}(jj,:)),'LineStyle','-.','Color',colors{n},'LineWidth',2)
%         stairs(time_plot,-q*sqrt(P_plot{n}(jj,:)),'LineStyle',':','Color',colors{n},'LineWidth',2)
%     end
%     xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
%     ylabel(labely{jj},'FontSize',14,'interpreter','latex')
%     grid on
%     set(gca,'FontSize',14);
%     set(gca,'TickLabelInterpreter','latex');
% end
