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

TestNames = {'SIM_1021_1040_FP4_HQTG1_GPSTC_M200',...
    'SIM_1021_1048_FP4_HQTG1_GPSTC_ALT_M200',...
    'SIM_1021_1051_FP4_HQTG1_GPSTC_ALT_MAG_M200',...
    'SIM_1105_0149_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_M200',...
    'SIM_1104_1927_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M200',...
    'SIM_1105_2017_FP4_HQTG1_ALT_MAG_CAMLMF_LPS_M200'}; %------> for ANEES

% TestNames = {'SIM_1108_1220_FP4_HQTG1_GPSTC_M200',...
%     'SIM_1108_0227_FP4_HQTG1_GPSTC_ALT_M200',...
%     'SIM_1107_0242_FP4_HQTG1_GPSTC_ALT_MAG_M200',...
%     'SIM_1105_0149_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_M200',...
%     'SIM_1104_1927_FP4_HQTG1_GPSTC_ALT_MAG_CAMLMF_LPS_M200',...
%     'SIM_1105_2017_FP4_HQTG1_ALT_MAG_CAMLMF_LPS_M200'};

NC = numel(TestNames);

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end

q = 1.96; % 95% double-sided chi-squared hypothesis test
% q = 2.58; % 99% double-sided chi-squared hypothesis test 

Ex_plot  = cell(NC,1);
stddev_sta_plot  = cell(NC,1);

x_plot = cell(NC,1);
stddev_hat_plot = cell(NC,1);

ANEES_TR_plot = cell(NC,1);
limInf_ANEES_TR_plot = cell(NC,1);
limSup_ANEES_TR_plot = cell(NC,1);

ANEES_NAV_plot = cell(NC,1);
limInf_ANEES_NAV_plot = cell(NC,1);
limSup_ANEES_NAV_plot = cell(NC,1);

ANIS_plot = cell(NC,1);
limInf_ANIS_plot = cell(NC,1);
limSup_ANIS_plot = cell(NC,1);

MCreal = [200 200 200 200 200 200];

for m=1:NC
    
    eval(['load ' TestNames{m} '.mat'])  
    
    % M = param.config.MCreal;
    M = MCreal(m);
    [nKF,N] = size(x_til);

    stddev_sta = NaN(nKF,N);
    stddev_hat = NaN(nKF,N);
    for i=1:N
        stddev_sta(:,i) = sqrt(diag(Cov_x(:,:,i)));
        stddev_hat(:,i) = sqrt(diag(P(:,:,i)));
    end
    
    nanvec = NaN(1,N);
    if nKF==17
        E_x = [E_x; nanvec];
        x_til = [x_til; nanvec];
        stddev_sta = [stddev_sta; nanvec];
        stddev_hat = [stddev_hat; nanvec];
    elseif nKF==16
        E_x = [E_x(1:end-1,:); nanvec; nanvec; E_x(end,:)];
        x_til = [x_til(1:end-1,:); nanvec; nanvec; x_til(end,:)];
        stddev_sta = [stddev_sta(1:end-1,:); nanvec; nanvec; stddev_sta(end,:)];
        stddev_hat = [stddev_hat(1:end-1,:); nanvec; nanvec; stddev_hat(end,:)];
    end
    
    Ex_plot{m} = E_x;
    x_plot{m} = x_til;
    stddev_sta_plot{m} = stddev_sta;
    stddev_hat_plot{m} = stddev_hat;
    
    time_plot = time/TS;
    
    % NEES_NAV
    ANEES_TR_plot{m} = ANEES_TR;
    dof_ANEES_TR = nKF*M;
    limInf_ANEES_TR_plot{m} = 1/2*(-q+sqrt(2*dof_ANEES_TR-1))^2/M;
    limSup_ANEES_TR_plot{m} = 1/2*(q+sqrt(2*dof_ANEES_TR-1))^2/M;
    
    % NEES_NAV
    ANEES_NAV_plot{m} = ANEES_NAV;
    dof_ANEES_NAV = 9*M;
    limInf_ANEES_NAV_plot{m} = 1/2*(-q+sqrt(2*dof_ANEES_NAV-1))^2/M;
    limSup_ANEES_NAV_plot{m} = 1/2*(q+sqrt(2*dof_ANEES_NAV-1))^2/M;
    
    % NIS
    dof_NIS = sum(DOF,2);
    NIS_bar = dof_NIS/M;
    limInf_NIS = 1/2*(-q+sqrt(2*dof_NIS-1)).^2/M;
    limSup_NIS = 1/2*(q+sqrt(2*dof_NIS-1)).^2/M;
    
    ANIS_plot{m} = ANIS./NIS_bar;
    limInf_ANIS_plot{m} = limInf_NIS./NIS_bar;
    limSup_ANIS_plot{m} = limSup_NIS./NIS_bar;

%     ANIS_plot{m} = ANIS;
%     NIS_bar_plot{m} = dof_NIS;
%     limInf_ANIS_plot{m} = limInf_NIS;
%     limSup_ANIS_plot{m} = limSup_NIS;
    
%     eval(['clear ' TestNames{m} '.mat'])
    
    S = whos('-file',[TestNames{m} '.mat']);
    for k = 1:length(S)
        eval(['clear ' S(k).name])
    end
end

load c1.mat
load c2.mat
load c3.mat
load c4.mat
load c5.mat
load c6.mat
autoc = {autocorr_c1,autocorr_c2,autocorr_c3,autocorr_c4,autocorr_c5,autocorr_c6};

%%

colors = {'#000000','#0072BD','#EDB120','#7E2F8E','#D95319','#77AC30'};
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
%%

% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
    hold on
    for n=1:NC
        semilogy(time_plot,abs(Ex_plot{n}(j,:))*convU(j),'Color',colors{n})
%         semilogy(time_plot,abs(x_plot{n}(j,:))*convU(j),'LineStyle','--','Color',colors{n})
%         semilogy(time_plot,q*stddev_sta_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         semilogy(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),'LineStyle',':','Color',colors{n},'LineWidth',2)
%         if j==1 || j==2 || j==3
%             yticks([1e-2 1e-1 1 10])
%         end
%         if j==4 || j==5 || j==6
%             yticks([1e-2 1e-1 1])
%         end
%         if j==7 || j==8 || j==9
%             yticks([1 10 100])
%         end %---> For covariance
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'yscale','log')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    xline([20 50 176.85 456.85 677],'--k') %456.85
    xlim([0 750])
    xticks([0 250 500 750])
    ylim([1e-4 10])
    yticks([1e-4 1e-3 1e-2 1e-1 1]) %---> For mean
%     ylim(limy{j})
end
legend('1','2','3','4','5','6','FontSize',10,'NumColumns',2)

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4covComp1.eps
% print -depsc2 FP4meanComp1.eps

%%
% nabla, epsilon, and bias
figure
for j=1:9
    subplot(3,3,j)
    jj = 9 + j;
    hold on
    for n=1:NC
        semilogy(time_plot,abs(Ex_plot{n}(jj,:))*convU(jj),'Color',colors{n})
%         semilogy(time_plot,abs(x_plot{n}(jj,:))*convU(jj),'LineStyle','--','Color',colors{n})
%         semilogy(time_plot,q*stddev_sta_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
%         semilogy(time_plot,q*stddev_hat_plot{n}(jj,:)*convU(jj),'LineStyle',':','Color',colors{n},'LineWidth',2)
%         if j==1 || j==2
%             yticks([1e-1 1 10])
%             ylim([1e-1 10])
%         end
%         if j==3
%             yticks([1e-1 1 10])
%         end
%         if j==4 || j==5 || j==6
%             yticks([1e-1 1 10])
%             ylim([1e-1 10])
%         end
%         if j==9
%             yticks([1e-3 1e-2])
%         end %---> For covariance
%         if j==7 || j==8 || j==9
%             yticks([1 10 100])
%         end
    end
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{jj},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'yscale','log')
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    xline([20 50 176.85 456.85 677],'--k')
    xlim([0 750])
    xticks([0 250 500 750])
    ylim([1e-4 10])
    yticks([1e-4 1e-3 1e-2 1e-1 1])
    % ylim([1e-4 1e-2]) %--> mean b_ALT
%     ylim(limy{jj})
end
legend('1','2','3','4','5','6','FontSize',10,'NumColumns',2)

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])
% pause
% print -depsc2 FP4covComp2.eps
% print -depsc2 FP4meanComp2.eps

%% --------------------
% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
    hold on
    for n=1:NC
%         stairs(time_plot,Ex_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,x_plot{n}(j,:)*convU(j),'LineStyle','--','Color',colors{n},'LineWidth',2)
        stairs(time_plot,q*stddev_sta_plot{n}(j,:)*convU(j),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*stddev_hat_plot{n}(j,:)*convU(j),'LineStyle','--','Color',colors{n},'LineWidth',2)
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
        stairs(time_plot,q*stddev_sta_plot{n}(jj,:)*convU(jj),'Color',colors{n},'LineWidth',2)
%         stairs(time_plot,q*stddev_hat_plot{n}(jj,:)*convU(jj),'LineStyle',':','Color',colors{n},'LineWidth',2)
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
    ylabel(['ANEES ' int2str(j)],'FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(6,2,2*j)
    hold on
    plot(time_plot,ANIS_plot{j},'.','Color','#0072BD','LineWidth',2)
%     plot(time_plot,NIS_bar_plot{j},'k','LineWidth',2)
%     plot(time_plot,ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
%     plot(time_plot,1.202*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,0.811*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,NIS_bar_plot{j}*ones(size(time_plot)),'k','LineWidth',2)
%     plot(time_plot,limInf_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,limSup_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(['$ANIS$ ' int2str(j)],'FontSize',14,'interpreter','latex')
    xlim([0 750])
%     ylim([0.5 1.5])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.1 0.1 0.5 1.1])
% pause
% print -depsc2 FP4_ANEESANIS.eps

%% --------------------
% ANEES_NAV autocorr
figure
for j=1:6
    subplot(6,2,2*j-1)
    hold on
    plot(time_plot,ANEES_NAV_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,9*ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(['ANEES ' int2str(j)],'FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(6,2,2*j)
    hold on
    stem(autoc{j},'filled','LineWidth',1)
    xlabel('Lag','FontSize',14,'interpreter','latex')
    ylabel(['Autocorr ' int2str(j)],'FontSize',14,'interpreter','latex')
%     xlim([0 750])
%     ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.1 0.1 0.5 1.1])
% pause
% print -depsc2 FP4_ANEESANIS.eps


%% --------------------
% ANEES_NAV ANIS AUTOCORR
figure
for j=1:6
    subplot(6,3,3*j-2)
    hold on
    plot(time_plot,ANEES_NAV_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,9*ones(size(time_plot)),'k','LineWidth',2)
    plot(time_plot,limInf_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,limSup_ANEES_NAV_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(['ANEES ' int2str(j)],'FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(6,3,3*j-1)
    hold on
    plot(time_plot,ANIS_plot{j},'.','Color','#0072BD','LineWidth',2)
%     plot(time_plot,NIS_bar_plot{j},'k','LineWidth',2)
    plot(time_plot,ones(size(time_plot)),'k','LineWidth',2)
%     plot(time_plot,limInf_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
%     plot(time_plot,limSup_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    plot(time_plot,1.202*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,0.761*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,NIS_bar_plot{j}*ones(size(time_plot)),'k','LineWidth',2)
%     plot(time_plot,limInf_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,limSup_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(['$ANIS$ ' int2str(j)],'FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([0.5 1.5])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(6,3,3*j)
    hold on
    stem(autoc{j},'LineWidth',1.5)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(['Autocorr ' int2str(j)],'FontSize',14,'interpreter','latex')
%     xlim([0 750])
%     ylim([6 12])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.1 0.1 0.7 1.2])
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
    xlim([0 750])
    ylim([15 20])
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
    subplot(6,3,3*j)
    hold on
    plot(time_plot,ANIS_plot{j},'.','Color','#0072BD','LineWidth',2)
    plot(time_plot,ones(size(time_plot)),'k','LineWidth',2)
%     plot(time_plot,limInf_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
%     plot(time_plot,limSup_ANIS_plot{j},'.','Color','#A2142F','LineWidth',2)
    plot(time_plot,1.202*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    plot(time_plot,0.811*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    % plot(time_plot,NIS_bar_plot{j}*ones(size(time_plot)),'k','LineWidth',2)
%     plot(time_plot,limInf_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
%     plot(time_plot,limSup_ANIS_plot{j}*ones(size(time_plot)),'Color','#A2142F','LineWidth',2)
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel('$ANIS$','FontSize',14,'interpreter','latex')
    xlim([0 750])
    ylim([0.5 1.5])
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.2 0.05 0.5 1])
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
