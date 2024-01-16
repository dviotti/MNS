close all

deg2rad = pi/180;
rad2deg = 180/pi;

g0 = 9.780327;            % [m/s2]

rad2arcsec = 180/pi*3600;
mps22mg = 1e3/g0;
radps2degph = 180/pi*3600;

% sample = 1;
type = 1;
switch type
    case 1 % all data
        idx = (1:length(time))';
    case 2 % update
        idx = ~isnan(kU(:,end));
end

TimeFormat = 1;
switch TimeFormat
    case 1
        TS = 1;
        TF = '($s$)';
    case 2
        TS = 60;
        TF = '($min$)';
end
time_sim = time(idx)/TS;

ErrorType = 2;
switch ErrorType
    case 1
        x_sim = x;
    case 2
        x_sim = x_til;
end

[nKF,N] = size(x_sim);

P_sim = NaN(nKF,N);
for i=1:N
    P_sim(:,i) = sqrt(diag(P(:,:,i)));
end

nanvec = NaN(1,N);
if nKF==17
    x_sim = [x_sim; nanvec];
    P_sim = [P_sim; nanvec];
end

TH = 1.96; % 95% double-sided chi-squared hypothesis test
% TH = 2.58; % 99% double-sided chi-squared hypothesis test 

%% Plot

labely = {'$\widetilde{\Delta R}_N$ ($m$)','$\widetilde{\Delta R}_E$ ($m$)','$\widetilde{\Delta R}_D$ ($m$)',...
    '$\widetilde{\Delta V}_N$ ($m/s$)','$\widetilde{\Delta V}_E$ ($m/s$)','$\widetilde{\Delta V}_D$ ($m/s$)',...
    '$\widetilde{\Psi}_N$ ($arcmin$)','$\widetilde{\Psi}_E$ ($arcmin$)','$\widetilde{\Psi}_D$ ($arcmin$)',...
    '$\widetilde{\nabla}_{b,x}$ ($mg$)','$\widetilde{\nabla}_{b,y}$ ($mg$)','$\widetilde{\nabla}_{b,z}$ ($mg$)',...
    '$\widetilde{\varepsilon}_{b,x}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,y}$ ($deg/h$)','$\widetilde{\varepsilon}_{b,z}$ ($deg/h$)',...
    '$\widetilde{b}_{GPS}$ ($m$)','$\widetilde{\dot{b}}_{GPS}$ ($m/s$)','$\widetilde{b}_{ALT}$'};
convU = [1 1 1 1 1 1 rad2deg*60 rad2deg*60 rad2deg*60 ...
    mps22mg mps22mg mps22mg radps2degph radps2degph radps2degph 1 1 1];

% DR, DV, and Psi
figure
for j=1:9
    subplot(3,3,j)
    hold on
    stairs(time_sim,x_sim(j,idx)*convU(j),'LineWidth',2)
    stairs(time_sim,TH*P_sim(j,idx)*convU(j),'LineWidth',2,'Color','#A2142F')
    stairs(time_sim,-TH*P_sim(j,idx)*convU(j),'LineWidth',2,'Color','#A2142F')
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{j},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])

%%
% nabla, epsilon, and bias
figure
for j=1:9
    subplot(3,3,j)
    jj = 9 + j;
    hold on
    stairs(time_sim,x_sim(jj,idx)*convU(jj),'LineWidth',2)
    stairs(time_sim,TH*P_sim(jj,idx)*convU(jj),'LineWidth',2,'Color','#A2142F')
    stairs(time_sim,-TH*P_sim(jj,idx)*convU(jj),'LineWidth',2,'Color','#A2142F')
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    xlabel(['Time ' TF],'FontSize',14,'interpreter','latex')
    ylabel(labely{jj},'FontSize',14,'interpreter','latex')
    grid on
    set(gca,'FontSize',14);
    set(gca,'TickLabelInterpreter','latex');
end

set(gcf(), 'Units', 'normalized');
set(gcf(), 'Position', [0.15 0.15 0.6 0.75])

