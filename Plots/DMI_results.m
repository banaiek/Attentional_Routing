%% Script to generate DMI results from source data files

clear
close all

%% Load source data
load('source_DMI_data.mat');
fprintf('Data loaded successfully!\n');


datasetNames = {'CP15','IR13','IR18','IR24','IR29','IR32','IR34','JH21','ST15','ST18','ST19','ST27'};

TOI1 = Peak_data.TOI1;
TOI2 = Peak_data.TOI2;

%% Figure 1: Individual Subject DMI - Task 1 (Subjects 1-7)
fprintf('Generating Task 1 DMI plots...\n');
figure('Name','Task 1 - Subjects 1-7')
ALL_Att_Task1 = [];
ALL_Tar_Task1 = [];

for iSubj = 1:7
    Att = find(Att_binary{iSubj});
    Tar = find(Tar_binary{iSubj});
    
    if ~isempty(Att) && ~isempty(Tar)
        % Att→Tar
        subplot(2,7,iSubj)
        Att_Tar = squeeze(nanmean(MI_In{iSubj,1}(Att,Tar,:,:),[1,2]));
        imagesc(window_time, time_lags, Att_Tar)
        hold on
        plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
        plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
        linePlt = nanmean(Att_Tar(:,TOI1),2);
        linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
        plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
        colormap turbo
        set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025],'YTick',-500:250:500)
        title(datasetNames{iSubj})
        if iSubj == 1
            ylabel('Time lag (ms)')
        end
        axis xy
        % Tar→Att
        subplot(2,7,iSubj+7)
        Tar_Att = squeeze(nanmean(MI_In{iSubj,1}(Tar,Att,:,:),[1,2]));
        imagesc(window_time, time_lags, Tar_Att)
        hold on
        plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
        plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
        linePlt = nanmean(Tar_Att(:,TOI2),2);
        linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
        plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
        colormap turbo
        set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025],'YTick',-500:250:500)
        if iSubj == 1
            ylabel('Time lag (ms)')
        end
        xlabel('Time to target (ms)')
        axis xy
        ALL_Att_Task1(iSubj,:,:) = squeeze(nanmean(MI_In{iSubj,1}(Att,Tar,:,:),[1,2]));
        ALL_Tar_Task1(iSubj,:,:) = squeeze(nanmean(MI_In{iSubj,1}(Tar,Att,:,:),[1,2]));
    end
end
set(gcf,'Position',[28 464 1350 304])

%% Figure 2: Individual Subject DMI - Task 2 (Subjects 8-12)
fprintf('Generating Task 2 DMI plots...\n');
figure('Name','Task 2 - Subjects 8-12')
ALL_Att_Task2 = [];
ALL_Tar_Task2 = [];

for idx = 1:5
    iSubj = idx + 7;
    Att = find(Att_binary{iSubj});
    Tar = find(Tar_binary{iSubj});
    
    if ~isempty(Att) && ~isempty(Tar)
        % Att→Tar
        subplot(2,5,idx)
        Att_Tar = squeeze(nanmean(MI_In{iSubj,1}(Att,Tar,:,:),[1,2]));
        imagesc(window_time, time_lags, Att_Tar)
        hold on
        plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
        plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
        linePlt = nanmean(Att_Tar(:,TOI1),2);
        linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
        plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
        colormap turbo
        set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025],'YTick',-500:250:500)
        title(datasetNames{iSubj})
        if idx == 1
            ylabel('Time lag (ms)')
        end
        axis xy
        % Tar→Att
        subplot(2,5,idx+5)
        Tar_Att = squeeze(nanmean(MI_In{iSubj,1}(Tar,Att,:,:),[1,2]));
        imagesc(window_time, time_lags, Tar_Att)
        hold on
        plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
        plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
        linePlt = nanmean(Tar_Att(:,TOI2),2);
        linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
        plot(linePlt*500, time_lags,'LineWidth',1,'Color','w','LineStyle','--')
        colormap turbo
        set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025],'YTick',-500:250:500)
        if idx == 1
            ylabel('Time lag (ms)')
        end
        xlabel('Time to target (ms)')
        axis xy
        ALL_Att_Task2(idx,:,:) = squeeze(nanmean(MI_In{iSubj,1}(Att,Tar,:,:),[1,2]));
        ALL_Tar_Task2(idx,:,:) = squeeze(nanmean(MI_In{iSubj,1}(Tar,Att,:,:),[1,2]));
    end
end
set(gcf,'Position',[28 464 1350 304])

%% Figure 3: Mean DMI - Task 1 vs Task 2
fprintf('Generating mean DMI comparison plots...\n');
figure('Name','Mean DMI - Task 1 vs Task 2')

% Task 1: Att→Tar mean
subplot(2,2,3)
Tar_Att_Mn_Task1 = squeeze(nanmean(ALL_Att_Task1./max(ALL_Att_Task1,[],[2,3])));
imagesc(window_time, time_lags, Tar_Att_Mn_Task1)
hold on
plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
linePlt = nanmean(Tar_Att_Mn_Task1(:,TOI1),2);
linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
colormap turbo
set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025])
title('Task 1: Att→Tar')
ylabel('Time lag (ms)')
axis xy
% Task 1: Tar→Att mean
subplot(2,2,4)
Att_Tar_Mn_Task1 = squeeze(nanmean(ALL_Tar_Task1./max(ALL_Tar_Task1,[],[2,3])));
imagesc(window_time, time_lags, Att_Tar_Mn_Task1)
hold on
plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
linePlt = nanmean(Att_Tar_Mn_Task1(:,TOI2),2);
linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
colormap turbo
set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025])
title('Task 1: Tar→Att')
axis xy
% Task 2: Att→Tar mean
subplot(2,2,1)
if ~isempty(ALL_Att_Task2)
    Tar_Att_Mn_Task2 = squeeze(nanmean(ALL_Att_Task2./max(ALL_Att_Task2,[],[2,3])));
    imagesc(window_time, time_lags, Tar_Att_Mn_Task2)
    hold on
    plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
    plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
    linePlt = nanmean(Tar_Att_Mn_Task2(:,TOI1),2);
    linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
    plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
    colormap turbo
    set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025])
end
title('Task 2: Att→Tar')
ylabel('Time lag (ms)')
xlabel('Time to target (ms)')
axis xy
% Task 2: Tar→Att mean
subplot(2,2,2)
if ~isempty(ALL_Tar_Task2)
    Att_Tar_Mn_Task2 = squeeze(nanmean(ALL_Tar_Task2./max(ALL_Tar_Task2,[],[2,3])));
    imagesc(window_time, time_lags, Att_Tar_Mn_Task2)
    hold on
    plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','w','LineStyle','--')
    plot([0 0],get(gca,'ylim'),'LineWidth',1,'Color','w','LineStyle','--')
    linePlt = nanmean(Att_Tar_Mn_Task2(:,TOI2),2);
    linePlt = (linePlt-min(linePlt))./(max(linePlt)-min(linePlt));
    plot(linePlt*500, time_lags,'LineWidth',1,'Color','w')
    colormap turbo
    set(gca,'TickDir','out','xlim',[-750 750],'TickLength',[0.025 0.025])
end
title('Task 2: Tar→Att')
xlabel('Time to target (ms)')
axis xy
set(gcf,'Position',[680 601 1000 600])

%% Figure 4: Peak Timing Analysis
fprintf('Generating peak timing analysis...\n');
figure('Name','Delay MI Analysis - Peak Timing')

subplot(1,2,1)
hold on
plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','k','LineStyle','--')

for iSubj = 1:7
    if ~isnan(Peak_data.AttToTar_lag(iSubj))
        % Att→Tar
        plot(Peak_data.AttToTar_time(iSubj), Peak_data.AttToTar_lag(iSubj), ...
             'o','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b')
        
        % Tar→Att
        plot(Peak_data.TarToAtt_time(iSubj), Peak_data.TarToAtt_lag(iSubj), ...
             'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r')
    end
end

box off
set(gca,'TickDir','out','xlim',[-500 500],'ylim',[-250 250],...
    'TickLength',[0.025 0.025],'YTick',-250:125:250)
xlabel('Window Time (ms)')
ylabel('Time Lag (ms)')
title('Task 1 (Subjects 1-7)')
legend({'','Att→Tar','Tar→Att'},'Location','best')

subplot(1,2,2)
hold on
plot(get(gca,'xlim'),[0 0],'LineWidth',1,'Color','k','LineStyle','--')

for iSubj = 8:12
    if ~isnan(Peak_data.AttToTar_lag(iSubj))
        % Att→Tar
        plot(Peak_data.AttToTar_time(iSubj), Peak_data.AttToTar_lag(iSubj), ...
             'o','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b')
        
        % Tar→Att
        plot(Peak_data.TarToAtt_time(iSubj), Peak_data.TarToAtt_lag(iSubj), ...
             'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r')
    end
end

box off
set(gca,'TickDir','out','xlim',[-500 500],'ylim',[-250 250],...
    'TickLength',[0.025 0.025],'YTick',-250:125:250)
xlabel('Window Time (ms)')
ylabel('Time Lag (ms)')
title('Task 2 (Subjects 8-12)')
legend({'','Att→Tar','Tar→Att'},'Location','best')

set(gcf,'Position',[680 400 1200 400])

