clear; close all; clc;

load('model_source_data.mat');

a1 = source_data.threshold.a1;
a2 = source_data.threshold.a2;

target_accuracy = 0.85;
plot_colors = lines(2);

%% threshold analysis
stats = @(x) struct('mean',mean(x),'std',std(x),'n',numel(x));
s1 = stats(a1);
s2 = stats(a2);

fprintf('=== Condition summaries (signal present) ===\n');
fprintf('a1: mean = %.3f,  SD = %.3f,  n = %d\n', s1.mean, s1.std, s1.n);
fprintf('a2: mean = %.3f,  SD = %.3f,  n = %d\n\n', s2.mean, s2.std, s2.n);

all_obs = sort([a1; a2]);
u = unique(all_obs);
thr_cand = [-inf; (u(1:end-1)+u(2:end))/2; inf];

det_rate_a2 = arrayfun(@(t) mean(a2 >= t), thr_cand);
[~, best_idx] = min(abs(det_rate_a2 - target_accuracy));
thr_opt = thr_cand(best_idx);

hit_a2 = det_rate_a2(best_idx);
hit_a1 = mean(a1 >= thr_opt);

fprintf('Optimal threshold : %.5g\n', thr_opt);
fprintf('Hit rate in a2    : %.2f %% (target ≈ 85 %%)\n', hit_a2*100);
fprintf('Hit rate in a1    : %.2f %%\n\n', hit_a1*100);

obs = [a1; a2];
resp = obs >= thr_opt;

glm_b = glmfit(obs, resp, 'binomial', 'link','logit');
xx = linspace(min(obs), max(obs), 400)';
p_det = glmval(glm_b, xx, 'logit');

%% threshold figure
figure('Color','w','Position',[50 200 1200 450]);

subplot(1,3,1); hold on;
[f1,x1] = ksdensity(a1);
[f2,x2] = ksdensity(a2);

area(x1, f1, 'FaceColor', plot_colors(1,:), 'FaceAlpha', 0.35, 'EdgeColor','none');
area(x2, f2, 'FaceColor', plot_colors(2,:), 'FaceAlpha', 0.35, 'EdgeColor','none');

idx_hit_a2 = x2 >= thr_opt;
area(x2(idx_hit_a2), f2(idx_hit_a2), 'FaceColor', plot_colors(2,:), 'FaceAlpha', 0.70, 'EdgeColor','none');

xline(thr_opt, 'k-', 'LineWidth', 2);
xlabel('Observation value');
ylabel('Probability density');
legend({'a1 density','a2 density','a2 hits','threshold'},'Location','best');
title('Density distributions with decision threshold');
box off; set(gca,"TickDir","out");

subplot(1,3,2); hold on;
hit_a1_all = arrayfun(@(t) mean(a1>=t), thr_cand);
plot(thr_cand, hit_a1_all*100, 'Color', plot_colors(1,:), 'LineWidth',1.2);
plot(thr_cand, det_rate_a2*100, 'Color', plot_colors(2,:), 'LineWidth',1.2);
plot(thr_opt, hit_a1*100, 'o', 'MarkerEdgeColor', plot_colors(1,:), 'MarkerFaceColor', plot_colors(1,:));
plot(thr_opt, hit_a2*100, 'o', 'MarkerEdgeColor', plot_colors(2,:), 'MarkerFaceColor', plot_colors(2,:));
yline(target_accuracy*100, 'k:', 'LineWidth',1.5);
xlabel('Threshold');
ylabel('Detection rate (%)');
title('Detection rate vs. threshold');
legend({'a1','a2',sprintf('%.1f %% (a1)',hit_a1*100),sprintf('%.1f %% (a2)',hit_a2*100), 'Target 85 %'},'Location','best');
grid on; box off; set(gca,"TickDir","out");

subplot(1,3,3); hold on;
plot(xx, p_det, 'k-', 'LineWidth',2);
rng(0);
scatter(obs(resp==1), resp(resp==1)+rand(sum(resp==1),1)*0.02-0.01, 15, plot_colors(1,:), 'filled', 'MarkerFaceAlpha',0.45);
scatter(obs(resp==0), resp(resp==0)+rand(sum(resp==0),1)*0.02-0.01, 15, plot_colors(2,:), 'filled', 'MarkerFaceAlpha',0.45);
plot([thr_opt thr_opt],[0 1], 'k--', 'LineWidth',1);
ylim([-0.05 1.05]); xlim([min(obs) max(obs)]);
xlabel('Observation value');
ylabel('P(detect)');
title('Logistic-regression fit');
legend({'logit fit','Detected','Missed','threshold'},'Location','southeast');
box off;

sgtitle(sprintf('Threshold chosen so that a2 hit rate ≈ %.1f %% (target 85 %%)', hit_a2*100), 'FontWeight','bold');

%% time course figure
cValid = [0.4940 0.1840 0.5560];
cInvalid = [0.0000 0.5000 0.0000];

figure('Color','w');

subplot(1,2,1); hold on;
shadedErrorBar(source_data.timecourse.target.time_ms, ...
    source_data.timecourse.target.valid_mean, ...
    source_data.timecourse.target.valid_sem, ...
    {'-','Color',cValid,'LineWidth',1.5}, 1);
hI = shadedErrorBar(source_data.timecourse.target.time_ms, ...
    source_data.timecourse.target.invalid_mean, ...
    source_data.timecourse.target.invalid_sem, ...
    {'-','Color',cInvalid,'LineWidth',1.5}, 1);
xlim([-700 700]); xticks([-500 0 500]);
plot([-700 700], [0 0], 'k--');
set(gca,'TickDir','out','Box','off');
xlabel('Time (ms)'); ylabel('Z-scored response');
legend({'Valid','Invalid'}, 'Location','Best');
title('Target response');

subplot(1,2,2); hold on;
shadedErrorBar(source_data.timecourse.cue.time_ms, ...
    source_data.timecourse.cue.valid_mean, ...
    source_data.timecourse.cue.valid_sem, ...
    {'-','Color',cValid,'LineWidth',1.5}, 1);
shadedErrorBar(source_data.timecourse.cue.time_ms, ...
    source_data.timecourse.cue.invalid_mean, ...
    source_data.timecourse.cue.invalid_sem, ...
    {'-','Color',cInvalid,'LineWidth',1.5}, 1);
xlim([-700 700]); xticks([-500 0 500]);
plot([-700 700], [0 0], 'k--');
set(gca,'TickDir','out','Box','off');
xlabel('Time (ms)'); ylabel('Z-scored response');
legend({'Valid','Invalid'}, 'Location','Best');
title('Cue response');

sgtitle('Valid vs. Invalid');