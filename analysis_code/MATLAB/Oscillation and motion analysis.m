clear
clc
close all;
javaaddpath 'D:\matlab 2018a\java\mij.jar'
javaaddpath 'D:\matlab 2018a\java\ij.jar'
MIJ.start;
%% import all the control movies
CtrlMov1 = mijread("DUP_30nm-48h-dc.lif - Series001 - 23nm-1.717s.avi");
save('DUP_150nm-dc-28h.lif - Series004_Processed001 - C=0-54nm-10s.mat');
% load('DUP_48h-dc-30nm.lif - Series007_Processed001 - C=0.mat');
Movies_ctrl = {CtrlMov1, CtrlMov1};
%%
%%%% need to import the Data_Contrl_Apr.mat for the following processing
%
% PodoCtrl_16pixels = {PodoCtrl1, PodoCtrl2, PodoCtrl3, PodoCtrl4};
% RoiCtrl_all = {RoiCtrl1, RoiCtrl2, RoiCtrl3, RoiCtrl10_2, PodoCtrlMean10_3, PodoCtrlMean10_4};
% Here stop, and run the ijm code manually, also put the enlarged image in front
Maxima = MIJ.getResultsTable;
Indexnum = 1;
Pixelratio = 23/1;                     % unit nm/pixel
dt = 1.717;                              % unit sKey_Movies_Ctrl = {};
i=1
% for i=1:size(Movies_ctrl,2)
%%
%img = mijread( '30k.tif' );
%roi=roipoly(img)
%cropped_img=imcrop(img, roi)
%%

RoiBound = MIJ.getRoi(1)
  % use the ROI as boundary to find all podosomes inside the ROI
    Maxima0 = Maxima(1:2,:);
    % RoiBound = RoiCtrl_all{i}; 
    index = inpolygon(Maxima0(1,:), Maxima0(2,:), RoiBound(1,:),RoiBound(2,:));
    
    Podo_dynamics = Maxima(3:end,index);

    Maxima = Maxima0(:,index);
    Maxima_coord = Maxima*Pixelratio;

%% 对log数据进行处理
logWindow = MIJ.getLog;
% 直接转换为字符向量
logText = char(logWindow);
% 分割处理
lines = strsplit(logText, '\n');
lines = lines(~cellfun('isempty', lines));
% 正则表达式提取数据
pattern = 'centroidX (-?\d+\.?\d*) centroidY (-?\d+\.?\d*) Slice (\d+) Oval (\d+)';
data = struct('Slice', {}, 'Oval', {}, 'X', {}, 'Y', {});
for i = 1:length(lines)
    line = lines{i};
    tokens = regexp(line, pattern, 'tokens');
    if ~isempty(tokens)
        data(end+1) = struct(...
            'Slice', str2double(tokens{1}{3}), ...
            'Oval', str2double(tokens{1}{4}), ...
            'X', str2double(tokens{1}{1}), ...
            'Y', str2double(tokens{1}{2}) ...
        );
    end
end
% 转换为排序表格
tableData = struct2table(data);
sortedTable = sortrows(tableData, {'Oval', 'Slice'});
time = (sortedTable.Slice - 1) * dt; % 计算实际时间
Xdis = sortedTable.X*Pixelratio;Ydis = sortedTable.Y*Pixelratio; % 计算物理位移
% 合并时间，位移到表格
sortedTable.Time = time;sortedTable.Xdis = Xdis;sortedTable.Ydis = Ydis;

%% 生成时间伪彩色轨迹图

% 生成颜色映射矩阵
cmap = colormap(jet); % 获取当前颜色映射矩阵（默认64行）
nColors = size(cmap, 1); % 颜色数量（通常为64）

% 计算时间归一化到颜色索引（1到nColors）
cmin = min(time);
cmax = max(time);
time_normalized = (time - cmin) / (cmax - cmin); % 归一化到 [0, 1]
color_indices = round(time_normalized * (nColors - 1) + 1); % 映射到 [1, nColors]

% 初始化图形
figure;
hold on;
colormap(jet);
h_colorbar = colorbar;
h_colorbar.Label.String = 'Time (s)';
caxis([cmin, cmax]); 
xlabel('X displacement');
ylabel('Y displacement');
title('podosome trajectory（时间伪彩）');

% 获取所有Oval编号
ovals = unique(sortedTable.Oval);
numOvals = length(ovals);
% 初始化 cell 数组 podotracks
podotracks = cell(numOvals, 1); % 每个元素对应一个 Oval 的轨迹
for i = 1:length(ovals)
    oval = ovals(i);
    idx = sortedTable.Oval == oval;
    
    % 提取当前Oval的轨迹数据（按时间排序）
    x = sortedTable.Xdis(idx);
    y = sortedTable.Ydis(idx);
    t = sortedTable.Time(idx);
    c_idx = color_indices(idx);
    % 保存为 [t, x, y] 矩阵
    podotracks{i} = [t, x/1000, y/1000];
    % 绘制轨迹线段（颜色随时间变化）
    for j = 1:length(x)-1
        % 计算线段颜色（基于起始点时间）
        color = cmap(c_idx(j), :); % 从颜色矩阵取对应行
        plot([x(j), x(j+1)], [y(j), y(j+1)], ...
            'Color', color, 'LineWidth', 0.5); % 示例颜色混合
    end
 
end
axis equal;

% 预存初始和最终坐标
startPoints = zeros(length(ovals), 2);endPoints = zeros(length(ovals), 2);   
% 遍历每个 oval
for i = 1:length(ovals)
    oval = ovals(i);
    % 筛选当前 oval 的数据
    ovalData = sortedTable(sortedTable.Oval == oval, :);
    % 提取初始位置（时间最早）
    startPoints(i, :) = [ovalData.Xdis(1), ovalData.Ydis(1)];
    % 提取最终位置（时间最晚）
    endPoints(i, :) = [ovalData.Xdis(end), ovalData.Ydis(end)];
end
% 绘制初始位置（蓝色圆圈）
scatter(startPoints(:,1), startPoints(:,2), 30, 'b', 'o', 'filled');
% 绘制最终位置（红色圆圈）
scatter(endPoints(:,1), endPoints(:,2), 30, 'r', 'o', 'filled');

%% 利用MSD analyzer
ma = msdanalyzer(2, 'um', 's');ma = ma.addAll(podotracks);ma.plotTracks;ma.labelPlotTracks;
%确认有没有漂移
ma = ma.computeDrift('velocity');
figure
ma.plotDrift
ma.labelPlotTracks
% 开始进行计算
ma = ma.computeMSD;figure
ma.plotMSD;
cla
ma.plotMeanMSD(gca, true)
% Estimating the mean diffusion coefficient
[fo, gof] = ma.fitMeanMSD(0.2);
plot(fo)
ma.labelPlotMSD;
legend off
%% Estimating the individual diffusion coefficient
ma = ma.fitMSD(0.2);
good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit))
%% velocity  
    v = ma.getVelocities;V = vertcat( v{:} );
hist(V(:, 2:end), 50) 
box off
xlabel([ 'Velocity ( um / s )' ])
ylabel('Frequency')
mean(V(:,2:end))

%% 提取速度分量
Vx = V(:,2);Vy = V(:,3);
%计算总速度nm/s
speed = sqrt(Vx.^2 + Vy.^2)*1000;
%绘制直方图
figure;
h = histogram(speed, 'BinWidth', 1);  % 可以调整BinWidth参数改变分箱宽度
xlabel('Speed (nm/s)');
ylabel('Frequency');
title('Total Speed Distribution');
grid on;
% 显示统计信息
fprintf('平均速度: %.2f nm/s\n', mean(speed));
fprintf('速度标准差: %.2f nm/s\n', std(speed));

% 拟合对数正态分布
pd_logn = fitdist(speed, 'Lognormal');log_speed = log(speed);
% 计算拟合优度R²
% 计算观测值的对数
y_obs = log_speed;
% 计算拟合值（对数尺度）
y_pred = pd_logn.mu + pd_logn.sigma * norminv((1:numel(speed))'/(numel(speed)+1));
% 计算总平方和
SST = sum((y_obs - mean(y_obs)).^2);
% 计算回归平方和
SSR = sum((y_pred - mean(y_obs)).^2);
% 计算R²
R2 = SSR / SST;

fprintf('=== PODOSOME速度分析 ===\n');
fprintf('对数正态分布参数:\n');
fprintf('μ = %.4f (对数均值)\n', pd_logn.mu);
fprintf('σ = %.4f (对数标准差)\n', pd_logn.sigma);
fprintf('拟合优度 R² = %.4f\n', R2);

fprintf('\n描述性统计:\n');
fprintf('原始速度均值 = %.4f nm/s\n', mean(speed));
fprintf('原始速度中位数 = %.4f nm/s\n', median(speed));
fprintf('对数速度均值 = %.4f\n', mean(log_speed));
fprintf('对数速度标准差 = %.4f\n', std(log_speed));
fprintf('偏度 = %.4f (正偏)\n', skewness(speed));
fprintf('峰度 = %.4f (重尾)\n', kurtosis(speed));

%% velocity autocorrelation
ma = ma.computeVCorr;
ma.plotMeanVCorr

% Motion type analysis through log-log fitting.
ma = ma.fitLogLogMSD(0.5);
mean(ma.loglogfit.alpha)

%%
r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;
R2LIMIT = 0.8;
% Remove bad fits
bad_fits = r2fits < R2LIMIT;
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas(bad_fits) = [];
% T-test
[htest, pval] = ttest(alphas, 1, 0.05, 'left');
if ~htest
    [htest, pval] = ttest(alphas, 1, 0.05);
end

% Prepare string
str = { [ '\alpha = ' sprintf('%.2f ± %.2f (mean ± std, N = %d)', mean(alphas), std(alphas), numel(alphas)) ] };

if htest
    str{2} = sprintf('Significantly below 1, with p = %.2g', pval);
else
    str{2} = sprintf('Not significantly differend from 1, with p = %.2g', pval);
end

figure
hist(alphas);
box off
xlabel('\alpha')
ylabel('#')
yl = ylim(gca);
xl = xlim(gca);
text(xl(2), yl(2)+2, str, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16)
title('\alpha values distribution', ...
    'FontSize', 20)
ylim([0 yl(2)+2])
%% 定义要删除的列索引直接删除指定列
cols_to_remove = [30,31];
Podo_dynamics(:, cols_to_remove) = [];
%%  数据平滑示意 Savitzky-Golay
window_size = 1;      % 调整窗口大小（需为奇数）
polynomial_order = 1;  % 多项式阶数
smoothed_data = smoothdata(Podo_dynamics, 1, 'sgolay', window_size);
figure;
subplot(2,1,1);
plot(Podo_dynamics(:,17), 'b');
title('Original Data');
subplot(2,1,2);
plot(smoothed_data(:,17), 'r');
title('Smoothed Data');

%% 显著度寻峰
MinPeakProminence1 = 25;MinPeakWidth1 = 16;
% 初始化存储变量
peakTimes = cell(size(Podo_dynamics, 2), 1);  % 各列峰值时间点
peakIntervals = cell(size(Podo_dynamics, 2), 1); % 相邻峰时间间隔
meanIntervals = zeros(size(Podo_dynamics, 2), 1); % 平均间隔
stdIntervals = zeros(size(Podo_dynamics, 2), 1);  % 间隔标准差
SD1list = zeros(size(Podo_dynamics, 2), 1);  % 数据的抖动
meanpeaks = zeros(size(Podo_dynamics, 2), 1);  % 波峰平均值
meanvalleys = zeros(size(Podo_dynamics, 2), 1);  % 波谷平均值
meansmoothed_data = zeros(size(Podo_dynamics, 2), 1);%数据时间平均
close all;
% 循环处理每一列
for col = 1:size(Podo_dynamics, 2)
    currentData = smoothed_data(:, col);
    currentData_neg = -smoothed_data(:, col);
    timeVector = (0:length(currentData)-1)*dt;  % 生成时间向量
    SD1 = std(currentData); SD1list(col,1) = SD1; 
    meansmoothed_data(col) = mean(currentData);
    % 检测峰值
    figure;
    findpeaks(currentData,timeVector,'MinPeakProminence',MinPeakProminence1,'MinPeakWidth',MinPeakWidth1,'Annotate','extents')
    [peaks, locs] = findpeaks(currentData,timeVector,'MinPeakProminence',MinPeakProminence1);
%   findpeaks(currentData_neg,timeVector,'MinPeakProminence',MinPeakProminence1,'Annotate','extents')
    [valleys_neg, valleys_locs] = findpeaks(currentData_neg,timeVector,'MinPeakProminence',MinPeakProminence1,'MinPeakWidth',MinPeakWidth1);
    valleys = -valleys_neg;
    hold on;
    plot(valleys_locs,valleys,'^','MarkerFaceColor', 'r','MarkerEdgeColor','r')
    % 记录峰值时间点
    peakTimes{col} = locs';

    % 计算相邻峰间隔
    if length(locs) >= 2
        intervals = diff(peakTimes{col});
        peakIntervals{col} = intervals;
        meanIntervals(col) = mean(intervals);
        meanpeaks(col) = mean(peaks);
        meanvalleys(col) = mean(valleys);
        stdIntervals(col) = std(intervals);
    else
        warning('Column %d: Insufficient peaks detected (found %d peaks)', col, length(locs));
        peakIntervals{col} = [];
        meanIntervals(col) = NaN;
        meanpeaks(col) = mean(peaks);
        meanvalleys(col) = mean(valleys);
        stdIntervals(col) = NaN;
    end
end

%% 如果根据SD滤去不需要数据
mask = (SD1list > 15) | (SD1list < 5);  % 逻辑3或4运算，满足任一条件则标记为删除
Podo_dynamics(:, mask) = [];  % 删除满足条件的列
SD1list = SD1list(~mask);  % 仅保留未删除的元素
smoothed_data(:, mask)= [];
%% 计算各足体振荡周期与幅度的平均值及方差
meanAveragePeriod = mean(meanIntervals);
STDAveragePeriod = std(meanIntervals);

%计算各足体对应的幅值
AverageAmplitude = (meanpeaks-meanvalleys)/2;
relativeAverageAmplitude = AverageAmplitude./meansmoothed_data;
MeanrelativeAverageAmplitude = mean(relativeAverageAmplitude);
stdrelativeAverageAmplitude = std(relativeAverageAmplitude);
% 输出结果
fprintf('Mean Average Period: %.4f seconds\n', meanAveragePeriod);
fprintf('Standard Deviation of Average Period: %.4f seconds^2\n', STDAveragePeriod);
fprintf('Mean Relative Amplitude: %.4f \n ', MeanrelativeAverageAmplitude);
fprintf('Standard Deviation of Amplitude: %.4f \n', stdrelativeAverageAmplitude);