% 清空工作区、关闭所有图形窗口、清除命令行窗口内容
clear;
close all;
clc;

% 加载Indian Pines数据集（这里假设数据集文件名为'Indian_pines_corrected.mat'，可根据实际情况修改）
load('Indian_pines_corrected.mat'); % 加载高光谱图像数据

% 提取高光谱数据
newpic = indian_pines_corrected;

% 获取图像的行数、列数和波段数
[Lines, Columns, Bands] = size(newpic);
N = Lines * Columns;    % 像元数

% 将图像数据重塑为N行Bands列的矩阵，每个像元为一列
iData = reshape(newpic, N, Bands)';


% 初始化端元数（可根据实际情况调整）
numEndmember = 4;  

% 波段降维操作
% 计算波段间的相关系数矩阵
corrMatrix = corrcoef(iData');
% 用于标记要保留的波段，初始化为全部保留
keepBands = true(1, Bands);  
% 设定相关性阈值，可根据实际情况调整，这里设为0.95
corrThreshold = 0.95;  
for i = 1:Bands - 1
    for j = i + 1:Bands
        if corrMatrix(i, j) >= corrThreshold
            % 如果相关性高于阈值，随机选择一个波段去除（也可根据其他规则选择）
            if rand() < 0.5
                keepBands(j) = false;
            else
                keepBands(i) = false;
            end
        end
    end
end
% 提取降维后的波段数据
reducedIData = iData(:, keepBands);

% 基于降维后的数据进行VCA端元提取
% 先获取降维后数据的维度信息
[reducedBands, numReducedBands] = size(reducedIData);
% 用于存储提取的端元矩阵（在降维后的波段维度下）
E_reduced = zeros(reducedBands, numEndmember);

% 寻找第一个端元
% 计算每个像元在所有波段上的投影长度平方和
projLengthSquared = sum(reducedIData.^2, 1);
% 找到投影长度平方和最大的像元的索引
[~, index] = max(projLengthSquared);
E_reduced(:, 1) = reducedIData(:, index);

% 循环寻找后续端元
for k = 2:numEndmember
    % 构建投影矩阵，基于已找到的前k - 1个端元
    projMatrix = E_reduced(:, 1:k - 1) / (E_reduced(:, 1:k - 1)' * E_reduced(:, 1:k - 1)) * E_reduced(:, 1:k - 1)';
    % 添加一个很小的正则化项避免数值不稳定（分母接近0的情况），这里正则化参数设为1e-10，可根据实际调整
    epsilon = 1e-10;  
    projMatrix = E_reduced(:, 1:k - 1) / (E_reduced(:, 1:k - 1)' * E_reduced(:, 1:k - 1) + epsilon * eye(k - 1)) * E_reduced(:, 1:k - 1)';
    % 计算每个像元投影到与已找到端元正交空间后的投影长度平方和
    projOrthogonalLengthSquared = sum((eye(reducedBands) - projMatrix) * reducedIData.^2, 1);
    % 找到投影长度平方和最大的像元的索引
    [~, index] = max(projOrthogonalLengthSquared);
    E_reduced(:, k) = reducedIData(:, index);
end

% 将提取的端元对应回原始波段维度（补充被去除波段处的值为0）
E = zeros(Bands, numEndmember);
colIndex = 1;
for i = 1:Bands
    if keepBands(i)
        E(i, :) = E_reduced(:, colIndex);
        colIndex = colIndex + 1;
    end
end

% 显示端元的1~5波段值，并添加注释说明
disp('以下是提取的端元矩阵E的前5个波段对应的数据：');
E(1:5,:)

% 绘制提取的端元光谱，添加图形标题、坐标轴标签等使图形更直观
figure;
plot(E);
title('提取的端元光谱');
xlabel('波段序号');
ylabel('光谱强度');