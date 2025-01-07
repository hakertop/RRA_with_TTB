function WaterHead = Search_WaterHead(FlowAcc, FlowDir, LeftX, DownY, cellsize, ...
    ouletX, ouletY, Ac, outputFilePath)
% Search_WaterHead 函数计算水头位置并将其保存至文件 输入参数:
%   FlowAcc - 流量累积矩阵 
%   FlowDir - 流向矩阵 
%   LeftX, DownY - 区域的左下角坐标 
%   cellsize - 网格的大小 
%   ouletX, ouletY - 出口点的坐标向量 
%   Ac - 汇水区阈值 
%   outputFilePath -输出文本路径(要求完整，包括文件及其后缀)

% 初始化相关变量
format long; % 保持数据精度
[m, n] = size(FlowAcc); % 获取流量累积矩阵的行列数
UpY = DownY + cellsize * m; % 计算区域的上边界Y坐标
Ac = floor(Ac/(cellsize^2)); % 将汇水面积转换成单元格数目

len_ouletX = length(ouletX);
stepPoint = [];
WaterHead = [];

num_ToPoint = 0; % 汇入该点的栅格数
num_slope = 0; % 汇入该点的栅格中，<Ac的栅格数目
% 当二者相等，说明该点为水头

visited = false(m, n);

for ilen = 1:len_ouletX
    % 起始搜索点的行、列号
    startRow = m - floor((ouletY(ilen) - DownY)/cellsize);
    startCol = floor((ouletX(ilen) - LeftX)/cellsize) + 1;

    % 检查索引值有效性
    if startRow >= 1 && startRow <= m && startCol >= 1 && startCol <= n
        stepPoint = [stepPoint; [startRow, startCol, FlowAcc(startRow, startCol), 0, 0]];
    else
        warning('起始点(%d, %d)/（%d, %d）超出范围，跳过此点', startRow, startCol, ouletX(ilen), ouletY(ilen));
        warning('UpY: %d, oultetY(ilen): %d', UpY, ouletY(ilen));
    end
    empty_stepPoint = 1;

    % 搜索点行、列号；该点的汇水面积；该点汇入区域的行、列号（0,0表示无汇入）
    while empty_stepPoint
        startRow = stepPoint(1, 1);
        startCol = stepPoint(1, 2);

        if visited(startRow, startCol)
            stepPoint(1, :) = []; % % 每一轮搜索完成后，原来的起始点 出栈, 第一个点为上一步的搜索点
            empty_stepPoint = ~isempty(stepPoint);
            continue;
        end
        visited(startRow, startCol) = true;

        if FlowAcc(startRow, startCol) >= Ac
            % 确保该点是河流点
            if startCol < n && FlowDir(startRow, startCol+1) == 2^4 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                % 该点右侧(>Ac)，入栈
                if FlowAcc(startRow, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow, startCol + 1, FlowAcc(startRow, startCol+1), startRow, startCol]];
                    % 该点右侧(<Ac)，该点有作为水头的可能
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && startCol < n && FlowDir(startRow+1, startCol+1) == 2^5 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow+1, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol + 1, FlowAcc(startRow+1, startCol+1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && FlowDir(startRow+1, startCol) == 2^6 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow+1, startCol) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol, FlowAcc(startRow+1, startCol), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && startCol > 1 && FlowDir(startRow+1, startCol-1) == 2^7 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow+1,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow+1, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol - 1, FlowAcc(startRow+1, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startCol > 1 && FlowDir(startRow, startCol-1) == 2^0 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow, startCol - 1, FlowAcc(startRow, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && startCol > 1 && FlowDir(startRow-1, startCol-1) == 2^1 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow-1, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol - 1, FlowAcc(startRow-1, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && FlowDir(startRow-1, startCol) == 2^2 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow-1, startCol) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol, FlowAcc(startRow-1, startCol), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && startCol < n && FlowDir(startRow-1, startCol+1) == 2^3 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow-1,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % 汇入该点
                if FlowAcc(startRow-1, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol + 1, FlowAcc(startRow-1, startCol+1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            % 汇入点数目 与 汇入点<Ac数目 相等，则该点为水头
            if num_slope == num_ToPoint
                WaterHead = [WaterHead; [startRow, startCol]];
            end
        end
        num_slope = 0;
        num_ToPoint = 0;
        stepPoint(1, :) = []; % % 每一轮搜索完成后，原来的起始点 出栈, 第一个点为上一步的搜索点
        empty_stepPoint = ~isempty(stepPoint);
    end
end
WaterHead(:, 1) = UpY - WaterHead(:, 1) .* cellsize + cellsize / 2; % row，Y坐标
WaterHead(:, 2) = LeftX + WaterHead(:, 2) .* cellsize - cellsize / 2; % Col，X坐标
length_WaterHead = length(WaterHead);

% 输出文件
fid = fopen(outputFilePath, 'w');
for j = 1:length_WaterHead
    fprintf(fid, '%f\t', WaterHead(j, 2)); % 先输入X坐标
    fprintf(fid, '%f\n', WaterHead(j, 1)); % 然后输入Y坐标
end
fclose(fid);
