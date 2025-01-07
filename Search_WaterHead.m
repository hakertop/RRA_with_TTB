function WaterHead = Search_WaterHead(FlowAcc, FlowDir, LeftX, DownY, cellsize, ...
    ouletX, ouletY, Ac, outputFilePath)
% Search_WaterHead ��������ˮͷλ�ò����䱣�����ļ� �������:
%   FlowAcc - �����ۻ����� 
%   FlowDir - ������� 
%   LeftX, DownY - ��������½����� 
%   cellsize - ����Ĵ�С 
%   ouletX, ouletY - ���ڵ���������� 
%   Ac - ��ˮ����ֵ 
%   outputFilePath -����ı�·��(Ҫ�������������ļ������׺)

% ��ʼ����ر���
format long; % �������ݾ���
[m, n] = size(FlowAcc); % ��ȡ�����ۻ������������
UpY = DownY + cellsize * m; % ����������ϱ߽�Y����
Ac = floor(Ac/(cellsize^2)); % ����ˮ���ת���ɵ�Ԫ����Ŀ

len_ouletX = length(ouletX);
stepPoint = [];
WaterHead = [];

num_ToPoint = 0; % ����õ��դ����
num_slope = 0; % ����õ��դ���У�<Ac��դ����Ŀ
% ��������ȣ�˵���õ�Ϊˮͷ

visited = false(m, n);

for ilen = 1:len_ouletX
    % ��ʼ��������С��к�
    startRow = m - floor((ouletY(ilen) - DownY)/cellsize);
    startCol = floor((ouletX(ilen) - LeftX)/cellsize) + 1;

    % �������ֵ��Ч��
    if startRow >= 1 && startRow <= m && startCol >= 1 && startCol <= n
        stepPoint = [stepPoint; [startRow, startCol, FlowAcc(startRow, startCol), 0, 0]];
    else
        warning('��ʼ��(%d, %d)/��%d, %d��������Χ�������˵�', startRow, startCol, ouletX(ilen), ouletY(ilen));
        warning('UpY: %d, oultetY(ilen): %d', UpY, ouletY(ilen));
    end
    empty_stepPoint = 1;

    % �������С��кţ��õ�Ļ�ˮ������õ����������С��кţ�0,0��ʾ�޻��룩
    while empty_stepPoint
        startRow = stepPoint(1, 1);
        startCol = stepPoint(1, 2);

        if visited(startRow, startCol)
            stepPoint(1, :) = []; % % ÿһ��������ɺ�ԭ������ʼ�� ��ջ, ��һ����Ϊ��һ����������
            empty_stepPoint = ~isempty(stepPoint);
            continue;
        end
        visited(startRow, startCol) = true;

        if FlowAcc(startRow, startCol) >= Ac
            % ȷ���õ��Ǻ�����
            if startCol < n && FlowDir(startRow, startCol+1) == 2^4 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                % �õ��Ҳ�(>Ac)����ջ
                if FlowAcc(startRow, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow, startCol + 1, FlowAcc(startRow, startCol+1), startRow, startCol]];
                    % �õ��Ҳ�(<Ac)���õ�����Ϊˮͷ�Ŀ���
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && startCol < n && FlowDir(startRow+1, startCol+1) == 2^5 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow+1, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol + 1, FlowAcc(startRow+1, startCol+1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && FlowDir(startRow+1, startCol) == 2^6 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow+1, startCol) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol, FlowAcc(startRow+1, startCol), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow < m && startCol > 1 && FlowDir(startRow+1, startCol-1) == 2^7 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow+1,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow+1, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow + 1, startCol - 1, FlowAcc(startRow+1, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startCol > 1 && FlowDir(startRow, startCol-1) == 2^0 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow, startCol - 1, FlowAcc(startRow, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && startCol > 1 && FlowDir(startRow-1, startCol-1) == 2^1 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol-1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow-1, startCol-1) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol - 1, FlowAcc(startRow-1, startCol-1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && FlowDir(startRow-1, startCol) == 2^2 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow-1, startCol) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol, FlowAcc(startRow-1, startCol), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            if startRow > 1 && startCol < n && FlowDir(startRow-1, startCol+1) == 2^3 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow-1,startCol+1)>Ac
                num_ToPoint = num_ToPoint + 1; % ����õ�
                if FlowAcc(startRow-1, startCol+1) > Ac
                    stepPoint = [stepPoint; [startRow - 1, startCol + 1, FlowAcc(startRow-1, startCol+1), startRow, startCol]];
                else
                    num_slope = num_slope + 1;
                end
            end
            % �������Ŀ �� �����<Ac��Ŀ ��ȣ���õ�Ϊˮͷ
            if num_slope == num_ToPoint
                WaterHead = [WaterHead; [startRow, startCol]];
            end
        end
        num_slope = 0;
        num_ToPoint = 0;
        stepPoint(1, :) = []; % % ÿһ��������ɺ�ԭ������ʼ�� ��ջ, ��һ����Ϊ��һ����������
        empty_stepPoint = ~isempty(stepPoint);
    end
end
WaterHead(:, 1) = UpY - WaterHead(:, 1) .* cellsize + cellsize / 2; % row��Y����
WaterHead(:, 2) = LeftX + WaterHead(:, 2) .* cellsize - cellsize / 2; % Col��X����
length_WaterHead = length(WaterHead);

% ����ļ�
fid = fopen(outputFilePath, 'w');
for j = 1:length_WaterHead
    fprintf(fid, '%f\t', WaterHead(j, 2)); % ������X����
    fprintf(fid, '%f\n', WaterHead(j, 1)); % Ȼ������Y����
end
fclose(fid);
