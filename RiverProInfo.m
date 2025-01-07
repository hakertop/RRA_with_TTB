function RiverProInfo(Elev, FlowDir, FlowAcc, cellsize, headX, headY, LeftX, DownY, eleWin, mVSn, outputWorkspace)
% 溯源提取河流纵剖面 RiverPro - long river profile：纵剖面点的平面坐标YX、溯源距离、高程、流向、汇水面积、Chi距离
%
% cellsize - cellsizetion;
% headX,headY, - origion of the river;
% LeftX,DownY - Down corner of the raster
% eleWin - elevation window
% Line75：河道高程剖面数据存储路径

A0 = 1; %无维度化汇水面积1m2
% rownumber and colnumber of the river origin
[m, n] = size(Elev);
lenheadX = length(headX); % RiverPro=cell(lenheadX,1);

% RasterToASCII工具中得到的txt的左下角坐标，但行列号开始的是左上
upY = DownY + m * cellsize;

% 创建一个profile文件夹
if ~exist(outputWorkspace, "dir")
    mkdir(outputWorkspace);
end

for j = 1:lenheadX
    % 起始点坐标
    midRow = floor((abs(upY-headY(j)))/cellsize) + 1;
    midCol = floor((abs(headX(j)-LeftX))/cellsize) + 1;

    RiverLen = 0;
    RiverPro1 = [];
    kk = 1;

    % 首先填充顺流矩阵 RiverPro1
    RiverPro1 = [RiverPro1; midRow, midCol, RiverLen, Elev(midRow, midCol), FlowDir(midRow, midCol), FlowAcc(midRow, midCol)];

    while (midRow > 1 && midRow < m && midCol > 1 && midCol < n) % the searching point is within the matrix

        if Elev(midRow, midCol) < 0
            break;
        end % elevation is positive

        switch FlowDir(midRow, midCol)
            case 1 % East
                midRow = midRow;
                midCol = midCol + 1;
                RiverLen = RiverLen + 1;
                if FlowDir(midRow, midCol) == 16
                    break;
                end
            case 2 % SE
                midRow = midRow + 1;
                midCol = midCol + 1;
                RiverLen = RiverLen + 1.414;
                if FlowDir(midRow, midCol) == 32
                    break;
                end
            case 4 % South
                midRow = midRow + 1;
                midCol = midCol;
                RiverLen = RiverLen + 1;
                if FlowDir(midRow, midCol) == 64
                    break;
                end
            case 8 % SW
                midRow = midRow + 1;
                midCol = midCol - 1;
                RiverLen = RiverLen + 1.414;
                if FlowDir(midRow, midCol) == 128
                    break;
                end
            case 16 % Weat
                midRow = midRow;
                midCol = midCol - 1;
                RiverLen = RiverLen + 1;
                if FlowDir(midRow, midCol) == 1
                    break;
                end
            case 32 % NW
                midRow = midRow - 1;
                midCol = midCol - 1;
                RiverLen = RiverLen + 1.414;
                if FlowDir(midRow, midCol) == 2
                    break;
                end
            case 64 % North
                midRow = midRow - 1;
                midCol = midCol;
                RiverLen = RiverLen + 1;
                if FlowDir(midRow, midCol) == 4
                    break;
                end
            case 128 % NE
                midRow = midRow - 1;
                midCol = midCol + 1;
                RiverLen = RiverLen + 1.414;
                if FlowDir(midRow, midCol) == 8
                    break;
                end
        end

        if Elev(midRow, midCol) < 0
            break;
        end % elevation is positive

        if (RiverPro1(kk, 4) - Elev(midRow, midCol) >= eleWin) % beyond the elevation window
            RiverPro1 = [RiverPro1; midRow, midCol, RiverLen, Elev(midRow, midCol), FlowDir(midRow, midCol), FlowAcc(midRow, midCol)];
            kk = kk + 1;
        end
    end

    % 溯流矩阵 upstream distance
    jRiverPro = zeros(kk, 7);
    minEle = RiverPro1(kk, 4);
    MaxRiverLen = RiverPro1(kk, 3);

    for i = 1:kk
        jRiverPro(i, 1) = RiverPro1(kk-i+1, 1); % Row Number
        %jRiverPro(i,1)=upY-cellsize.*jRiverPro(i,1); %实际的地理坐标Y
        %jRiverPro(i,1)=upY-cellsize.*RiverPro{j,1}(i,1)
        jRiverPro(i, 2) = RiverPro1(kk-i+1, 2); % Col Number
        %jRiverPro(i,2)=LeftX+cellsize.*jRiverPro(i,2); %实际的地理坐标X
        %jRiverPro(i,2)=LeftX+cellsize.*RiverPro{j,1}(i,2)

        jRiverPro(i, 3) = (MaxRiverLen - RiverPro1(kk-i+1, 3)); % river length
        jRiverPro(i, 4) = RiverPro1(kk-i+1, 4); % elevation
        jRiverPro(i, 5) = RiverPro1(kk-i+1, 5); % Flow Direction
        jRiverPro(i, 6) = (RiverPro1(kk-i+1, 6)); % drainage area
    end

    jRiverPro(:, 1) = upY - cellsize .* jRiverPro(:, 1);
    jRiverPro(:, 2) = LeftX + cellsize .* jRiverPro(:, 2);
    jRiverPro(:, 3) = jRiverPro(:, 3) .* cellsize; % river length
    jRiverPro(:, 4) = jRiverPro(:, 4); % -minEle; % elevation
    jRiverPro(:, 6) = jRiverPro(:, 6) .* (cellsize^2);

    %计算Chi距离 (Perron and Royden,2012)
    for i = 2:kk
        jRiverPro(i, 7) = jRiverPro(i-1, 7) + ((A0 ./ jRiverPro(i-1, 6)).^mVSn) .* (jRiverPro(i, 3) - jRiverPro(i-1, 3));
    end

    % 输出
    outputFilePath = strcat(outputWorkspace, "profile_", num2str(j), ".txt");
    fid = fopen(outputFilePath, 'w');
    for i = 1:kk
        for k = 1:6
            fprintf(fid, '%f ', jRiverPro(i, k));
        end
        fprintf(fid, '%f\n', jRiverPro(i, 7));
    end
    fclose(fid);
end
