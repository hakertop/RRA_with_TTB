function Auto_ksMap(DEMFilePath, Elev, FlowDir, ChiMap, cellsize, headX, headY, LeftX, DownY, LenWin, outputWorkspace)

%% 直接生成矢量与栅格的，按照高程间隔计算ksMap;起始点 与 ChiPlot的起始点一样
% - DEMFilePath             DEM路径，主要用于生成结果矢量的坐标系信息
% - Elev                    高程二维矩阵
% - FlowDir                 流向二维矩阵
% - ChiMap                  ChiMap二维矩阵
% - cellsize                空间分辨率
% - headX, headY            计算pntKs的起始点<河流起点>（不是 ChiMap的起始点<河流终点>）
% - LeftX, DownY            左下角的地理坐标
% - LenWin                  局部河段长度, 即以此分割河道, 每一段的各个节点被赋予相同的陡峭指数
% - outputWorkspace         文件输出工作空间

%%
[m, n] = size(Elev);
lenheadX = length(headX); % RiverPro=cell(lenheadX,1);
ksMap = zeros(m, n); % ks分布数据
ksValue = [];
% RasterToASCII工具中得到的txt的左下角坐标，但行列号开始的是左上
upY = DownY + m * cellsize;
max_ksn = 1000; % 当超过最大值，则全部赋这个值
% MulNum=zeros(lenheadX,3);

Sfea_ksnMap = struct();
Sfea_num = 1;

%%
for j = 1:lenheadX

    %% midMtrix存储搜索到的栅格，信息：行号，列号，高程，chi；"kk"记录midMtrix行数；midMtrix_len记录该河段长度；kk_sfea记录河段中，实际纳入矢量显示的行数
    midMtrix = [];
    kk = 0;
    midMtrix_len = 0;
    kk_sfea = 0;

    % 从midRow,midCol开始寻找
    midRow = floor((abs(upY-headY(j)))/cellsize) + 1;
    midCol = floor((abs(headX(j)-LeftX))/cellsize) + 1;

    % 起始的搜索点越界，直接退出本次循环
    if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n
        continue;
    end

    % midRow,midCol虽然不越界，但这些点有可能ChiMap为0 （点在河道上，但面积有限，未被计算ChiMap）
    while (midRow > 1 && midRow < m && midCol > 1 && midCol < n && ChiMap(midRow, midCol) == 0) % 找到ChiMap的起始点，即ChiMap~=0，开始记录
        switch FlowDir(midRow, midCol)
            case 1
                midRow = midRow;
                midCol = midCol + 1; % East
            case 2
                midRow = midRow + 1;
                midCol = midCol + 1; % SE
            case 4
                midRow = midRow + 1;
                midCol = midCol; % South
            case 8
                midRow = midRow + 1;
                midCol = midCol - 1; % SW
            case 16
                midRow = midRow;
                midCol = midCol - 1; % West
            case 32
                midRow = midRow - 1;
                midCol = midCol - 1; % NW
            case 64
                midRow = midRow - 1;
                midCol = midCol; % North
            case 128
                midRow = midRow - 1;
                midCol = midCol + 1; % NE
        end
    end

    % 从开始搜索ChiMap非0点，但一直到越界了，都没找到，直接退出本次循环
    if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n
        continue;
    end

    % 初始搜索点
    midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
    kk = kk + 1;
    kk_sfea = kk_sfea + 1;

    %%
    while (midRow > 1 && midRow < m && midCol > 1 && midCol < n && ChiMap(midRow, midCol) > 0) % 该点是否值得搜索：必须不越矩阵的界，且有ChiMap值
        % while中的判断，可以保证，即便在行列号追加之后，还能不越界；因为，河流出山口的ChiMap是0值，所以这里需要的判别条件是ChiMap>0

        if FlowDir(midRow, midCol) == 1
            midRow = midRow;
            midCol = midCol + 1; % East
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize; % 不管ksnmap是否被计算，都纳入计算局部ksn的河段
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end % 但是，只有当ksnmap未被计算，才纳入矢量显示
        elseif FlowDir(midRow, midCol) == 2
            midRow = midRow + 1;
            midCol = midCol + 1; % SE
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize * sqrt(2);
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 4
            midRow = midRow + 1;
            midCol = midCol; % South
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize;
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 8
            midRow = midRow + 1;
            midCol = midCol - 1; % SW
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize * sqrt(2);
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 16
            midRow = midRow;
            midCol = midCol - 1; % West
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize;
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 32
            midRow = midRow - 1;
            midCol = midCol - 1; % NW
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize * sqrt(2);
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 64
            midRow = midRow - 1;
            midCol = midCol; % North
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize;
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        elseif FlowDir(midRow, midCol) == 128
            midRow = midRow - 1;
            midCol = midCol + 1; % NE
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize * sqrt(2);
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end
        end

        %%
        if midMtrix_len >= LenWin % 搜索点未过界，就已经超过距离限制
            % 计算河段的局部ksn值；midMtrix存储搜索到的栅格，信息：行号，列号，高程，chi；
            [Ks, Ks_error] = regress(midMtrix(:, 3), [ones(kk, 1), midMtrix(:, 4)]);
            if abs(Ks(2)) > max_ksn
                ksValue_1 = max_ksn;
            else
                ksValue_1 = abs(Ks(2));
            end
            ksValue = [ksValue; [ksValue_1, abs(Ks_error(2, 2)-Ks(2))]];

            % 生成栅格文件
            for i_kk = 1:kk_sfea
                ksMap(midMtrix(i_kk, 1), midMtrix(i_kk, 2)) = ksValue_1;
            end

            % 生成矢量文件  midMtrix 行号，列号，高程，chi
            Sfea_ksnMap(Sfea_num).Geometry = 'Line';
            Sfea_ksnMap(Sfea_num).ID = Sfea_num;
            if kk > kk_sfea % 在距离限制内，栅格被计算的ksn值
                CorX = midMtrix(1:(kk_sfea + 1), 2);
                CorX = CorX .* cellsize - cellsize / 2 + LeftX;
                Sfea_ksnMap(Sfea_num).X = [CorX; NaN]; %如果出错，就去掉"+1"
                CorY = midMtrix(1:(kk_sfea + 1), 1);
                CorY = upY - CorY .* cellsize + cellsize / 2;
                Sfea_ksnMap(Sfea_num).Y = [CorY; NaN];
                Sfea_ksnMap(Sfea_num).ksnValue = ksValue_1;
                Sfea_num = Sfea_num + 1;
                midMtrix = [];
                kk = 0;
                midMtrix_len = 0;
                kk_sfea = 0;
                break; % 当在距离限制内，下段栅格已经被计算ksn，一般出现在河流汇合处，退出由这个初始点开始while语句的河道搜索
            else % 在距离限制内，栅格未被计算的ksn值
                CorX = midMtrix(1:(kk_sfea), 2);
                CorX = CorX .* cellsize - cellsize / 2 + LeftX;
                Sfea_ksnMap(Sfea_num).X = [CorX; NaN];
                CorY = midMtrix(1:(kk_sfea), 1);
                CorY = upY - CorY .* cellsize + cellsize / 2;
                Sfea_ksnMap(Sfea_num).Y = [CorY; NaN];
                Sfea_ksnMap(Sfea_num).ksnValue = ksValue_1;
                Sfea_num = Sfea_num + 1;
                midMtrix = [];
                kk = 0;
                midMtrix_len = 0;
                kk_sfea = 0;
            end
            %             midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)];
            %             kk=kk+1; % midMtrix_len=midMtrix_len+cellsize*sqrt(2);

        else % 搜索点没有超过距离限制，进行是否越界判别
            if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n %|| ChiMap(midRow,midCol)<=0
                ChiMap_text = ChiMap(midRow, midCol);
                if ChiMap_text == 0 % 搜索点还没有超过距离限制，就已经过界
                    % 生成栅格文件
                    len_ksValue = length(ksValue(:, 1));
                    ksValue_1 = ksValue(len_ksValue, 1);
                    for i_kk = 1:kk_sfea
                        ksMap(midMtrix(i_kk, 1), midMtrix(i_kk, 2)) = ksValue_1;
                    end

                    % 生成矢量文件  midMtrix 行号，列号，高程，chi
                    Sfea_ksnMap(Sfea_num).Geometry = 'Line';
                    Sfea_ksnMap(Sfea_num).ID = Sfea_num;
                    CorX = midMtrix(1:(kk_sfea), 2);
                    CorX = CorX .* cellsize - cellsize / 2 + LeftX;
                    Sfea_ksnMap(Sfea_num).X = [CorX; NaN];
                    CorY = midMtrix(1:(kk_sfea), 1);
                    CorY = upY - CorY .* cellsize + cellsize / 2;
                    Sfea_ksnMap(Sfea_num).Y = [CorY; NaN];
                    Sfea_ksnMap(Sfea_num).ksnValue = ksValue_1;
                    Sfea_num = Sfea_num + 1;
                    midMtrix = [];
                    kk = 0;
                    midMtrix_len = 0;
                    kk_sfea = 0;
                    break;
                end
            end
        end

        %         if (ksMap(midRow,midCol)~=0 || Elev(midRow,midCol)<0)
        %             break; % 该点是否值得搜索：ksMap（空值=0）已经被计算 或者
        %             越高程界；那么下一个不必再寻找，直接退出while循环，寻找下一个起始点
        %         end if ChiMap(midRow,midCol)>0 %
        %         判断搜索出的点，是否满足要求（ChiMap非空，也是针对起始点的问题而言的）
        %             midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)];kk=kk+1;
        %         end
    end
end

shapefileOutputPath = strcat(outputWorkspace, 'ksn.shp');
shapewrite(Sfea_ksnMap, shapefileOutputPath);

% Zhaolong Dai 添加：补充空间坐标系信息
info = georasterinfo(DEMFilePath);
p = info.CoordinateReferenceSystem;
wkt = wktstring(p, 'Version', "wkt1");
fid = fopen(strcat(outputWorkspace, "ksn.prj"), "w");
fprintf(fid, "%s", wkt);
fclose(fid);

asciiGridOutputPath = strcat(outputWorkspace, 'ras_ksn.txt');
fid = fopen(asciiGridOutputPath, 'w');
for i = 1:m
    for j = 1:n - 1
        fprintf(fid, '%f ', ksMap(i, j));
    end
    fprintf(fid, '%f\n', ksMap(i, n));
end
fclose(fid);
