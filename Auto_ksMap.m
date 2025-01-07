function Auto_ksMap(DEMFilePath, Elev, FlowDir, ChiMap, cellsize, headX, headY, LeftX, DownY, LenWin, outputWorkspace)

%% ֱ������ʸ����դ��ģ����ո̼߳������ksMap;��ʼ�� �� ChiPlot����ʼ��һ��
% - DEMFilePath             DEM·������Ҫ�������ɽ��ʸ��������ϵ��Ϣ
% - Elev                    �̶߳�ά����
% - FlowDir                 �����ά����
% - ChiMap                  ChiMap��ά����
% - cellsize                �ռ�ֱ���
% - headX, headY            ����pntKs����ʼ��<�������>������ ChiMap����ʼ��<�����յ�>��
% - LeftX, DownY            ���½ǵĵ�������
% - LenWin                  �ֲ��Ӷγ���, ���Դ˷ָ�ӵ�, ÿһ�εĸ����ڵ㱻������ͬ�Ķ���ָ��
% - outputWorkspace         �ļ���������ռ�

%%
[m, n] = size(Elev);
lenheadX = length(headX); % RiverPro=cell(lenheadX,1);
ksMap = zeros(m, n); % ks�ֲ�����
ksValue = [];
% RasterToASCII�����еõ���txt�����½����꣬�����кſ�ʼ��������
upY = DownY + m * cellsize;
max_ksn = 1000; % ���������ֵ����ȫ�������ֵ
% MulNum=zeros(lenheadX,3);

Sfea_ksnMap = struct();
Sfea_num = 1;

%%
for j = 1:lenheadX

    %% midMtrix�洢��������դ����Ϣ���кţ��кţ��̣߳�chi��"kk"��¼midMtrix������midMtrix_len��¼�úӶγ��ȣ�kk_sfea��¼�Ӷ��У�ʵ������ʸ����ʾ������
    midMtrix = [];
    kk = 0;
    midMtrix_len = 0;
    kk_sfea = 0;

    % ��midRow,midCol��ʼѰ��
    midRow = floor((abs(upY-headY(j)))/cellsize) + 1;
    midCol = floor((abs(headX(j)-LeftX))/cellsize) + 1;

    % ��ʼ��������Խ�磬ֱ���˳�����ѭ��
    if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n
        continue;
    end

    % midRow,midCol��Ȼ��Խ�磬����Щ���п���ChiMapΪ0 �����ںӵ��ϣ���������ޣ�δ������ChiMap��
    while (midRow > 1 && midRow < m && midCol > 1 && midCol < n && ChiMap(midRow, midCol) == 0) % �ҵ�ChiMap����ʼ�㣬��ChiMap~=0����ʼ��¼
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

    % �ӿ�ʼ����ChiMap��0�㣬��һֱ��Խ���ˣ���û�ҵ���ֱ���˳�����ѭ��
    if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n
        continue;
    end

    % ��ʼ������
    midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
    kk = kk + 1;
    kk_sfea = kk_sfea + 1;

    %%
    while (midRow > 1 && midRow < m && midCol > 1 && midCol < n && ChiMap(midRow, midCol) > 0) % �õ��Ƿ�ֵ�����������벻Խ����Ľ磬����ChiMapֵ
        % while�е��жϣ����Ա�֤�����������к�׷��֮�󣬻��ܲ�Խ�磻��Ϊ��������ɽ�ڵ�ChiMap��0ֵ������������Ҫ���б�������ChiMap>0

        if FlowDir(midRow, midCol) == 1
            midRow = midRow;
            midCol = midCol + 1; % East
            midMtrix = [midMtrix; midRow, midCol, Elev(midRow, midCol), ChiMap(midRow, midCol)];
            kk = kk + 1;
            midMtrix_len = midMtrix_len + cellsize; % ����ksnmap�Ƿ񱻼��㣬���������ֲ�ksn�ĺӶ�
            if ksMap(midRow, midCol) == 0
                kk_sfea = kk_sfea + 1;
            end % ���ǣ�ֻ�е�ksnmapδ�����㣬������ʸ����ʾ
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
        if midMtrix_len >= LenWin % ������δ���磬���Ѿ�������������
            % ����Ӷεľֲ�ksnֵ��midMtrix�洢��������դ����Ϣ���кţ��кţ��̣߳�chi��
            [Ks, Ks_error] = regress(midMtrix(:, 3), [ones(kk, 1), midMtrix(:, 4)]);
            if abs(Ks(2)) > max_ksn
                ksValue_1 = max_ksn;
            else
                ksValue_1 = abs(Ks(2));
            end
            ksValue = [ksValue; [ksValue_1, abs(Ks_error(2, 2)-Ks(2))]];

            % ����դ���ļ�
            for i_kk = 1:kk_sfea
                ksMap(midMtrix(i_kk, 1), midMtrix(i_kk, 2)) = ksValue_1;
            end

            % ����ʸ���ļ�  midMtrix �кţ��кţ��̣߳�chi
            Sfea_ksnMap(Sfea_num).Geometry = 'Line';
            Sfea_ksnMap(Sfea_num).ID = Sfea_num;
            if kk > kk_sfea % �ھ��������ڣ�դ�񱻼����ksnֵ
                CorX = midMtrix(1:(kk_sfea + 1), 2);
                CorX = CorX .* cellsize - cellsize / 2 + LeftX;
                Sfea_ksnMap(Sfea_num).X = [CorX; NaN]; %���������ȥ��"+1"
                CorY = midMtrix(1:(kk_sfea + 1), 1);
                CorY = upY - CorY .* cellsize + cellsize / 2;
                Sfea_ksnMap(Sfea_num).Y = [CorY; NaN];
                Sfea_ksnMap(Sfea_num).ksnValue = ksValue_1;
                Sfea_num = Sfea_num + 1;
                midMtrix = [];
                kk = 0;
                midMtrix_len = 0;
                kk_sfea = 0;
                break; % ���ھ��������ڣ��¶�դ���Ѿ�������ksn��һ������ں�����ϴ����˳��������ʼ�㿪ʼwhile���ĺӵ�����
            else % �ھ��������ڣ�դ��δ�������ksnֵ
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

        else % ������û�г����������ƣ������Ƿ�Խ���б�
            if midRow <= 1 || midRow >= m || midCol <= 1 || midCol >= n %|| ChiMap(midRow,midCol)<=0
                ChiMap_text = ChiMap(midRow, midCol);
                if ChiMap_text == 0 % �����㻹û�г����������ƣ����Ѿ�����
                    % ����դ���ļ�
                    len_ksValue = length(ksValue(:, 1));
                    ksValue_1 = ksValue(len_ksValue, 1);
                    for i_kk = 1:kk_sfea
                        ksMap(midMtrix(i_kk, 1), midMtrix(i_kk, 2)) = ksValue_1;
                    end

                    % ����ʸ���ļ�  midMtrix �кţ��кţ��̣߳�chi
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
        %             break; % �õ��Ƿ�ֵ��������ksMap����ֵ=0���Ѿ������� ����
        %             Խ�߳̽磻��ô��һ��������Ѱ�ң�ֱ���˳�whileѭ����Ѱ����һ����ʼ��
        %         end if ChiMap(midRow,midCol)>0 %
        %         �ж��������ĵ㣬�Ƿ�����Ҫ��ChiMap�ǿգ�Ҳ�������ʼ���������Եģ�
        %             midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)];kk=kk+1;
        %         end
    end
end

shapefileOutputPath = strcat(outputWorkspace, 'ksn.shp');
shapewrite(Sfea_ksnMap, shapefileOutputPath);

% Zhaolong Dai ��ӣ�����ռ�����ϵ��Ϣ
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
