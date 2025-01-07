function RiverProAnalysis(workspace)

%% 可以直接生成矢量文件 较原来程序 增添 界面选点,可以从不同图上选点 汇水面积识别功能，参数也增多
% 输出文件：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，max,min），折点高程(max,min)、最高折点的溯源距离、chi值
% 读取河流剖面数据，选择折点 RiverPro - long river
% profile：纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离 Line 18（高程剖面存储路径）；294,
% 298河段shapefile文件存储路径；298河段分段点存储路径

close all;
format long;

% 绘图时候的坐标显示约束
minLim_Area = 10^1;
maxLim_Area = 10^9;
% minLim_Len=0;maxLim_Len=5e4;
minLim_Slope = 10^(-3);
maxLim_Slope = 10^0;
MaxChi = 20;
MinEle = 900;
MaxEle = MinEle + 1100; %(floor(max(Ele)/100)+1)*100; 区域相对高差最大5000，而河流最低处约为500

% ******
str1 = strcat(workspace, 'profile/', "profile_");
river_i = input('河流编号： ');
str2 = num2str(river_i);
str3 = '.txt';
str_i = strcat(str1, str2, str3);
riverPro_old = textread(str_i);

% 计算logSA，新的riverPro 除了原先的，还包括平滑高程、局部ksn、坡度
SmoothLen = 1000;
ResampleWin = 20; %单位 m，平滑距离，高程重采样距离
riverPro = CalLogSA(riverPro_old, SmoothLen, ResampleWin);
riverPro_old = [];
[m, n] = size(riverPro);

% 存储新的河流剖面: 纵剖面点的平面坐标Y、X、溯源距离(3)、高程(4)、流向、汇水面积(6)、Chi距离、平滑高程、局部ksn、坡度(10)
str_iNew = strcat(str1, str2, 'New', str3);
fid1 = fopen(str_iNew, 'w');
for i = 1:m
    for j = 1:n - 1
        fprintf(fid1, '%f\t', riverPro(i, j));
    end
    fprintf(fid1, '%f', riverPro(i, n));
    fprintf(fid1, '\n');
end
fclose(fid1);

% ***** plot figures  ******
h1 = figure;
set(h1, 'units', 'normalized', 'position', [0.05, 0.05, 0.6, 0.85]); % Fig1; x-log(A);x-log(S);log(S-A)
% x-log(Area)
subplot(3, 1, 1);
set(gca, 'position', [0.1, 0.71, 0.85, 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 6));
% axis([minLim_Len,maxLim_Len,minLim_Area,maxLim_Area]);
ylim([minLim_Area, maxLim_Area]);
xlabel('Length (m)');
ylabel('Area (m^2)');
hold on;
% x-log(Slope)
subplot(3, 1, 2);
set(gca, 'position', [0.1, 0.39, 0.85, 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 10), '+');
% axis([minLim_Len,maxLim_Len,minLim_Slope,maxLim_Slope]);
ylim([minLim_Slope, maxLim_Slope]);
xlabel('Length (m)');
ylabel('Slope');
hold on;
% log(Area-Slope)
h13 = subplot(3, 1, 3);
set(gca, 'position', [0.1, 0.07, 0.85, 0.25]);
loglog(riverPro(:, 6), riverPro(:, 10), '+');
axis([minLim_Area, maxLim_Area, minLim_Slope, maxLim_Slope]);
xlabel('Area (m^2)');
ylabel('Slope');
hold on;
% subplot(2,1,1); plot(riverPro(:,3),riverPro(:,4));hold on;
% subplot(2,1,2);loglog(riverPro(:,6),riverPro(:,10),'+');axis([10^3,10^11,10^(-4),10^0]);hold
% on; % 既定范围 Wobus et al., 2006; Whipple et al., 2007

%% 根据汇水面积识别河道部分，对于崩积河道 剔除
h2 = figure;
set(h2, 'units', 'normalized', 'position', [0.05, 0.05, 0.85, 0.5]); % Fig2; x-Elevation; chi-Elevation
% x-Elevation;
h21 = subplot(1, 2, 1);
set(gca, 'position', [0.07, 0.1, 0.4, 0.85]);
plot(riverPro(:, 3), riverPro(:, 4));
% axis([minLim_Len,maxLim_Len,MinEle,MaxEle]);
ylim([MinEle, MaxEle]);
xlabel('Length (m)');
ylabel('Elevation (m)');
hold on;

disp('Choose A_cr,according to area')
figure(h1);
subplot(h13);
[Ac, Angle_limit] = ginput(1);
%Ac=input('汇水面积临界值m^2：'); %汇水面积临界值，超过这个值域，才是河道
disp(strcat('Acr=', num2str(Ac/1e6), 'km_2'));
i = m;
while riverPro(i, 6) < Ac
    i = i - 1;
end

figure(2);
subplot(h21);
plot(riverPro(i, 3), riverPro(i, 4), '*');
hold on; %在UpLen-z上标记河道顶点
riverPro = riverPro(1:i, :); % 对于崩积河道 剔除

% Y、X坐标，溯源距离、高程、流向、汇水面积、Chi距离
Chi = riverPro(:, 7);
Ele = riverPro(:, 4);
m = length(Chi);

figure(h2);
h22 = subplot(1, 2, 2);
set(gca, 'position', [0.55, 0.1, 0.4, 0.85]);
plot(Chi, Ele, 'g', 'linewidth', 3);
axis([0, MaxChi, MinEle, MaxEle]);

%% Fig3，用于DW校正后的chi-z
h3 = figure; % plot(Chi,Ele,'r','linewidth',3);axis([0,MaxChi,MinEle,MaxEle]);hold on;

% %% Fig4，UpLen-area; area-Elevation figure;
% subplot(2,1,1);semilogy(riverPro(:,3),riverPro(:,6)); hold on;
% subplot(2,1,2); semilogx(riverPro(:,6),Ele,'+');
% axis([10^3,10^11,MinEle,MaxEle]);hold on;

%% 根据高程剔除
% 1.根据高程剔除小的高程erase minor elevation
minorEle = input('erase minor elevation -1 or else：minorEle=');
while minorEle > 0
    for iErase = 1:m
        if Ele(iErase) > minorEle
            riverPro = riverPro(iErase:m, :);
            break;
        end % 找到大于设置的小值的全部留下；结束循环
    end
    minorEle = input('erase minor elevation -1 or else：minorEle=');
end
% 剔除小的高程后新的河流纵剖面
Chi = riverPro(:, 7);
Ele = riverPro(:, 4);
m = length(Chi);
figure(h2);
subplot(h22);
plot(Chi, Ele, 'm', 'linewidth', 3);
axis([0, MaxChi, MinEle, MaxEle]); % Fig2

% 2.根据高程剔除大的高程erase larger elevation
largeEle = input('erase larger elevation -1 or elese：largeEle=');
while largeEle > 0
    for jErase = 1:m
        if Ele(jErase) > largeEle
            riverPro = riverPro(1:(jErase - 1), :);
            break;
        end
    end
    largeEle = input('erase larger elevation -1 or elese：largeEle=');
end

% 剔除设置的高程后新的河流纵剖面 % set(gca,'ytick',MinEle:100:MaxEle);%
% text(Chi(1:2:m),Ele(1:2:m),num2str(Ele(1:2:m)));
Chi = riverPro(:, 7);
Ele = riverPro(:, 4);
Area = riverPro(:, 6);
Slope = riverPro(:, 10);
m = length(Chi);
figure(h2);
subplot(h22);
plot(Chi, Ele, 'r', 'linewidth', 3);
axis([0, MaxChi, MinEle, MaxEle]);
xlabel('Chi (m)');
ylabel('Elevation (m)');
hold on;

% figure; % Fig 5，根据高程剔除后的Length_Ele; Area_slope figure(5); subplot(2,1,1);
% plot(riverPro(:,3),riverPro(:,4));hold on; figure(5);
% subplot(2,1,2);loglog(riverPro(:,6),riverPro(:,10),'+');axis([10^3,10^10,10^(-3),10^0]);hold
% on; % 既定范围 Wobus et al., 2006; Whipple et al., 2007

%% 程序主干
PntKs = [];
PntSum = input('Number of divide pnts：PntSum=');

%% 选择找寻分割点的方法：1—从logSA；2—从chi-z图
if PntSum > 0
    choose_fig = input('选择找寻分割点的方法：1—从logSA；2—从chi-z图。 ');

    if choose_fig == 1
        disp(strcat('Upstream choose ', num2str(PntSum), ' points (according to Area):'))
        figure(h1);
        subplot(h13);
        [PntArea, PntSlope] = ginput(PntSum);
        PntArea = sort(PntArea, 'descend');
        %         for i=1:PntSum; disp(strcat(num2str(PntArea(i)/1e6),'km2,
        %         '));end;
        PntNum = [];
        PntNum = [PntNum; 1]; % 分界点在高程剖面对应的序号，包含第一个点、最后一个起始点

        for di = 1:PntSum %第di个分界点
            for k = 1:m - 1
                if (PntArea(di) <= Area(k) && PntArea(di) >= Area(k+1))
                    PntNum = [PntNum; k];
                    disp(strcat(num2str(Area(k)/1e6), 'km2； ', num2str(Ele(k)), 'm'));
                    break;
                end
            end
        end

        PntNum = [PntNum; m]; % 高程间隔的编号,PntNum至少包含首尾2个
        PntSum = length(PntNum) - 2; % 实际的分割点数目，PntNum是包含首尾的，因此，自然比分割点多2个

    elseif choose_fig == 2
        disp(strcat('Upstream choose ', num2str(PntSum), 'points (according to chi):'))
        figure(h2);
        subplot(h22); % PntEle=[]; %存储分界点大致高程段
        [PntChi, PntEle] = ginput(PntSum);
        PntChi = sort(PntChi);
        %         for i=1:PntSum; disp(strcat(num2str(PntEle(i)),', '));end;
        % for i=1:PntSum;
        % PntEle_i=input(strcat('从小到大，输入第',num2str(i),'个分界点高程:'));
        % PntEle=[PntEle;PntEle_i]; end 寻找确切的分界点编号
        PntNum = [];
        PntNum = [PntNum; 1]; % 分界点在高程剖面对应的序号，包含第一个点、最后一个起始点
        for j = 1:PntSum %第j个分界点的高程
            for k = 1:m - 1
                if (PntChi(j) >= Chi(k) && PntChi(j) <= Chi(k+1))
                    PntNum = [PntNum; k];
                    disp(strcat(num2str(Area(k)/1e6), 'km2； ', num2str(Ele(k)), 'm'));
                    break;
                end
            end
        end
        PntNum = [PntNum; m]; % 高程间隔的编号
        PntSum = length(PntNum) - 2; % PntNum是包含首尾的，因此，自然比分割点多2个
    end
    disp(strcat('最终，从大到小，实际分界点', num2str(PntSum), '个'))
end

%%
if PntSum == 0
    [R2, z0, ksn, std_ksn] = Bi_Regress(Ele, Chi); % 线性回归 Y=a0+a1*X
    figure(h2);
    subplot(h22);
    plot(Chi, Chi.*ksn+z0, 'k');
    hold on;
    text(max(Chi)/2, max(Ele), ['R=', num2str(R2), 'z=Chi.*', num2str(ksn), '±', num2str(std_ksn), '+', num2str(z0)]); % 这个是在Fig2上面画的
    [DW_r, DW_Ele, DW_Chi] = DW_test(Ele, Chi, z0, ksn); %DW检验
    [R2, DWz0, DWksn, DWstd_ksn] = Bi_Regress(DW_Ele, DW_Chi);
    figure(h3);
    plot(DW_Chi, DW_Ele, 'bo');
    hold on;
    plot(DW_Chi, DW_Chi.*DWksn+DWz0, 'k');
    text(max(DW_Chi)/2, max(DW_Ele), ['ro=', num2str(DW_r), 'R=', num2str(R2), 'ksn=', num2str(DWksn), '±', num2str(DWstd_ksn)]);
    [R2, log_ks, theta, std_theta] = Bi_Regress(log(Slope), log(Area));
    figure(h1);
    subplot(h13);
    loglog(Area, (Area.^(theta)).*(exp(log_ks)), 'k');
    hold on;
    text(max(Area)/2, 0.1, ['R=', num2str(R2), 'θ=', num2str(theta), '±', num2str(std_theta)]);

    %% 输出分界点信息：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
    PntKs = [PntKs; [river_i, riverPro(m, 2), riverPro(m, 1), theta, std_theta, ksn, std_ksn, DWksn, DWstd_ksn, (riverPro(m, 6)) ./ (1e6), (riverPro(1, 6)) ./ (1e6), riverPro(m, 4), riverPro(1, 4), riverPro(m, 3), riverPro(m, 7)]];

    %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
    Sfea = struct();

    Sfea.Geometry = 'Line';
    Xcor = riverPro(:, 2);
    Xcor = Xcor';
    Sfea.X = [Xcor, NaN];
    Ycor = riverPro(:, 1);
    Ycor = Ycor';
    Sfea.Y = [Ycor, NaN];
    Sfea.ID = 0;
    Sfea.RiverNo = river_i;
    Sfea.theta = theta;
    Sfea.std_theta = std_theta;
    Sfea.ksn = ksn;
    Sfea.std_ksn = std_ksn;
    Sfea.DWksn = DWksn;
    Sfea.DWstd_ksn = DWstd_ksn;
    Sfea.Amin_km = (riverPro(m, 6)) ./ (1e6);
    Sfea.Amax_km = (riverPro(1, 6)) ./ (1e6);
    Sfea.Emax_m = riverPro(m, 4);
    Sfea.Emin_m = riverPro(1, 4);
    Sfea.Lenmax_m = riverPro(m, 3);
    Sfea.Chimax = riverPro(m, 7);

    shapewrite(Sfea, strcat(workspace, '\Sfeature', num2str(river_i), '.shp'));

else % >=1个分界点

    %% 从哪个图上找分界点
    %     choose_fig=input('选择找寻分割点的方法：1—从logSA；2—从chi-z图。 ');
    %
    %     if choose_fig==1
    %         disp(strcat('从大到小，输入',num2str(PntSum),'个分界点面积:'))
    %         figure(5);subplot(2,1,2);  [PntArea,PntSlope]=ginput(PntSum);
    %         PntArea=sort(PntArea,'descend');
    % %         for i=1:PntSum; disp(strcat(num2str(PntArea(i)/1e6),'km2,
    % '));end;
    %         PntNum=[]; PntNum=[PntNum;1]; % 分界点在高程剖面对应的序号，包含第一个点、最后一个起始点
    %
    %         for di=1:PntSum %第di个分界点
    %             for k=1:m-1
    %                 if (PntArea(di)<=Area(k) && PntArea(di)>=Area(k+1))
    %                     PntNum=[PntNum;k];
    %                     disp(strcat(num2str(Area(k)/1e6),'km2；
    %                     ',num2str(Ele(k)),'m')); break;
    %                 end
    %             end
    %         end
    %
    %         PntNum=[PntNum;m]; % 高程间隔的编号,PntNum至少包含首尾2个
    %         PntSum=length(PntNum)-2; % PntNum是包含首尾的，因此，自然比分割点多2个
    %
    %     elseif choose_fig==2
    %         disp(strcat('从小到大，输入',num2str(PntSum),'个分界点高程:')) figure(2); %
    %         PntEle=[]; %存储分界点大致高程段 [PntChi,PntEle]=ginput(PntSum);
    %         PntEle=sort(PntEle);
    % %         for i=1:PntSum; disp(strcat(num2str(PntEle(i)),', '));end;
    %     % for i=1:PntSum;
    %     PntEle_i=input(strcat('从小到大，输入第',num2str(i),'个分界点高程:'));
    %     PntEle=[PntEle;PntEle_i]; end % 寻找确切的分界点编号
    %         PntNum=[]; PntNum=[PntNum;1]; % 分界点在高程剖面对应的序号，包含第一个点、最后一个起始点 for
    %         j=1:PntSum %第j个分界点的高程
    %             for k=1:m-1;
    %                 if (PntEle(j)>=Ele(k) && PntEle(j)<=Ele(k+1));
    %                     PntNum=[PntNum;k];
    %                     disp(strcat(num2str(Area(k)/1e6),'km2；
    %                     ',num2str(Ele(k)),'m')); break;
    %                 end;
    %             end
    %         end PntNum=[PntNum;m]; % 高程间隔的编号 PntSum=length(PntNum)-2; %
    %         PntNum是包含首尾的，因此，自然比分割点多2个
    %     end

    %% 填充 PntKs    % 第一段
    Ele_1 = Ele(PntNum(1):PntNum(2));
    Chi_1 = Chi(PntNum(1):PntNum(2));
    Area_1 = Area(PntNum(1):PntNum(2));
    Slope_1 = Slope(PntNum(1):PntNum(2));
    figure(h2);
    subplot(h21);
    plot(riverPro(PntNum(2), 3), riverPro(PntNum(2), 4), 'kx');
    hold on; %在UpLen-z上标记河道裂点
    [R2, z0_1, ksn_1, std_ksn_1] = Bi_Regress(Ele_1, Chi_1); % 线性回归 Y=a0+a1*X
    figure(h2);
    subplot(h22);
    plot(Chi_1, Chi_1.*ksn_1+z0_1, 'k');
    hold on;
    text(max(Chi_1)/2, max(Ele_1), ['R1=', num2str(R2), 'z1=Chi.*', num2str(ksn_1), '±', num2str(std_ksn_1), '+', num2str(z0_1)]); % 这个是在Fig2上面画的
    [DW_r_1, DW_Ele_1, DW_Chi_1] = DW_test(Ele_1, Chi_1, z0_1, ksn_1); %DW检验
    [R2, DWz0_1, DWksn_1, DWstd_ksn_1] = Bi_Regress(DW_Ele_1, DW_Chi_1);
    figure(h3);
    plot(DW_Chi_1, DW_Ele_1, 'bo');
    hold on;
    plot(DW_Chi_1, DW_Chi_1.*DWksn_1+DWz0_1, 'k');
    text(max(DW_Chi_1)/2, max(DW_Ele_1), ['ro1=', num2str(DW_r_1), 'R=', num2str(R2), 'ksn=', num2str(DWksn_1), '±', num2str(DWstd_ksn_1)]);
    [R2, log_ks_1, theta_1, std_theta_1] = Bi_Regress(log(Slope_1), log(Area_1));
    figure(h1);
    subplot(h13);
    loglog(Area_1, (Area_1.^(theta_1)).*(exp(log_ks_1)), 'k');
    hold on;
    text(max(Area_1)/2, max(Slope_1)/2, ['R1=', num2str(R2), 'θ=', num2str(theta_1), '±', num2str(std_theta_1)]);

    %% 输出分界点信息：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
    PntKs = [PntKs; [river_i, riverPro(PntNum(2), 2), riverPro(PntNum(2), 1), theta_1, std_theta_1, ksn_1, std_ksn_1, DWksn_1, DWstd_ksn_1, (riverPro(PntNum(2), 6)) ./ (1e6), (riverPro(PntNum(1), 6)) ./ (1e6), ...
        riverPro(PntNum(2), 4), riverPro(PntNum(1), 4), riverPro(PntNum(2), 3), riverPro(PntNum(2), 7)]];

    %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
    Sfea = struct();

    Sfea(1).Geometry = 'Line';
    Xcor = riverPro(PntNum(1):PntNum(2), 2);
    Xcor = Xcor';
    Sfea(1).X = [Xcor, NaN];
    Ycor = riverPro(PntNum(1):PntNum(2), 1);
    Ycor = Ycor';
    Sfea(1).Y = [Ycor, NaN];
    Sfea(1).ID = 0;
    Sfea(1).RiverNo = river_i;
    Sfea(1).theta = theta_1;
    Sfea(1).std_theta = std_theta_1;
    Sfea(1).ksn = ksn_1;
    Sfea(1).std_ksn = std_ksn_1;
    Sfea(1).DWksn = DWksn_1;
    Sfea(1).DWstd_ksn = DWstd_ksn_1;
    Sfea(1).Amin_km = (riverPro(PntNum(2), 6)) ./ (1e6);
    Sfea(1).Amax_km = (riverPro(PntNum(1), 6)) ./ (1e6);
    Sfea(1).Emax_m = riverPro(PntNum(2), 4);
    Sfea(1).Emin_m = riverPro(PntNum(1), 4);
    Sfea(1).Lenmax_m = riverPro(PntNum(2), 3);
    Sfea(1).Chimax = riverPro(PntNum(2), 7);

    %% ************剩余的河段*************************************************
    for j = 2:(PntSum + 1)
        Ele_j = Ele((PntNum(j) + 1):PntNum(j+1));
        Chi_j = Chi((PntNum(j) + 1):PntNum(j+1));
        Area_j = Area((PntNum(j) + 1):PntNum(j+1));
        Slope_j = Slope((PntNum(j) + 1):PntNum(j+1));
        figure(h2);
        subplot(h21);
        plot(riverPro(PntNum(j+1), 3), riverPro(PntNum(j+1), 4), 'kx');
        hold on; %在UpLen-z上标记河道裂点
        [R2, z0_j, ksn_j, std_ksn_j] = Bi_Regress(Ele_j, Chi_j); % 线性回归 Y=a0+a1*X
        figure(h2);
        subplot(h22);
        plot(Chi_j, Chi_j.*ksn_j+z0_j, 'k');
        hold on;
        text(max(Chi_j)/2, max(Ele_j), ['R', num2str(j), '=', num2str(R2), 'z=Chi.*', num2str(ksn_j), '±', num2str(std_ksn_j), '+', num2str(z0_j)]); % 这个是在Fig2上面画的
        [DW_r_j, DW_Ele_j, DW_Chi_j] = DW_test(Ele_j, Chi_j, z0_j, ksn_j); %DW检验
        [R2, DWz0_j, DWksn_j, DWstd_ksn_j] = Bi_Regress(DW_Ele_j, DW_Chi_j);
        figure(h3);
        plot(DW_Chi_j, DW_Ele_j, 'bo');
        hold on;
        plot(DW_Chi_j, DW_Chi_j.*DWksn_j+DWz0_j, 'k');
        text(max(DW_Chi_j)/2, max(DW_Ele_j), ['ro', num2str(j), '=', num2str(DW_r_j), 'R=', num2str(R2), 'ksn=', num2str(DWksn_j), '±', num2str(DWstd_ksn_j)]);
        %         figure (3); plot(Chi_j,Chi_j.*DWksn_j+DWz0_j/(1-DW_r_j),'k');hold
        %         on;
        %         text(max(Chi_j)/2,max(Ele_j)/2,['zj=Chi.*',num2str(DWksn_j),'±',num2str(DWstd_ksn_j),'+',num2str(DWz0_j/(1-DW_r_j))]);
        %         % 这个是在Fig3上面画的
        [R2, log_ks_j, theta_j, std_theta_j] = Bi_Regress(log(Slope_j), log(Area_j));
        figure(h1);
        subplot(h13);
        loglog(Area_j, (Area_j.^(theta_j)).*(exp(log_ks_j)), 'k');
        hold on;
        text(max(Area_j)/2, max(Slope_j)/2, ['R', num2str(j), '=', num2str(R2), 'θ=', num2str(theta_j), '±', num2str(std_theta_j)]);

        %% 输出分界点信息：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
        PntKs = [PntKs; [river_i, riverPro(PntNum(j+1), 2), riverPro(PntNum(j+1), 1), theta_j, std_theta_j, ksn_j, std_ksn_j, DWksn_j, DWstd_ksn_j, (riverPro(PntNum(j+1), 6)) ./ (1e6), (riverPro(PntNum(j), 6)) ./ (1e6), ...
            riverPro(PntNum(j+1), 4), riverPro(PntNum(j), 4), riverPro(PntNum(j+1), 3), riverPro(PntNum(j+1), 7)]];

        %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
        Sfea(j).Geometry = 'Line';
        Xcor = riverPro((PntNum(j) + 1):PntNum(j+1), 2);
        Xcor = Xcor';
        Sfea(j).X = [Xcor, NaN];
        Ycor = riverPro((PntNum(j) + 1):PntNum(j+1), 1);
        Ycor = Ycor';
        Sfea(j).Y = [Ycor, NaN];
        Sfea(j).ID = 0;
        Sfea(j).RiverNo = river_i;
        Sfea(j).theta = theta_j;
        Sfea(j).std_theta = std_theta_j;
        Sfea(j).ksn = ksn_j;
        Sfea(j).std_ksn = std_ksn_j;
        Sfea(j).DWksn = DWksn_j;
        Sfea(j).DWstd_ksn = DWstd_ksn_j;
        Sfea(j).Amin_km = (riverPro(PntNum(j+1), 6)) ./ (1e6);
        Sfea(j).Amax_km = (riverPro(PntNum(j), 6)) ./ (1e6);
        Sfea(j).Emax_m = riverPro(PntNum(j+1), 4);
        Sfea(j).Emin_m = riverPro(PntNum(j), 4);
        Sfea(j).Lenmax_m = riverPro(PntNum(j+1), 3);
        Sfea(j).Chimax = riverPro(PntNum(j+1), 7);
    end

    shapewrite(Sfea, strcat(workspace, 'Sfeature', num2str(river_i), '.shp'));
end

%% PntKs 到 D:\ErosionQilian\MatlabWork\basin_lons1west\result_lons1\PntKs.txt
[Pnt_m, Pnt_n] = size(PntKs);
str_PntKs = strcat(workspace, 'Pnt_divide.txt');
fid = fopen(str_PntKs, 'a+');
for i = 1:Pnt_m
    for j = 1:Pnt_n - 1
        fprintf(fid, '%f\t', PntKs(i, j));
    end
    fprintf(fid, '%f', PntKs(i, Pnt_n));
    fprintf(fid, '\n');
end
fclose(fid);

%% 计算坡度
function riverPro_f = CalLogSA(riverPro_Old_f, SmoothLen_f, ResampleWin_f)
% riverPro_i：纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离 新的到的riverPro
% 除了原先的，还包括平滑高程、局部ksn、坡度
UpLen = riverPro_Old_f(:, 3);
Ele = riverPro_Old_f(:, 4);
Chi = riverPro_Old_f(:, 7);
m = length(UpLen);

%计算平滑后的高程
SmoothedEle = Ele;
Localksn = zeros(m, 1);
Slope = zeros(m, 1);

for i = 1:m
    if abs(UpLen(i)-UpLen(1)) <= SmoothLen_f / 2 % 当点i与出水口距离小于SmoothLen_f/2，则一直向后寻找直至达到阈值为止
        j = 1;
        while i + j <= m && abs(UpLen(i+j)-UpLen(1)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if (i + j) - 1 < 2
            j = 3 - i;
        end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        SmoothedEle(i) = interp1(UpLen(1:i+j), Ele(1:i+j), UpLen(i));
        b = regress(Ele(1:i+j), [ones(size(Chi(1:i+j))), Chi(1:i+j)]); %regress函数用法，y=a0+a1*x
        Localksn(i) = abs(b(2));
    elseif abs(UpLen(m)-UpLen(i)) <= SmoothLen_f / 2 % 当点i与水头距离小于SmoothLen_f/2，则一直向前寻找直至达到阈值为止
        j = 1;
        while i - j >= 1 && abs(UpLen(m)-UpLen(i-j)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if (i - j) > (m - 2)
            j = i - (m - 2);
        end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        SmoothedEle(i) = interp1(UpLen(i-j:m), Ele(i-j:m), UpLen(i));
        b = regress(Ele(i-j:m), [ones(size(Chi(i-j:m))), Chi(i-j:m)]);
        Localksn(i) = abs(b(2));
    elseif abs(UpLen(i)-UpLen(1)) > SmoothLen_f / 2 && abs(UpLen(m)-UpLen(i)) > SmoothLen_f / 2 %点与出水口、水头距离都超过SmoothLen_f/2
        j = 1;
        while i - j >= 1 && i + j <= m && abs(UpLen(i+j)-UpLen(i-j)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if j == 0
            j = 1;
        end %当搜索点的左右相邻点的距离＜SmoothLen_f，强制选择相邻的点，确保插值、回归中有三个点。
        SmoothedEle(i) = interp1(UpLen(i-j:i+j), Ele(i-j:i+j), UpLen(i));
        b = regress(Ele(i-j:i+j), [ones(size(Chi(i-j:i+j))), Chi(i-j:i+j)]);
        Localksn(i) = abs(b(2));
    end
end

%计算重采样后的Localksn和Slope
for i = 1:m
    if abs(Ele(i)-Ele(1)) <= ResampleWin_f / 2 % 当点i与出水口高程小于ResampleWin_f/2，则一直向后寻找直至达到阈值为止
        j = 1;
        while i + j <= m && abs(Ele(i+j)-Ele(1)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if (i + j) - 1 < 2
            j = 3 - i;
        end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        b = regress(Ele(1:i+j), [ones(size(UpLen(1:i+j))), UpLen(1:i+j)]); %regress函数用法，y=a0+a1*x
        Slope(i) = abs(b(2));
    elseif abs(Ele(m)-Ele(i)) <= ResampleWin_f / 2 % 当点i与水头高程小于ResampleWin_f/2，则一直向前寻找直至达到阈值为止
        j = 1;
        while i - j >= 1 && abs(Ele(m)-Ele(i-j)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if (i - j) > (m - 2)
            j = i - (m - 2);
        end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        b = regress(Ele(i-j:m), [ones(size(UpLen(i-j:m))), UpLen(i-j:m)]);
        Slope(i) = abs(b(2));
    elseif abs(Ele(i)-Ele(1)) > ResampleWin_f / 2 && abs(Ele(m)-Ele(i)) > ResampleWin_f / 2 %点与出水口、水头高程都超过ResampleWin_f
        j = 1;
        while i - j >= 1 && i + j <= m && abs(Ele(i+j)-Ele(i-j)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if j == 0
            j = 1;
        end %当搜索点的左右相邻点的高差＜ResampleWin_f，强制选择相邻的点，确保插值、回归中有三个点。
        b = regress(Ele(i-j:i+j), [ones(size(UpLen(i-j:i+j))), UpLen(i-j:i+j)]);
        Slope(i) = abs(b(2));
    end
end

riverPro_f = [riverPro_Old_f, SmoothedEle, Localksn, Slope];

%% 统计检验
function [DW_r_f, DW_Ele_f, DW_Chi_f] = DW_test(Ele_f, Chi_f, z0_f, ksn_f)
DW_e = Ele_f - (z0_f + ksn_f .* Chi_f); % 残差
m = length(DW_e);
DW_e1 = DW_e(1:m-1);
DW_e2 = DW_e(2:m);
DW_r_f = sum(DW_e1.*DW_e2) / sqrt(sum(DW_e1.^2)) / sqrt(sum(DW_e2.^2));
DW_Chi_f = Chi_f(2:m) - DW_r_f .* Chi_f(1:m-1);
DW_Ele_f = Ele_f(2:m) - DW_r_f .* Ele_f(1:m-1);

%%
function [R2, a0, a1, std_a1] = Bi_Regress(Y, X)
% 线性回归 Y=a0+a1*X
Length_X = length(X);
if Length_X >= 3
    R2 = corrcoef(Y, X);
    if length(R2) > 1
        R2 = R2(1, 2);
        R2 = R2^2;
    else
        R2 = R2^2;
    end
    mean_X = mean(X);
    mean_Y = mean(Y);
    Sxx = sum((X - mean_X).^2);
    Syy = sum((Y - mean_Y).^2);
    Sxy = sum((X - mean_X).*(Y - mean_Y));
    a1 = Sxy / Sxx;
    a0 = mean_Y - a1 * mean_X;
    sigma = sqrt((Syy - a1 * Sxy)/(Length_X - 2));
    std_a1 = sigma / sqrt(Sxx);
elseif Length_X == 2
    R2 = 1;
    a1 = (Y(1) - Y(2)) ./ (X(1) - X(2));
    std_a1 = 0;
    a0 = Y(1) - a1 * X(1);
else
    %回归数据不够，不可模拟
    R2 = 9999;
    a0 = 9999;
    a1 = 9999;
    std_a1 = 9999;
end

%%
