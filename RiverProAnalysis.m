function RiverProAnalysis(workspace)

%% ����ֱ������ʸ���ļ� ��ԭ������ ���� ����ѡ��,���ԴӲ�ͬͼ��ѡ�� ��ˮ���ʶ���ܣ�����Ҳ����
% ����ļ���������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��max,min�����۵�߳�(max,min)������۵����Դ���롢chiֵ
% ��ȡ�����������ݣ�ѡ���۵� RiverPro - long river
% profile����������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���� Line 18���߳�����洢·������294,
% 298�Ӷ�shapefile�ļ��洢·����298�Ӷηֶε�洢·��

close all;
format long;

% ��ͼʱ���������ʾԼ��
minLim_Area = 10^1;
maxLim_Area = 10^9;
% minLim_Len=0;maxLim_Len=5e4;
minLim_Slope = 10^(-3);
maxLim_Slope = 10^0;
MaxChi = 20;
MinEle = 900;
MaxEle = MinEle + 1100; %(floor(max(Ele)/100)+1)*100; ������Ը߲����5000����������ʹ�ԼΪ500

% ******
str1 = strcat(workspace, 'profile/', "profile_");
river_i = input('������ţ� ');
str2 = num2str(river_i);
str3 = '.txt';
str_i = strcat(str1, str2, str3);
riverPro_old = textread(str_i);

% ����logSA���µ�riverPro ����ԭ�ȵģ�������ƽ���̡߳��ֲ�ksn���¶�
SmoothLen = 1000;
ResampleWin = 20; %��λ m��ƽ�����룬�߳��ز�������
riverPro = CalLogSA(riverPro_old, SmoothLen, ResampleWin);
riverPro_old = [];
[m, n] = size(riverPro);

% �洢�µĺ�������: ��������ƽ������Y��X����Դ����(3)���߳�(4)�����򡢻�ˮ���(6)��Chi���롢ƽ���̡߳��ֲ�ksn���¶�(10)
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
% on; % �ȶ���Χ Wobus et al., 2006; Whipple et al., 2007

%% ���ݻ�ˮ���ʶ��ӵ����֣����ڱ����ӵ� �޳�
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
%Ac=input('��ˮ����ٽ�ֵm^2��'); %��ˮ����ٽ�ֵ���������ֵ�򣬲��Ǻӵ�
disp(strcat('Acr=', num2str(Ac/1e6), 'km_2'));
i = m;
while riverPro(i, 6) < Ac
    i = i - 1;
end

figure(2);
subplot(h21);
plot(riverPro(i, 3), riverPro(i, 4), '*');
hold on; %��UpLen-z�ϱ�Ǻӵ�����
riverPro = riverPro(1:i, :); % ���ڱ����ӵ� �޳�

% Y��X���꣬��Դ���롢�̡߳����򡢻�ˮ�����Chi����
Chi = riverPro(:, 7);
Ele = riverPro(:, 4);
m = length(Chi);

figure(h2);
h22 = subplot(1, 2, 2);
set(gca, 'position', [0.55, 0.1, 0.4, 0.85]);
plot(Chi, Ele, 'g', 'linewidth', 3);
axis([0, MaxChi, MinEle, MaxEle]);

%% Fig3������DWУ�����chi-z
h3 = figure; % plot(Chi,Ele,'r','linewidth',3);axis([0,MaxChi,MinEle,MaxEle]);hold on;

% %% Fig4��UpLen-area; area-Elevation figure;
% subplot(2,1,1);semilogy(riverPro(:,3),riverPro(:,6)); hold on;
% subplot(2,1,2); semilogx(riverPro(:,6),Ele,'+');
% axis([10^3,10^11,MinEle,MaxEle]);hold on;

%% ���ݸ߳��޳�
% 1.���ݸ߳��޳�С�ĸ߳�erase minor elevation
minorEle = input('erase minor elevation -1 or else��minorEle=');
while minorEle > 0
    for iErase = 1:m
        if Ele(iErase) > minorEle
            riverPro = riverPro(iErase:m, :);
            break;
        end % �ҵ��������õ�Сֵ��ȫ�����£�����ѭ��
    end
    minorEle = input('erase minor elevation -1 or else��minorEle=');
end
% �޳�С�ĸ̺߳��µĺ���������
Chi = riverPro(:, 7);
Ele = riverPro(:, 4);
m = length(Chi);
figure(h2);
subplot(h22);
plot(Chi, Ele, 'm', 'linewidth', 3);
axis([0, MaxChi, MinEle, MaxEle]); % Fig2

% 2.���ݸ߳��޳���ĸ߳�erase larger elevation
largeEle = input('erase larger elevation -1 or elese��largeEle=');
while largeEle > 0
    for jErase = 1:m
        if Ele(jErase) > largeEle
            riverPro = riverPro(1:(jErase - 1), :);
            break;
        end
    end
    largeEle = input('erase larger elevation -1 or elese��largeEle=');
end

% �޳����õĸ̺߳��µĺ��������� % set(gca,'ytick',MinEle:100:MaxEle);%
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

% figure; % Fig 5�����ݸ߳��޳����Length_Ele; Area_slope figure(5); subplot(2,1,1);
% plot(riverPro(:,3),riverPro(:,4));hold on; figure(5);
% subplot(2,1,2);loglog(riverPro(:,6),riverPro(:,10),'+');axis([10^3,10^10,10^(-3),10^0]);hold
% on; % �ȶ���Χ Wobus et al., 2006; Whipple et al., 2007

%% ��������
PntKs = [];
PntSum = input('Number of divide pnts��PntSum=');

%% ѡ����Ѱ�ָ��ķ�����1����logSA��2����chi-zͼ
if PntSum > 0
    choose_fig = input('ѡ����Ѱ�ָ��ķ�����1����logSA��2����chi-zͼ�� ');

    if choose_fig == 1
        disp(strcat('Upstream choose ', num2str(PntSum), ' points (according to Area):'))
        figure(h1);
        subplot(h13);
        [PntArea, PntSlope] = ginput(PntSum);
        PntArea = sort(PntArea, 'descend');
        %         for i=1:PntSum; disp(strcat(num2str(PntArea(i)/1e6),'km2,
        %         '));end;
        PntNum = [];
        PntNum = [PntNum; 1]; % �ֽ���ڸ߳������Ӧ����ţ�������һ���㡢���һ����ʼ��

        for di = 1:PntSum %��di���ֽ��
            for k = 1:m - 1
                if (PntArea(di) <= Area(k) && PntArea(di) >= Area(k+1))
                    PntNum = [PntNum; k];
                    disp(strcat(num2str(Area(k)/1e6), 'km2�� ', num2str(Ele(k)), 'm'));
                    break;
                end
            end
        end

        PntNum = [PntNum; m]; % �̼߳���ı��,PntNum���ٰ�����β2��
        PntSum = length(PntNum) - 2; % ʵ�ʵķָ����Ŀ��PntNum�ǰ�����β�ģ���ˣ���Ȼ�ȷָ���2��

    elseif choose_fig == 2
        disp(strcat('Upstream choose ', num2str(PntSum), 'points (according to chi):'))
        figure(h2);
        subplot(h22); % PntEle=[]; %�洢�ֽ����¸̶߳�
        [PntChi, PntEle] = ginput(PntSum);
        PntChi = sort(PntChi);
        %         for i=1:PntSum; disp(strcat(num2str(PntEle(i)),', '));end;
        % for i=1:PntSum;
        % PntEle_i=input(strcat('��С���������',num2str(i),'���ֽ��߳�:'));
        % PntEle=[PntEle;PntEle_i]; end Ѱ��ȷ�еķֽ����
        PntNum = [];
        PntNum = [PntNum; 1]; % �ֽ���ڸ߳������Ӧ����ţ�������һ���㡢���һ����ʼ��
        for j = 1:PntSum %��j���ֽ��ĸ߳�
            for k = 1:m - 1
                if (PntChi(j) >= Chi(k) && PntChi(j) <= Chi(k+1))
                    PntNum = [PntNum; k];
                    disp(strcat(num2str(Area(k)/1e6), 'km2�� ', num2str(Ele(k)), 'm'));
                    break;
                end
            end
        end
        PntNum = [PntNum; m]; % �̼߳���ı��
        PntSum = length(PntNum) - 2; % PntNum�ǰ�����β�ģ���ˣ���Ȼ�ȷָ���2��
    end
    disp(strcat('���գ��Ӵ�С��ʵ�ʷֽ��', num2str(PntSum), '��'))
end

%%
if PntSum == 0
    [R2, z0, ksn, std_ksn] = Bi_Regress(Ele, Chi); % ���Իع� Y=a0+a1*X
    figure(h2);
    subplot(h22);
    plot(Chi, Chi.*ksn+z0, 'k');
    hold on;
    text(max(Chi)/2, max(Ele), ['R=', num2str(R2), 'z=Chi.*', num2str(ksn), '��', num2str(std_ksn), '+', num2str(z0)]); % �������Fig2���滭��
    [DW_r, DW_Ele, DW_Chi] = DW_test(Ele, Chi, z0, ksn); %DW����
    [R2, DWz0, DWksn, DWstd_ksn] = Bi_Regress(DW_Ele, DW_Chi);
    figure(h3);
    plot(DW_Chi, DW_Ele, 'bo');
    hold on;
    plot(DW_Chi, DW_Chi.*DWksn+DWz0, 'k');
    text(max(DW_Chi)/2, max(DW_Ele), ['ro=', num2str(DW_r), 'R=', num2str(R2), 'ksn=', num2str(DWksn), '��', num2str(DWstd_ksn)]);
    [R2, log_ks, theta, std_theta] = Bi_Regress(log(Slope), log(Area));
    figure(h1);
    subplot(h13);
    loglog(Area, (Area.^(theta)).*(exp(log_ks)), 'k');
    hold on;
    text(max(Area)/2, 0.1, ['R=', num2str(R2), '��=', num2str(theta), '��', num2str(std_theta)]);

    %% ����ֽ����Ϣ��������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
    PntKs = [PntKs; [river_i, riverPro(m, 2), riverPro(m, 1), theta, std_theta, ksn, std_ksn, DWksn, DWstd_ksn, (riverPro(m, 6)) ./ (1e6), (riverPro(1, 6)) ./ (1e6), riverPro(m, 4), riverPro(1, 4), riverPro(m, 3), riverPro(m, 7)]];

    %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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

else % >=1���ֽ��

    %% ���ĸ�ͼ���ҷֽ��
    %     choose_fig=input('ѡ����Ѱ�ָ��ķ�����1����logSA��2����chi-zͼ�� ');
    %
    %     if choose_fig==1
    %         disp(strcat('�Ӵ�С������',num2str(PntSum),'���ֽ�����:'))
    %         figure(5);subplot(2,1,2);  [PntArea,PntSlope]=ginput(PntSum);
    %         PntArea=sort(PntArea,'descend');
    % %         for i=1:PntSum; disp(strcat(num2str(PntArea(i)/1e6),'km2,
    % '));end;
    %         PntNum=[]; PntNum=[PntNum;1]; % �ֽ���ڸ߳������Ӧ����ţ�������һ���㡢���һ����ʼ��
    %
    %         for di=1:PntSum %��di���ֽ��
    %             for k=1:m-1
    %                 if (PntArea(di)<=Area(k) && PntArea(di)>=Area(k+1))
    %                     PntNum=[PntNum;k];
    %                     disp(strcat(num2str(Area(k)/1e6),'km2��
    %                     ',num2str(Ele(k)),'m')); break;
    %                 end
    %             end
    %         end
    %
    %         PntNum=[PntNum;m]; % �̼߳���ı��,PntNum���ٰ�����β2��
    %         PntSum=length(PntNum)-2; % PntNum�ǰ�����β�ģ���ˣ���Ȼ�ȷָ���2��
    %
    %     elseif choose_fig==2
    %         disp(strcat('��С��������',num2str(PntSum),'���ֽ��߳�:')) figure(2); %
    %         PntEle=[]; %�洢�ֽ����¸̶߳� [PntChi,PntEle]=ginput(PntSum);
    %         PntEle=sort(PntEle);
    % %         for i=1:PntSum; disp(strcat(num2str(PntEle(i)),', '));end;
    %     % for i=1:PntSum;
    %     PntEle_i=input(strcat('��С���������',num2str(i),'���ֽ��߳�:'));
    %     PntEle=[PntEle;PntEle_i]; end % Ѱ��ȷ�еķֽ����
    %         PntNum=[]; PntNum=[PntNum;1]; % �ֽ���ڸ߳������Ӧ����ţ�������һ���㡢���һ����ʼ�� for
    %         j=1:PntSum %��j���ֽ��ĸ߳�
    %             for k=1:m-1;
    %                 if (PntEle(j)>=Ele(k) && PntEle(j)<=Ele(k+1));
    %                     PntNum=[PntNum;k];
    %                     disp(strcat(num2str(Area(k)/1e6),'km2��
    %                     ',num2str(Ele(k)),'m')); break;
    %                 end;
    %             end
    %         end PntNum=[PntNum;m]; % �̼߳���ı�� PntSum=length(PntNum)-2; %
    %         PntNum�ǰ�����β�ģ���ˣ���Ȼ�ȷָ���2��
    %     end

    %% ��� PntKs    % ��һ��
    Ele_1 = Ele(PntNum(1):PntNum(2));
    Chi_1 = Chi(PntNum(1):PntNum(2));
    Area_1 = Area(PntNum(1):PntNum(2));
    Slope_1 = Slope(PntNum(1):PntNum(2));
    figure(h2);
    subplot(h21);
    plot(riverPro(PntNum(2), 3), riverPro(PntNum(2), 4), 'kx');
    hold on; %��UpLen-z�ϱ�Ǻӵ��ѵ�
    [R2, z0_1, ksn_1, std_ksn_1] = Bi_Regress(Ele_1, Chi_1); % ���Իع� Y=a0+a1*X
    figure(h2);
    subplot(h22);
    plot(Chi_1, Chi_1.*ksn_1+z0_1, 'k');
    hold on;
    text(max(Chi_1)/2, max(Ele_1), ['R1=', num2str(R2), 'z1=Chi.*', num2str(ksn_1), '��', num2str(std_ksn_1), '+', num2str(z0_1)]); % �������Fig2���滭��
    [DW_r_1, DW_Ele_1, DW_Chi_1] = DW_test(Ele_1, Chi_1, z0_1, ksn_1); %DW����
    [R2, DWz0_1, DWksn_1, DWstd_ksn_1] = Bi_Regress(DW_Ele_1, DW_Chi_1);
    figure(h3);
    plot(DW_Chi_1, DW_Ele_1, 'bo');
    hold on;
    plot(DW_Chi_1, DW_Chi_1.*DWksn_1+DWz0_1, 'k');
    text(max(DW_Chi_1)/2, max(DW_Ele_1), ['ro1=', num2str(DW_r_1), 'R=', num2str(R2), 'ksn=', num2str(DWksn_1), '��', num2str(DWstd_ksn_1)]);
    [R2, log_ks_1, theta_1, std_theta_1] = Bi_Regress(log(Slope_1), log(Area_1));
    figure(h1);
    subplot(h13);
    loglog(Area_1, (Area_1.^(theta_1)).*(exp(log_ks_1)), 'k');
    hold on;
    text(max(Area_1)/2, max(Slope_1)/2, ['R1=', num2str(R2), '��=', num2str(theta_1), '��', num2str(std_theta_1)]);

    %% ����ֽ����Ϣ��������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
    PntKs = [PntKs; [river_i, riverPro(PntNum(2), 2), riverPro(PntNum(2), 1), theta_1, std_theta_1, ksn_1, std_ksn_1, DWksn_1, DWstd_ksn_1, (riverPro(PntNum(2), 6)) ./ (1e6), (riverPro(PntNum(1), 6)) ./ (1e6), ...
        riverPro(PntNum(2), 4), riverPro(PntNum(1), 4), riverPro(PntNum(2), 3), riverPro(PntNum(2), 7)]];

    %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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

    %% ************ʣ��ĺӶ�*************************************************
    for j = 2:(PntSum + 1)
        Ele_j = Ele((PntNum(j) + 1):PntNum(j+1));
        Chi_j = Chi((PntNum(j) + 1):PntNum(j+1));
        Area_j = Area((PntNum(j) + 1):PntNum(j+1));
        Slope_j = Slope((PntNum(j) + 1):PntNum(j+1));
        figure(h2);
        subplot(h21);
        plot(riverPro(PntNum(j+1), 3), riverPro(PntNum(j+1), 4), 'kx');
        hold on; %��UpLen-z�ϱ�Ǻӵ��ѵ�
        [R2, z0_j, ksn_j, std_ksn_j] = Bi_Regress(Ele_j, Chi_j); % ���Իع� Y=a0+a1*X
        figure(h2);
        subplot(h22);
        plot(Chi_j, Chi_j.*ksn_j+z0_j, 'k');
        hold on;
        text(max(Chi_j)/2, max(Ele_j), ['R', num2str(j), '=', num2str(R2), 'z=Chi.*', num2str(ksn_j), '��', num2str(std_ksn_j), '+', num2str(z0_j)]); % �������Fig2���滭��
        [DW_r_j, DW_Ele_j, DW_Chi_j] = DW_test(Ele_j, Chi_j, z0_j, ksn_j); %DW����
        [R2, DWz0_j, DWksn_j, DWstd_ksn_j] = Bi_Regress(DW_Ele_j, DW_Chi_j);
        figure(h3);
        plot(DW_Chi_j, DW_Ele_j, 'bo');
        hold on;
        plot(DW_Chi_j, DW_Chi_j.*DWksn_j+DWz0_j, 'k');
        text(max(DW_Chi_j)/2, max(DW_Ele_j), ['ro', num2str(j), '=', num2str(DW_r_j), 'R=', num2str(R2), 'ksn=', num2str(DWksn_j), '��', num2str(DWstd_ksn_j)]);
        %         figure (3); plot(Chi_j,Chi_j.*DWksn_j+DWz0_j/(1-DW_r_j),'k');hold
        %         on;
        %         text(max(Chi_j)/2,max(Ele_j)/2,['zj=Chi.*',num2str(DWksn_j),'��',num2str(DWstd_ksn_j),'+',num2str(DWz0_j/(1-DW_r_j))]);
        %         % �������Fig3���滭��
        [R2, log_ks_j, theta_j, std_theta_j] = Bi_Regress(log(Slope_j), log(Area_j));
        figure(h1);
        subplot(h13);
        loglog(Area_j, (Area_j.^(theta_j)).*(exp(log_ks_j)), 'k');
        hold on;
        text(max(Area_j)/2, max(Slope_j)/2, ['R', num2str(j), '=', num2str(R2), '��=', num2str(theta_j), '��', num2str(std_theta_j)]);

        %% ����ֽ����Ϣ��������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
        PntKs = [PntKs; [river_i, riverPro(PntNum(j+1), 2), riverPro(PntNum(j+1), 1), theta_j, std_theta_j, ksn_j, std_ksn_j, DWksn_j, DWstd_ksn_j, (riverPro(PntNum(j+1), 6)) ./ (1e6), (riverPro(PntNum(j), 6)) ./ (1e6), ...
            riverPro(PntNum(j+1), 4), riverPro(PntNum(j), 4), riverPro(PntNum(j+1), 3), riverPro(PntNum(j+1), 7)]];

        %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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

%% PntKs �� D:\ErosionQilian\MatlabWork\basin_lons1west\result_lons1\PntKs.txt
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

%% �����¶�
function riverPro_f = CalLogSA(riverPro_Old_f, SmoothLen_f, ResampleWin_f)
% riverPro_i����������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���� �µĵ���riverPro
% ����ԭ�ȵģ�������ƽ���̡߳��ֲ�ksn���¶�
UpLen = riverPro_Old_f(:, 3);
Ele = riverPro_Old_f(:, 4);
Chi = riverPro_Old_f(:, 7);
m = length(UpLen);

%����ƽ����ĸ߳�
SmoothedEle = Ele;
Localksn = zeros(m, 1);
Slope = zeros(m, 1);

for i = 1:m
    if abs(UpLen(i)-UpLen(1)) <= SmoothLen_f / 2 % ����i���ˮ�ھ���С��SmoothLen_f/2����һֱ���Ѱ��ֱ���ﵽ��ֵΪֹ
        j = 1;
        while i + j <= m && abs(UpLen(i+j)-UpLen(1)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if (i + j) - 1 < 2
            j = 3 - i;
        end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i) = interp1(UpLen(1:i+j), Ele(1:i+j), UpLen(i));
        b = regress(Ele(1:i+j), [ones(size(Chi(1:i+j))), Chi(1:i+j)]); %regress�����÷���y=a0+a1*x
        Localksn(i) = abs(b(2));
    elseif abs(UpLen(m)-UpLen(i)) <= SmoothLen_f / 2 % ����i��ˮͷ����С��SmoothLen_f/2����һֱ��ǰѰ��ֱ���ﵽ��ֵΪֹ
        j = 1;
        while i - j >= 1 && abs(UpLen(m)-UpLen(i-j)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if (i - j) > (m - 2)
            j = i - (m - 2);
        end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i) = interp1(UpLen(i-j:m), Ele(i-j:m), UpLen(i));
        b = regress(Ele(i-j:m), [ones(size(Chi(i-j:m))), Chi(i-j:m)]);
        Localksn(i) = abs(b(2));
    elseif abs(UpLen(i)-UpLen(1)) > SmoothLen_f / 2 && abs(UpLen(m)-UpLen(i)) > SmoothLen_f / 2 %�����ˮ�ڡ�ˮͷ���붼����SmoothLen_f/2
        j = 1;
        while i - j >= 1 && i + j <= m && abs(UpLen(i+j)-UpLen(i-j)) <= SmoothLen_f
            j = j + 1;
        end
        j = j - 1;
        if j == 0
            j = 1;
        end %����������������ڵ�ľ��룼SmoothLen_f��ǿ��ѡ�����ڵĵ㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i) = interp1(UpLen(i-j:i+j), Ele(i-j:i+j), UpLen(i));
        b = regress(Ele(i-j:i+j), [ones(size(Chi(i-j:i+j))), Chi(i-j:i+j)]);
        Localksn(i) = abs(b(2));
    end
end

%�����ز������Localksn��Slope
for i = 1:m
    if abs(Ele(i)-Ele(1)) <= ResampleWin_f / 2 % ����i���ˮ�ڸ߳�С��ResampleWin_f/2����һֱ���Ѱ��ֱ���ﵽ��ֵΪֹ
        j = 1;
        while i + j <= m && abs(Ele(i+j)-Ele(1)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if (i + j) - 1 < 2
            j = 3 - i;
        end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        b = regress(Ele(1:i+j), [ones(size(UpLen(1:i+j))), UpLen(1:i+j)]); %regress�����÷���y=a0+a1*x
        Slope(i) = abs(b(2));
    elseif abs(Ele(m)-Ele(i)) <= ResampleWin_f / 2 % ����i��ˮͷ�߳�С��ResampleWin_f/2����һֱ��ǰѰ��ֱ���ﵽ��ֵΪֹ
        j = 1;
        while i - j >= 1 && abs(Ele(m)-Ele(i-j)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if (i - j) > (m - 2)
            j = i - (m - 2);
        end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        b = regress(Ele(i-j:m), [ones(size(UpLen(i-j:m))), UpLen(i-j:m)]);
        Slope(i) = abs(b(2));
    elseif abs(Ele(i)-Ele(1)) > ResampleWin_f / 2 && abs(Ele(m)-Ele(i)) > ResampleWin_f / 2 %�����ˮ�ڡ�ˮͷ�̶߳�����ResampleWin_f
        j = 1;
        while i - j >= 1 && i + j <= m && abs(Ele(i+j)-Ele(i-j)) <= ResampleWin_f
            j = j + 1;
        end
        j = j - 1;
        if j == 0
            j = 1;
        end %����������������ڵ�ĸ߲ResampleWin_f��ǿ��ѡ�����ڵĵ㣬ȷ����ֵ���ع����������㡣
        b = regress(Ele(i-j:i+j), [ones(size(UpLen(i-j:i+j))), UpLen(i-j:i+j)]);
        Slope(i) = abs(b(2));
    end
end

riverPro_f = [riverPro_Old_f, SmoothedEle, Localksn, Slope];

%% ͳ�Ƽ���
function [DW_r_f, DW_Ele_f, DW_Chi_f] = DW_test(Ele_f, Chi_f, z0_f, ksn_f)
DW_e = Ele_f - (z0_f + ksn_f .* Chi_f); % �в�
m = length(DW_e);
DW_e1 = DW_e(1:m-1);
DW_e2 = DW_e(2:m);
DW_r_f = sum(DW_e1.*DW_e2) / sqrt(sum(DW_e1.^2)) / sqrt(sum(DW_e2.^2));
DW_Chi_f = Chi_f(2:m) - DW_r_f .* Chi_f(1:m-1);
DW_Ele_f = Ele_f(2:m) - DW_r_f .* Ele_f(1:m-1);

%%
function [R2, a0, a1, std_a1] = Bi_Regress(Y, X)
% ���Իع� Y=a0+a1*X
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
    %�ع����ݲ���������ģ��
    R2 = 9999;
    a0 = 9999;
    a1 = 9999;
    std_a1 = 9999;
end

%%
