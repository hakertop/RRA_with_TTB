function linear_uplift_uniformChi(profile_workspace, min_fit_range, K_values)

format long;
% close all;
% 均匀间隔 线性模拟山体隆升历史，去掉最小、最大高程
% Line13 河道高程剖面存储路径

%% 计算侵蚀系数
% K1=(1.08e-6); K2=(1.08e-6)-2.85e-7; K3=(1.08e-6)+2.85e-7;
% std=2.85e-7; 4.42e-7 侵蚀系数，所有的需要换算成国际单位制 m/a
% 1.47中,1.06东；1.10西，0.891中

%% 输入需要合并的河流个数
Path1 = strcat(profile_workspace, "profile_");
% 读取河流纵剖面的数量
profile_num = length(dir(strcat(profile_workspace, "*.txt")));
%所有的河流纵剖面都存在这个文件夹中
Path3 = '.txt';
% disp("河流共有 " + num2str(profile_num) + " 条河流");
% startRiver=input('输入开始计算的河流编号：startRiver='); % 这里从5开始  注意 Inyo山从5开始
% endRiver=input('输入最终计算的河流编号：endRiver=');

startRiver = 1;
endRiver = profile_num;

sumRiver = endRiver - startRiver + 1; %河流个数
AllRiverPro = cell(sumRiver, 1); % 利用cell存储所有的河流纵剖面矩阵
lengthChi = [];

fig1 = figure();

% 逐个读取河流剖面数据
% 遍历读取每个剖面文件，将其chi距离和相对高程存入AllRiverPro列表，并记录每条河流的chi最大值到lengthChi。
for i = startRiver:endRiver
    % numRiver_i=input('输入需要合并的河流编号：');
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % Y、X坐标，溯源距离、高程、流向、汇水面积、Chi距离
    chiLen_i = RiverPro_i(:, 7);
    Ele_i = RiverPro_i(:, 4);
    minEle_i = min(Ele_i);
    Ele_i = Ele_i - minEle_i;
    plot(chiLen_i, Ele_i, 'm');
    hold on; % 上述河流已经是相对高程
    AllRiverPro{i-startRiver+1, 1} = [chiLen_i, Ele_i]; % 利用cell存储所有的河流纵剖面矩阵
    lengthChi = [lengthChi; max(chiLen_i)]; % 记录每一河流长度 最大chhi
end

%% 将河流合并为一条

% 合并河流剖面数据，将所有河流的chi距离分成相等的间隔段，生成合并后的chi间隔interChi和对应的高程数组interEle。

% numRange_1=input('输入最小的chi值间隔'); % 1.3;   %input('输入时间节点间隔：'); % 每段chi的长度
numRange_1 = min_fit_range;

MaxChi = max(lengthChi);
lengthChi = []; % 寻找最长的河流长度

n_localksn = ceil(MaxChi/numRange_1); % 向上取整，chi空间分段数目；n_localksn=floor(MaxChi/numRange_1);
localksn = zeros(2, n_localksn); %Line1: 每一段的local_ksn值; Line2: 每一段计算了几个ksn

interChi = 0:numRange_1:MaxChi; % 合并之后的chi空间，等分间隔
if max(interChi) < MaxChi
    interChi = [interChi, MaxChi];
end % interChi长度比localksn多1

interEle = zeros(1, length(interChi)); % interEle=zeros(1,n_localksn+1);

% 计算各段的local_ksn
% 根据不同chi区间的分段长度计算每段的local_ksn，累计各段ksn值并计算其均值。
for j = 1:sumRiver
    RiverProj = AllRiverPro{j, 1};
    chiLen_j = RiverProj(:, 1);
    Ele_j = RiverProj(:, 2); % 都是列向量
    max_chiLen_j = max(chiLen_j);

    i = 1; % 第1段，chi坐标（1,1+1）
    while i <= n_localksn % 对于每一条河chi-z，按照分段标准，判别其chi空间最大值 是否在小段内
        if max_chiLen_j >= interChi(i+1) % 该河道chi空间，完全包含 第(i,i+1)段
            %             (interChi(i),interChi(i+1))段的local_ksn，编写function函数(chiLen_j,Ele_j,interChi(i),interChi(i+1))
            local_ksn_i = cal_local_ksn(chiLen_j, Ele_j, interChi(i), interChi(i+1));
            localksn(1, i) = localksn(1, i) + local_ksn_i; % 这一段累计ksn
            localksn(2, i) = localksn(2, i) + 1; % 这一段累计次数
            i = i + 1; % 下一段
        elseif max_chiLen_j > interChi(i) && max_chiLen_j <= interChi(i+1) % 第(i,i+1)段包含局部
            %             (interChi(i),max_chiLen_j)段的local_ksn
            local_ksn_i = cal_local_ksn(chiLen_j, Ele_j, interChi(i), max_chiLen_j);
            localksn(1, i) = localksn(1, i) + local_ksn_i;
            localksn(2, i) = localksn(2, i) + 1;
            break;
        elseif max_chiLen_j < interChi(i) % 第(i,i+1)段不包含
            break;
        end
    end
end
% figure;plot(localksn(1,:),localksn(2,:),'o')
localksn = localksn(1, :) ./ localksn(2, :);

% 绘制合并的河流纵剖面
for i = 1:n_localksn
    interEle(i+1) = interEle(i) + (interChi(i+1) - interChi(i)) .* localksn(i);
end

%%

% interChi=0:numRange_1:MaxChi;
% if max(interChi) < MaxChi; interChi=[interChi,MaxChi]; end
interChi = interChi';
interEle = interEle';

%% % 以上步骤得到合并后的河流纵剖面
plot(interChi, interEle, 'r', 'linewidth', 3);
hold on;
max_xlabel = ceil(ceil(max(interChi))/2) * 2;
xlim([0, max_xlabel]);
xlabel('Chi (m)');
ylabel('Relative elevation (m)');

saveas(fig1, strcat(profile_workspace, "../chi_plot.png"));

%% 计算无维度的隆升历史
len_path = length(interChi);
t_path = interChi;
U_path = (interEle(2:len_path) - interEle(1:len_path-1)) ./ (t_path(2:len_path) - t_path(1:len_path-1));

U_path = [U_path; U_path(length(U_path))];

fig2 = figure();
subplot(2, 1, 1);
stairs(t_path, U_path);
% axis([0,20,0,300]);
xlabel('t* (m)');
ylabel('U*');
% set(gca,'ytick',0:50:300);
% subplot(2,1,1);stairs(t_path(1:upLen-1),U_path);

%%
% upRate=(U_path).*K1*1e3; % 显示速率转换为mm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K1./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'r','LineWidth',2.5); hold on;
%
% upRate=(U_path).*K2*1e3; % 显示速率转换为mm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K2./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'g','LineWidth',1.5); hold on;
%
% upRate=(U_path).*K3*1e3; % 显示速率转换为mm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K3./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'b','LineWidth',1.5); hold on;

% 戴钊龙修改 使用for循环遍历所有K
colors = lines(length(K_values)); % 自动分配颜色
legend_labels = cell(1, length(K_values));
for i = 1:length(K_values)
    K = K_values(i);
    color = colors(i, :); % 获取第 i 个颜色

    % 计算有量纲的隆升速率和时间
    upRate = U_path .* K * 1e3; % 转换为 mm/a
    upTime = t_path ./ K ./ 1e6; % 转换为 Ma

    % 绘制阶梯图
    subplot(2, 1, 2);
    stairs(upTime, upRate, 'Color', color, 'LineWidth', 1.5);
    hold on;

    % 生成 legend 标签
    legend_labels{i} = sprintf('K%d = %.2e', i, K); % 自动生成 K 值标签
end

xlabel('t (Ma)');
ylabel('U (mm/a)');
legend(legend_labels, 'Location', 'best');
axis auto;
% xlim([0,20]);
% set(gca,'xtick',0:2:20);

saveas(fig2, strcat(profile_workspace, "../构造隆升历史.png"));

%%
function local_ksn = cal_local_ksn(chiLen_j, Ele_j, interChi_i, interChi_i1)
% 计算local_ksn，河道整体chi-z，待计算的chi小段(interChi_i,interChi_i1)
% 都是列向量，这个小段一定在河道整个chi空间之内

m_len = length(chiLen_j);
chi_cal = [];
Ele_cal = [];
for i = 1:m_len
    if chiLen_j(i) >= interChi_i && chiLen_j(i) <= interChi_i1
        chi_cal = [chi_cal; chiLen_j(i)];
        Ele_cal = [Ele_cal; Ele_j(i)];
    elseif chiLen_j(i) > interChi_i1
        break;
    end
end

if length(chi_cal) < 2
    chi_cal = [];
    Ele_cal = [];
    chi_cal = [interChi_i; interChi_i1];
    Ele_cal = interp1(chiLen_j, Ele_j, chi_cal);
end

chi_cal = [ones(length(chi_cal), 1), chi_cal];
localksn_ij = regress(Ele_cal, chi_cal);
local_ksn = localksn_ij(2);
