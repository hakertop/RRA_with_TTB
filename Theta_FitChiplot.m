function min_theta = Theta_FitChiplot(workspace, profile_num)
% Mudd et al (2018)方法计算最佳凹度θ Line10 河流高程剖面存储路径；Line54 存储θ-高程偏差图

%%
% 初始化参数
format long;
close all;
theta = 0.1:0.05:0.9;
len_theta = length(theta);

%% 读取所有河流，并计算不同θ对应的chi-z
Path1 = strcat(workspace, "profile_");
%所有的河流纵剖面都存在这个文件夹中
Path3 = '.txt';
startRiver = 1; % input('输入开始计算的河流编号：startRiver='); % 这里从5开始  注意 Inyo山从5开始
% endRiver=input('输入最终计算的河流编号：endRiver=');
endRiver = profile_num;
sumRiver = endRiver - startRiver + 1; %河流个数
AllRiverPro = cell(sumRiver, 1); % 利用cell存储所有的河流纵剖面矩阵
lengthChi = [];

for i = startRiver:endRiver
    % numRiver_i=input('输入需要合并的河流编号：');
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % Y、X坐标，溯源距离、高程、流向、汇水面积、Chi距离
    Ele_i = RiverPro_i(:, 4);
    Area_i = RiverPro_i(:, 6);
    upLen_i = RiverPro_i(:, 3);

    chiLen_i = cal_chiLen_i(Area_i, upLen_i, theta); % θ为一行θ值，每个θ对应一列chiLen_i，所以，chiLen_i是一个矩阵

    minEle_i = min(Ele_i);
    Ele_i = Ele_i - minEle_i;
    %     plot(chiLen_i,Ele_i,'m');hold on; % 上述河流已经是相对高程
    AllRiverPro{i-startRiver+1, 1} = [chiLen_i, Ele_i]; % 利用cell存储所有的河流纵剖面矩阵
    %     lengthChi=[lengthChi;max(chiLen_i)]; % 记录每一河流长度 最大chhi
end

%% 对cell里每条河流，选定一列Chi（同一个θ）和Ele，单拿出来，寻找最大Chi对应的河道，然后每条河流依次计算misFit。
MisFit_i = zeros(1, len_theta);
for i = 1:len_theta
    RiverPro_Chosen = cell(sumRiver, 1);
    for j = 1:sumRiver
        ChiZ_ij = AllRiverPro{j, 1};
        RiverPro_Chosen{j, 1} = [ChiZ_ij(:, i), ChiZ_ij(:, len_theta+1)]; % 对应θ下的Chi-z
    end
    MisFit_i(i) = co_stem(RiverPro_Chosen); %RiverPro_Chosen是cell，每个元素都是chi-z
end
fig = figure;
plot(theta, MisFit_i);
hold on;
xlabel('theta');
ylabel('Elevation scatter (m)');

[min_MisFit_i, num_MisFit_i] = min(MisFit_i);
min_theta = theta(num_MisFit_i);
title(strcat(num2str(min_MisFit_i), ',theta=', num2str(min_theta)));
max_MisFit_i = max(MisFit_i);
ylim_down = floor(min_MisFit_i/10) * 10;
ylim_up = ceil(max_MisFit_i/10) * 10;
ylim([ylim_down, ylim_up]);

print theta_scatter_1-3.eps -depsc2 -r300; % **************************

% 保存θ和最小高程scatter
saveas(fig, strcat(workspace, "../", "figure.jpg"), "jpeg");
fid = fopen(strcat(workspace, "../", 'theta-min_Scatter_rather.txt'), 'w');
fprintf(fid, '%f\t', min_MisFit_i);
fprintf(fid, '%f\n', min_theta);
fclose(fid);

%%  计算不同θ值对应的Chi
function chiLen_i = cal_chiLen_i(Area_new, upLen_new, concavity)

chiLen_i = [];

len_con = length(concavity);
mLen_new = length(Area_new);
for i = 1:len_con
    Chi = [];
    sum_i = 0;
    Chi = [Chi; sum_i];
    for j = 2:mLen_new
        sum_i = sum_i + (1 / Area_new(j-1)).^(concavity(i)) * (upLen_new(j) - upLen_new(j-1));
        Chi = [Chi; sum_i]; % Chi是列向量
    end
    chiLen_i = [chiLen_i, Chi];
end

%% 计算Chi-z的misFit,最接近主干道
function MisFit_i = co_stem(RiverPro_Chosen)
%RiverPro_Chosen是cell，每个元素都是chi-z

sumRiver = length(RiverPro_Chosen);

% 寻找stem，即Chi最长的河段
max_Chi = 0;
for i = 1:sumRiver
    RiverPro_i = RiverPro_Chosen{i, 1};
    Chi_i = RiverPro_i(:, 1);
    max_Chi_i = max(Chi_i);
    if max_Chi_i > max_Chi
        max_Chi = max_Chi_i;
        num_max_Chi = i;
    end
end
RiverPro_Stem = RiverPro_Chosen{num_max_Chi, 1};
Stem_Chi = RiverPro_Stem(:, 1);
Stem_z = RiverPro_Stem(:, 2);

% 计算MisFit_i
num_node = 0;
MisFit_i = 0;
for j = 1:sumRiver
    if j ~= num_max_Chi
        RiverProj = RiverPro_Chosen{j, 1};
        Chi_j = RiverProj(:, 1);
        z_j = RiverProj(:, 2); % 都是列向量
        z_j_recal = interp1(Stem_Chi, Stem_z, Chi_j);
        MisFit_i = MisFit_i + sum((z_j - z_j_recal).^2);
        num_node = num_node + length(Chi_j);
    end
end
MisFit_i = sqrt(MisFit_i/num_node);

%%
