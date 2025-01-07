function min_theta = Theta_FitChiplot(workspace, profile_num)
% Mudd et al (2018)����������Ѱ��Ȧ� Line10 �����߳�����洢·����Line54 �洢��-�߳�ƫ��ͼ

%%
% ��ʼ������
format long;
close all;
theta = 0.1:0.05:0.9;
len_theta = length(theta);

%% ��ȡ���к����������㲻ͬ�ȶ�Ӧ��chi-z
Path1 = strcat(workspace, "profile_");
%���еĺ��������涼��������ļ�����
Path3 = '.txt';
startRiver = 1; % input('���뿪ʼ����ĺ�����ţ�startRiver='); % �����5��ʼ  ע�� Inyoɽ��5��ʼ
% endRiver=input('�������ռ���ĺ�����ţ�endRiver=');
endRiver = profile_num;
sumRiver = endRiver - startRiver + 1; %��������
AllRiverPro = cell(sumRiver, 1); % ����cell�洢���еĺ������������
lengthChi = [];

for i = startRiver:endRiver
    % numRiver_i=input('������Ҫ�ϲ��ĺ�����ţ�');
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % Y��X���꣬��Դ���롢�̡߳����򡢻�ˮ�����Chi����
    Ele_i = RiverPro_i(:, 4);
    Area_i = RiverPro_i(:, 6);
    upLen_i = RiverPro_i(:, 3);

    chiLen_i = cal_chiLen_i(Area_i, upLen_i, theta); % ��Ϊһ�Ц�ֵ��ÿ���ȶ�Ӧһ��chiLen_i�����ԣ�chiLen_i��һ������

    minEle_i = min(Ele_i);
    Ele_i = Ele_i - minEle_i;
    %     plot(chiLen_i,Ele_i,'m');hold on; % ���������Ѿ�����Ը߳�
    AllRiverPro{i-startRiver+1, 1} = [chiLen_i, Ele_i]; % ����cell�洢���еĺ������������
    %     lengthChi=[lengthChi;max(chiLen_i)]; % ��¼ÿһ�������� ���chhi
end

%% ��cell��ÿ��������ѡ��һ��Chi��ͬһ���ȣ���Ele�����ó�����Ѱ�����Chi��Ӧ�ĺӵ���Ȼ��ÿ���������μ���misFit��
MisFit_i = zeros(1, len_theta);
for i = 1:len_theta
    RiverPro_Chosen = cell(sumRiver, 1);
    for j = 1:sumRiver
        ChiZ_ij = AllRiverPro{j, 1};
        RiverPro_Chosen{j, 1} = [ChiZ_ij(:, i), ChiZ_ij(:, len_theta+1)]; % ��Ӧ���µ�Chi-z
    end
    MisFit_i(i) = co_stem(RiverPro_Chosen); %RiverPro_Chosen��cell��ÿ��Ԫ�ض���chi-z
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

% ����Ⱥ���С�߳�scatter
saveas(fig, strcat(workspace, "../", "figure.jpg"), "jpeg");
fid = fopen(strcat(workspace, "../", 'theta-min_Scatter_rather.txt'), 'w');
fprintf(fid, '%f\t', min_MisFit_i);
fprintf(fid, '%f\n', min_theta);
fclose(fid);

%%  ���㲻ͬ��ֵ��Ӧ��Chi
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
        Chi = [Chi; sum_i]; % Chi��������
    end
    chiLen_i = [chiLen_i, Chi];
end

%% ����Chi-z��misFit,��ӽ����ɵ�
function MisFit_i = co_stem(RiverPro_Chosen)
%RiverPro_Chosen��cell��ÿ��Ԫ�ض���chi-z

sumRiver = length(RiverPro_Chosen);

% Ѱ��stem����Chi��ĺӶ�
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

% ����MisFit_i
num_node = 0;
MisFit_i = 0;
for j = 1:sumRiver
    if j ~= num_max_Chi
        RiverProj = RiverPro_Chosen{j, 1};
        Chi_j = RiverProj(:, 1);
        z_j = RiverProj(:, 2); % ����������
        z_j_recal = interp1(Stem_Chi, Stem_z, Chi_j);
        MisFit_i = MisFit_i + sum((z_j - z_j_recal).^2);
        num_node = num_node + length(Chi_j);
    end
end
MisFit_i = sqrt(MisFit_i/num_node);

%%
