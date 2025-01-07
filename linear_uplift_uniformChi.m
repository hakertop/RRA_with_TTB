function linear_uplift_uniformChi(profile_workspace, min_fit_range, K_values)

format long;
% close all;
% ���ȼ�� ����ģ��ɽ��¡����ʷ��ȥ����С�����߳�
% Line13 �ӵ��߳�����洢·��

%% ������ʴϵ��
% K1=(1.08e-6); K2=(1.08e-6)-2.85e-7; K3=(1.08e-6)+2.85e-7;
% std=2.85e-7; 4.42e-7 ��ʴϵ�������е���Ҫ����ɹ��ʵ�λ�� m/a
% 1.47��,1.06����1.10����0.891��

%% ������Ҫ�ϲ��ĺ�������
Path1 = strcat(profile_workspace, "profile_");
% ��ȡ���������������
profile_num = length(dir(strcat(profile_workspace, "*.txt")));
%���еĺ��������涼��������ļ�����
Path3 = '.txt';
% disp("�������� " + num2str(profile_num) + " ������");
% startRiver=input('���뿪ʼ����ĺ�����ţ�startRiver='); % �����5��ʼ  ע�� Inyoɽ��5��ʼ
% endRiver=input('�������ռ���ĺ�����ţ�endRiver=');

startRiver = 1;
endRiver = profile_num;

sumRiver = endRiver - startRiver + 1; %��������
AllRiverPro = cell(sumRiver, 1); % ����cell�洢���еĺ������������
lengthChi = [];

fig1 = figure();

% �����ȡ������������
% ������ȡÿ�������ļ�������chi�������Ը̴߳���AllRiverPro�б�����¼ÿ��������chi���ֵ��lengthChi��
for i = startRiver:endRiver
    % numRiver_i=input('������Ҫ�ϲ��ĺ�����ţ�');
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % Y��X���꣬��Դ���롢�̡߳����򡢻�ˮ�����Chi����
    chiLen_i = RiverPro_i(:, 7);
    Ele_i = RiverPro_i(:, 4);
    minEle_i = min(Ele_i);
    Ele_i = Ele_i - minEle_i;
    plot(chiLen_i, Ele_i, 'm');
    hold on; % ���������Ѿ�����Ը߳�
    AllRiverPro{i-startRiver+1, 1} = [chiLen_i, Ele_i]; % ����cell�洢���еĺ������������
    lengthChi = [lengthChi; max(chiLen_i)]; % ��¼ÿһ�������� ���chhi
end

%% �������ϲ�Ϊһ��

% �ϲ������������ݣ������к�����chi����ֳ���ȵļ���Σ����ɺϲ����chi���interChi�Ͷ�Ӧ�ĸ߳�����interEle��

% numRange_1=input('������С��chiֵ���'); % 1.3;   %input('����ʱ��ڵ�����'); % ÿ��chi�ĳ���
numRange_1 = min_fit_range;

MaxChi = max(lengthChi);
lengthChi = []; % Ѱ����ĺ�������

n_localksn = ceil(MaxChi/numRange_1); % ����ȡ����chi�ռ�ֶ���Ŀ��n_localksn=floor(MaxChi/numRange_1);
localksn = zeros(2, n_localksn); %Line1: ÿһ�ε�local_ksnֵ; Line2: ÿһ�μ����˼���ksn

interChi = 0:numRange_1:MaxChi; % �ϲ�֮���chi�ռ䣬�ȷּ��
if max(interChi) < MaxChi
    interChi = [interChi, MaxChi];
end % interChi���ȱ�localksn��1

interEle = zeros(1, length(interChi)); % interEle=zeros(1,n_localksn+1);

% ������ε�local_ksn
% ���ݲ�ͬchi����ķֶγ��ȼ���ÿ�ε�local_ksn���ۼƸ���ksnֵ���������ֵ��
for j = 1:sumRiver
    RiverProj = AllRiverPro{j, 1};
    chiLen_j = RiverProj(:, 1);
    Ele_j = RiverProj(:, 2); % ����������
    max_chiLen_j = max(chiLen_j);

    i = 1; % ��1�Σ�chi���꣨1,1+1��
    while i <= n_localksn % ����ÿһ����chi-z�����շֶα�׼���б���chi�ռ����ֵ �Ƿ���С����
        if max_chiLen_j >= interChi(i+1) % �úӵ�chi�ռ䣬��ȫ���� ��(i,i+1)��
            %             (interChi(i),interChi(i+1))�ε�local_ksn����дfunction����(chiLen_j,Ele_j,interChi(i),interChi(i+1))
            local_ksn_i = cal_local_ksn(chiLen_j, Ele_j, interChi(i), interChi(i+1));
            localksn(1, i) = localksn(1, i) + local_ksn_i; % ��һ���ۼ�ksn
            localksn(2, i) = localksn(2, i) + 1; % ��һ���ۼƴ���
            i = i + 1; % ��һ��
        elseif max_chiLen_j > interChi(i) && max_chiLen_j <= interChi(i+1) % ��(i,i+1)�ΰ����ֲ�
            %             (interChi(i),max_chiLen_j)�ε�local_ksn
            local_ksn_i = cal_local_ksn(chiLen_j, Ele_j, interChi(i), max_chiLen_j);
            localksn(1, i) = localksn(1, i) + local_ksn_i;
            localksn(2, i) = localksn(2, i) + 1;
            break;
        elseif max_chiLen_j < interChi(i) % ��(i,i+1)�β�����
            break;
        end
    end
end
% figure;plot(localksn(1,:),localksn(2,:),'o')
localksn = localksn(1, :) ./ localksn(2, :);

% ���ƺϲ��ĺ���������
for i = 1:n_localksn
    interEle(i+1) = interEle(i) + (interChi(i+1) - interChi(i)) .* localksn(i);
end

%%

% interChi=0:numRange_1:MaxChi;
% if max(interChi) < MaxChi; interChi=[interChi,MaxChi]; end
interChi = interChi';
interEle = interEle';

%% % ���ϲ���õ��ϲ���ĺ���������
plot(interChi, interEle, 'r', 'linewidth', 3);
hold on;
max_xlabel = ceil(ceil(max(interChi))/2) * 2;
xlim([0, max_xlabel]);
xlabel('Chi (m)');
ylabel('Relative elevation (m)');

saveas(fig1, strcat(profile_workspace, "../chi_plot.png"));

%% ������ά�ȵ�¡����ʷ
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
% upRate=(U_path).*K1*1e3; % ��ʾ����ת��Ϊmm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K1./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'r','LineWidth',2.5); hold on;
%
% upRate=(U_path).*K2*1e3; % ��ʾ����ת��Ϊmm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K2./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'g','LineWidth',1.5); hold on;
%
% upRate=(U_path).*K3*1e3; % ��ʾ����ת��Ϊmm/a
% % upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
% upTime=t_path./K3./1e6;
% subplot(2,1,2);stairs(upTime,upRate,'b','LineWidth',1.5); hold on;

% �������޸� ʹ��forѭ����������K
colors = lines(length(K_values)); % �Զ�������ɫ
legend_labels = cell(1, length(K_values));
for i = 1:length(K_values)
    K = K_values(i);
    color = colors(i, :); % ��ȡ�� i ����ɫ

    % ���������ٵ�¡�����ʺ�ʱ��
    upRate = U_path .* K * 1e3; % ת��Ϊ mm/a
    upTime = t_path ./ K ./ 1e6; % ת��Ϊ Ma

    % ���ƽ���ͼ
    subplot(2, 1, 2);
    stairs(upTime, upRate, 'Color', color, 'LineWidth', 1.5);
    hold on;

    % ���� legend ��ǩ
    legend_labels{i} = sprintf('K%d = %.2e', i, K); % �Զ����� K ֵ��ǩ
end

xlabel('t (Ma)');
ylabel('U (mm/a)');
legend(legend_labels, 'Location', 'best');
axis auto;
% xlim([0,20]);
% set(gca,'xtick',0:2:20);

saveas(fig2, strcat(profile_workspace, "../����¡����ʷ.png"));

%%
function local_ksn = cal_local_ksn(chiLen_j, Ele_j, interChi_i, interChi_i1)
% ����local_ksn���ӵ�����chi-z���������chiС��(interChi_i,interChi_i1)
% ���������������С��һ���ںӵ�����chi�ռ�֮��

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
