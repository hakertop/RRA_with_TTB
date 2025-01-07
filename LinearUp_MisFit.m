function MisFit=LinearUp_MisFit(chiDivision, workspace)
% MainQinL_Cal_MisFit�������򣬸�������ģ��ȼ������������Ԥ��chi-z��ʵ�ʶ�����ƫ��̶�
% ������chi�ռ������ļ���·�����ļ����к�����Ŀ
%%
format long; 
close all;

Path1=strcat(workspace, 'profile_'); 
%% ������Ҫ�ϲ��ĺ�������
% Path1='E:\NorthQinlin_GDEMV2 30M\Matlab_work\New_Inversion\basins_1-3\basin_1\chi_z45\RiverPro'; 
%���еĺ��������涼��������ļ�����
Path3='.txt';
% ��ȡprofile���ļ��� 
startRiver = 1;
endRiver = length(dir(strcat(workspace, "*.txt")));
sumRiver = endRiver - startRiver + 1; % ��������
AllRiverPro=cell(sumRiver,1); % ����cell�洢���еĺ������������
lengthChi=[];
for i=startRiver:endRiver
    % numRiver_i=input('������Ҫ�ϲ��ĺ�����ţ�');
    Path2=num2str(i);
    Path_i=strcat(Path1,Path2,Path3);
    RiverPro_i=textread(Path_i); % Y��X���꣬��Դ���롢�̡߳����򡢻�ˮ�����Chi����
    chiLen_i=RiverPro_i(:,7); Ele_i=RiverPro_i(:,4);
    minEle_i=min(Ele_i); Ele_i=Ele_i-minEle_i;
%     plot(chiLen_i,Ele_i,'m');hold on; % ���������Ѿ�����Ը߳�
    AllRiverPro{i-startRiver+1,1}=[chiLen_i,Ele_i];% ����cell�洢���еĺ������������
    lengthChi=[lengthChi;max(chiLen_i)]; % ��¼ÿһ�������� ���chhi
end

%% �������ϲ�Ϊһ��
% chiDivision=input('����ʱ��ڵ�����'); % ÿ��chi�ĳ���
MaxChi=max(lengthChi); 
lengthChi=[]; % Ѱ����ĺ�������

n_localksn=ceil(MaxChi/chiDivision); % ����ȡ����chi�ռ�ֶ���Ŀ��n_localksn=floor(MaxChi/chiDivision);
localksn=zeros(2,n_localksn); %Line1: ÿһ�ε�local_ksnֵ; Line2: ÿһ�μ����˼���ksn

localksn_mat=zeros(sumRiver,n_localksn); %% ÿһ�У���ʾһ��������ÿһ�У�������ͬ��Chi��

interChi=0:chiDivision:MaxChi; % �ϲ�֮���chi�ռ䣬�ȷּ��
if max(interChi) < MaxChi
    interChi=[interChi,MaxChi]; 
end  % interChi���ȱ�localksn��1

interEle=zeros(1,length(interChi)); % interEle=zeros(1,n_localksn+1);

for j=1:sumRiver
    RiverProj=AllRiverPro{j,1}; 
    chiLen_j=RiverProj(:,1); Ele_j=RiverProj(:,2); % ����������
    max_chiLen_j=max(chiLen_j);
    
    i=1; % ��1�Σ�chi���꣨1,1+1��
    while i <= n_localksn % ����ÿһ����chi-z�����շֶα�׼���б���chi�ռ����ֵ �Ƿ���С����
        if max_chiLen_j >= interChi(i+1) % �úӵ�chi�ռ䣬��ȫ���� ��(i,i+1)��
%             (interChi(i),interChi(i+1))�ε�local_ksn����дfunction����(chiLen_j,Ele_j,interChi(i),interChi(i+1))
            local_ksn_i=cal_local_ksn(chiLen_j,Ele_j,interChi(i),interChi(i+1));
            localksn(1,i)=localksn(1,i)+local_ksn_i; % ��һ���ۼ�ksn
            localksn(2,i)=localksn(2,i)+1; % ��һ���ۼƴ���
            i=i+1; % ��һ��
            localksn_mat(j,i)=local_ksn_i;
        elseif max_chiLen_j > interChi(i) && max_chiLen_j <= interChi(i+1) % ��(i,i+1)�ΰ����ֲ�
%             (interChi(i),max_chiLen_j)�ε�local_ksn
            local_ksn_i=cal_local_ksn(chiLen_j,Ele_j,interChi(i),max_chiLen_j);
            localksn(1,i)=localksn(1,i)+local_ksn_i;
            localksn(2,i)=localksn(2,i)+1;
            localksn_mat(j,i)=local_ksn_i;
            break;
        elseif max_chiLen_j < interChi(i) % ��(i,i+1)�β�����
            break;
        end
    end
end
% figure;plot(localksn(1,:),localksn(2,:),'o')
localksn=localksn(1,:)./localksn(2,:); 

for i=1:n_localksn
    interEle(i+1)=interEle(i)+(interChi(i+1)-interChi(i)).*localksn(i);
end

%% 
interChi=interChi';
interEle=interEle';
% plot(interChi,interEle,'g','linewidth',3);

% ������ƽ��chi-z�Ĳ�MisFit
MisFit=0;
for j=1:sumRiver
    RiverProj=AllRiverPro{j,1}; 
    chiLen_j=RiverProj(:,1); Ele_j=RiverProj(:,2); % ����������
    yi = interp1(interChi,interEle,chiLen_j);
    MisFit=MisFit+sqrt(mean((Ele_j-yi).^2));
end
MisFit=MisFit/sumRiver;




%% ��¼ÿһ��������ÿһ�������local_ksn
% % ����localksn_mat(sumRiver,n_localksn)
% Chi_ksnmat=chiDivision:chiDivision:chiDivision*n_localksn;
% mean_localksn=zeros(1,n_localksn);
% std_localksn=zeros(1,n_localksn);
% for kk=1:n_localksn
%     a=localksn_mat(:,kk);
%     b=a(a~=0);
%     mean_localksn(kk)=mean(b);
%     if length(b)>1 % ����һ��ֵ���Ͳ�����std��
%         std_localksn(kk)=std(b);
%     end
% end
% figure; plot(Chi_ksnmat,mean_localksn); hold on; title('mean_loacal_ksn');
% figure; plot(Chi_ksnmat,std_localksn); hold on; title('std_loacal_ksn')
% mean(std_localksn(std_localksn~=0)); % �Է����std�����ֵ



function local_ksn=cal_local_ksn(chiLen_j,Ele_j,interChi_i,interChi_i1)
% ����local_ksn���ӵ�����chi-z���������chiС��(interChi_i,interChi_i1)
% ���������������С��һ���ںӵ�����chi�ռ�֮��

m_len=length(chiLen_j);
chi_cal=[];Ele_cal=[];
for i=1:m_len
    if chiLen_j(i) >= interChi_i && chiLen_j(i) <= interChi_i1
        chi_cal=[chi_cal;chiLen_j(i)];
        Ele_cal=[Ele_cal;Ele_j(i)];
    elseif chiLen_j(i) > interChi_i1
        break;
    end
end

if length(chi_cal) < 2
    chi_cal=[];Ele_cal=[];
    chi_cal=[interChi_i;interChi_i1];
    Ele_cal=interp1(chiLen_j,Ele_j,chi_cal);
end

chi_cal=[ones(length(chi_cal),1),chi_cal];
localksn_ij=regress(Ele_cal,chi_cal);
local_ksn=localksn_ij(2);
















