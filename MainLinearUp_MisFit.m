function minMisFitRange = MainLinearUp_MisFit(workspace)
profile_workspace = strcat(workspace, "profile\");
outputFilePath = strcat(workspace, 'MisFit_chiz_river.txt');
fid = fopen(outputFilePath, 'w');

% 找到高程偏差最小的
minMisFit = 999999;

for i = 1:99
    numRange_1 = i / 10; % chi值间隔
    MisFit = LinearUp_MisFit(numRange_1, profile_workspace); % 直接调用LinearUp_MisFit函数
    fprintf(fid, '%f\t', numRange_1);
    fprintf(fid, '%f\n', MisFit);

    if MisFit < minMisFit
        minMisFitRange = numRange_1;
        minMisFit = MisFit;
    end
end
fclose(fid);
end