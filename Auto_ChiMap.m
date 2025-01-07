function ChiMap=Auto_ChiMap(DEMFilePath, FlowAcc,FlowDir,LeftX,DownY,cellsize,startX,startY,Ac,mVSn,workspace)
% 计算chi-map,可以生成矢量文件（willitte 2014 Science） 
% - DEMFilePath         DEM文件路径，主要用于生成结果矢量的空间坐标系
% - FlowAcc             汇水面积矩阵； 
% - FlowDir             流向矩阵；
% - LeftX, DownY        分别为栅格图像起始左下点坐标 
% - cellsize            栅格图像步长； 
% - startX, startY      分别为搜索起始点坐标；
% - Ac                  汇水面积阈值 m^2； 
% - mVSn                参考凹度
% - workspace           输出空间

format long;
[m,n]=size(FlowAcc);
ChiMap=zeros(m,n);
UpY=DownY+cellsize*m;
Ac=floor(Ac/(cellsize^2)); % Ac位汇水面积临界值

len_startX=length(startX);
stepPoint=[];

Sfea_ChiMap=struct(); %生成polyline的结构体
Sfea_num=1; % chimap矢量文件的累积增加

visited = false(m, n);

for ilen=1:len_startX
    % 起始搜索点的行、列号
%     startRow=floor((UpY-startY(ilen))/cellsize)+1; 
%     startCol=floor((startX(ilen)-LeftX)/cellsize)+1;
    startRow = m - floor((startY(ilen) - DownY)/cellsize);
    startCol = floor((startX(ilen) - LeftX)/cellsize) + 1;

    ChiMap(startRow,startCol)=0;
    stepPoint=[stepPoint;[startRow,startCol,FlowAcc(startRow,startCol),0,0]];
    empty_stepPoint=1;
    
    % 搜索点行、列号；该点的汇水面积；该点汇入区域的行、列号（0,0表示无汇入）
    while empty_stepPoint
        startRow=stepPoint(1,1);
        startCol=stepPoint(1,2);

        if visited(startRow, startCol)
            stepPoint(1, :) = []; % % 每一轮搜索完成后，原来的起始点 出栈, 第一个点为上一步的搜索点
            empty_stepPoint = ~isempty(stepPoint);
            continue;
        end
        visited(startRow, startCol) = true;

        
        if startCol<n && FlowDir(startRow,startCol+1)==2^4 && FlowAcc(startRow,startCol+1)>Ac
            % 该点 右侧 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow,startCol+1,FlowAcc(startRow,startCol+1),startRow,startCol]];
            ChiMap(startRow,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow];   CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow,startCol+1);
            Sfea_num=Sfea_num+1;
        end  
        if startRow<m && startCol<n && FlowDir(startRow+1,startCol+1)==2^5 && FlowAcc(startRow+1,startCol+1)>Ac
            % 该点 右下方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow+1,startCol+1,FlowAcc(startRow+1,startCol+1),startRow,startCol]];
            ChiMap(startRow+1,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
             
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol+1);
            Sfea_num=Sfea_num+1;

        end
        if startRow<m && FlowDir(startRow+1,startCol)==2^6 && FlowAcc(startRow+1,startCol)>Ac
            % 该点 下方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow+1,startCol,FlowAcc(startRow+1,startCol),startRow,startCol]];
            ChiMap(startRow+1,startCol)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol];   CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol);
            Sfea_num=Sfea_num+1;

        end
        if startRow<m && startCol>1 && FlowDir(startRow+1,startCol-1)==2^7 && FlowAcc(startRow+1,startCol-1)>Ac
            % 该点 左下方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow+1,startCol-1,FlowAcc(startRow+1,startCol-1),startRow,startCol]];
            ChiMap(startRow+1,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol-1);
            Sfea_num=Sfea_num+1;
        end
        if startCol>1 && FlowDir(startRow,startCol-1)==2^0 && FlowAcc(startRow,startCol-1)>Ac
            % 该点 左方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow,startCol-1,FlowAcc(startRow,startCol-1),startRow,startCol]];
            ChiMap(startRow,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow];   CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow,startCol-1);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && startCol>1 && FlowDir(startRow-1,startCol-1)==2^1 && FlowAcc(startRow-1,startCol-1)>Ac
            % 该点 左上方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow-1,startCol-1,FlowAcc(startRow-1,startCol-1),startRow,startCol]];
            ChiMap(startRow-1,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol-1);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && FlowDir(startRow-1,startCol)==2^2 && FlowAcc(startRow-1,startCol)>Ac
            % 该点 上方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow-1,startCol,FlowAcc(startRow-1,startCol),startRow,startCol]];
            ChiMap(startRow-1,startCol)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol];   CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && startCol<n && FlowDir(startRow-1,startCol+1)==2^3 && FlowAcc(startRow-1,startCol+1)>Ac
            % 该点 右上方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow-1,startCol+1,FlowAcc(startRow-1,startCol+1),startRow,startCol]];
            ChiMap(startRow-1,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol+1);
            Sfea_num=Sfea_num+1;

        end
        % 每一轮搜索完成后，原来的起始点 出栈
        stepPoint(1,:)=[]; % 第一个点为上一步的搜索点
        empty_stepPoint=~isempty(stepPoint);
    end
end

% 输出矢量信息
output_shapefile_path = strcat(workspace, "ChiMap.shp");
shapewrite(Sfea_ChiMap, output_shapefile_path);

% Zhaolong Dai 添加：补充空间坐标系信息
info = georasterinfo(DEMFilePath);
p = info.CoordinateReferenceSystem;
wkt = wktstring(p, 'Version', "wkt1");
fid = fopen(strcat(workspace, "ChiMap.prj"), "w");
fprintf(fid, "%s", wkt);
fclose(fid);

output_ras_path = strcat(workspace, 'Ras_ChiMap.txt');
fid=fopen(output_ras_path,'w');
for j=1:m
    for k=1:n-1
        fprintf(fid,'%f\t',ChiMap(j,k));
    end
    fprintf(fid,'%f\n',ChiMap(j,n));
end
fclose(fid);























