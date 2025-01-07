function ttb_export_selected_basins(dem_file_path, output_workspace)
    %% 流域分析，选取所需的流域导出为 TIF
    %
    % 输入参数：
    %   dem_file_path - 研究区 DEM 路径 
    %   output_workspace - 输出文件夹
    %
    % 输出数据：output_workspace 中生成
    %   basins_id_list - 流域id列表 
    %   basins\basin_(id).tif - 流域范围DEM，其中id为其编号

    % 检查输入文件和输出文件夹
    if ~isfile(dem_file_path)
        error('DEM 文件路径无效：%s', dem_file_path);
    end
    if ~isfolder(output_workspace)
        mkdir(output_workspace);
    end
    output_workspace = strcat(output_workspace, "basins\");
    mkdir(output_workspace);
    
    % 加载 DEM 并进行填洼、高斯滤波
    dem = fillsinks(GRIDobj(dem_file_path));
    dem.Z = imgaussfilt(dem.Z, 2);
    
    % 流向分析
    flow_direction = FLOWobj(dem);
    
    % 流域分析
    basins = drainagebasins(flow_direction);
    [~, L_x, L_y] = GRIDobj2polygon(basins);
    
    % 可视化
    figure;
    imageschs(dem, basins, 'colormap', lines);
    hold on;
    plot(L_x, L_y, 'Color', 'w', 'LineWidth', 2);
    hold off;
    title('流域分区（单击右键选择感兴趣的流域）');
    
    % 用户交互选择流域
    disp('请点击感兴趣的流域（单击右键完成选择）');
    selected_basins = ginput();
    % 将地理坐标转换为矩阵索引
    [row, col] = dem.coord2sub(selected_basins(:, 1), selected_basins(:, 2));
    % 根据矩阵索引提取流域 ID
    selected_ids = unique(basins.Z(sub2ind(size(basins.Z), row, col)));
    selected_ids(selected_ids == 0) = []; % 排除背景值
    
    % id列表存为txt
    fid = fopen(strcat(output_workspace, "../basins_id_list.txt"), "w+");

    % 导出选定流域为 TIF
    for id = selected_ids'
        fprintf(fid, "%d\n", id);
        % 生成该流域的文件夹
        basin_output_file = strcat(output_workspace, "\basin_", num2str(id), "\");
        mkdir(basin_output_file);

        mask = basins == id;
        basin_elevation = dem;
        basin_elevation.Z(~mask.Z) = NaN; % 非流域区域设置为 NaN
        bound = ttb_extract_valid_bound(basin_elevation);
        basin_elevation = crop(basin_elevation, [bound(1), bound(3)], [bound(2), bound(4)]);
        output_file = fullfile(basin_output_file, sprintf('basin_%d.tif', id));
        basin_elevation.GRIDobj2geotiff(output_file);
         % 展示导出的DEM
         
         figure;
         imageschs(basin_elevation);
         title(strcat("流域 ", num2str(id), " DEM"));
        fprintf('流域 %d 导出到 %s\n', id, output_file);
    end
    
    fclose(fid);
    
    disp('流域导出完成！');
end

function rectangle = ttb_extract_valid_bound(gridObj)
    %% 提取TopoToolBox GRIDobj对象的有效范围（忽略四边NaN值）
    %
    % 输入参数：
    %   - gridObj: GRIDobj 对象
    %
    % 输出参数：
    %   - geoBounds: 包含有效范围的地理坐标，格式为 [minX, minY, maxX, maxY]
    %
    % 主要是由于TopoToolBox生成的GRIDobj周边包含了大量的NaN无法去除，因此使用此函数获取有效范围
    % 得到有效范围后可以用来限制可视化时的范围

    % 检查输入类型
    if ~isa(gridObj, 'GRIDobj')
        error('输入必须是GRIDobj对象');
    end

    % 提取有效区域的逻辑索引
    validMask = ~isnan(gridObj.Z);

    % 如果整个 GRIDobj 都是 NaN，直接返回空
    if all(~validMask, 'all')
        warning('GRIDobj 对象中没有有效数据');
        rectangle = [];
        return;
    end

    % 找到包含有效数据的行和列范围
    rowIndices = find(any(validMask, 2));
    colIndices = find(any(validMask, 1));

    minRow = rowIndices(1);
    maxRow = rowIndices(end);
    minCol = colIndices(1);
    maxCol = colIndices(end);

    % 转换为地理坐标
    [minX, minY] = gridObj.sub2coord(maxRow, minCol);
    [maxX, maxY] = gridObj.sub2coord(minRow, maxCol);

    % 返回地理范围
    rectangle = [minX, minY, maxX, maxY];
end
