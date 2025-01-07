function visualize_shapefile(dem, shapefile_name, field_name)
    % visualize_shapefile - 可视化 shapefile 文件，支持 chi-map 和 ksn-map
    %
    % 输入参数:
    %   dem            - dem（GRIDobj）
    %   shapefile_name - shapefile 文件名（字符串）
    %   field_name     - 用于颜色映射的字段名（字符串，例如 'ChiValue' 或 'ksn'）
    %
    % 示例:
    %   visualize_shapefile('chimap.shp', 'ChiValue');
    %   visualize_shapefile('ksnmap.shp', 'ksn');

    % 读取 shapefile 文件
    shapefile_data = shaperead(shapefile_name);
    
    % 检查字段是否存在
    if ~isfield(shapefile_data, field_name)
        error('字段 "%s" 不存在于 shapefile 中，请检查字段名。', field_name);
    end
    
    % 提取用于颜色映射的字段值
    field_values = [shapefile_data.(field_name)];
    
    % 归一化字段值到 [0, 1]
    field_values_norm = (field_values - min(field_values)) / (max(field_values) - min(field_values));
    
    % 创建颜色映射矩阵
    color_map = parula(256); % 使用 Parula 色彩方案
    colors = interp1(linspace(0, 1, size(color_map, 1)), color_map, field_values_norm);

    % 绘制图形
    figure;
    imageschs(dem, [], 'colormap', [.9 .9 .9], 'colorbar',false);
    hold on;
    for i = 1:length(shapefile_data)
        % 获取每个 shapefile 对象的几何形状
        if strcmp(shapefile_data(i).Geometry, 'Line') || strcmp(shapefile_data(i).Geometry, 'Polyline')
            % 线状对象
            plot(shapefile_data(i).X, shapefile_data(i).Y, 'Color', colors(i, :), 'LineWidth', 1.5);
        elseif strcmp(shapefile_data(i).Geometry, 'Polygon')
            % 多边形对象
            fill(shapefile_data(i).X, shapefile_data(i).Y, colors(i, :), 'EdgeColor', 'none');
        end
    end
    hold off;
    
    % 添加颜色条
    colormap(color_map);
    colorbar;
    title(sprintf('Shapefile Visualized by %s', field_name));
end