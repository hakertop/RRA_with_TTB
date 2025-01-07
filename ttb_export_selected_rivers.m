function selected_river_head = ttb_export_selected_rivers(workspace, dem, stream)
    %% 基于TopoToolBox展示流域信息，挑选出所需要分析的河流的源头点坐标数据后导出
    % 输入数据：
    %   - workspace
    %   - dem
    %   - stream
    %
    % 输出数据：
    %   - selected_river_head
    %   - workspace/selected_river_head.txt
    %
    % 注意：要求之前分析过河流源头点数据，在workspace中存有river_head.txt文件

    % 加载前面得到的河流源头数据
    river_head_path = strcat(workspace, "river_head.txt");
    if ~isfile(river_head_path)
        error('源头点文件不存在：%s', river_head_path);
    end
    river_head = load(river_head_path);
    river_head_x = river_head(:, 1);
    river_head_y = river_head(:, 2);

    % 可视化河流及其源头点
    figure;
    imageschs(dem);
    hold on;
    plot(stream, 'LineWidth', 1, 'Color', 'b');                                         % 绘制河流
    scatter(river_head_x, river_head_y, 50, 'g', 'filled', 'MarkerEdgeColor', 'k');     % 绘制源头点
    title('河流及其源头点（左键点选选取，Enter结束选取）');
    xlabel('X 坐标');
    ylabel('Y 坐标');

    % 用户选择需要分析的河道源头点
    disp('请点击选择感兴趣的河流源头点（左键点选选取，Enter结束选取）');
    selected_river_head = [];
    selected_indices = false(size(river_head_x)); % 用于记录已选择的点索引

    while true
        % 获取用户点击的点
        [x, y, button] = ginput(1);

        if isempty(button) || button ~= 1 % 右键或完成选择
            k = msgbox('选择完毕，将退出！');

            pause(2)
            if ishandle(k)
                delete(k)
            end

            break;
        end

        % 查找距离用户点击点最近的源头点
        distances = sqrt((river_head_x - x).^2 + (river_head_y - y).^2);
        [min_distance, idx] = min(distances);

        if min_distance < 5 * dem.cellsize % 如果点击点在阈值范围内
            if ~selected_indices(idx) % 如果尚未选择
                % 标记为已选择，并改变颜色
                selected_indices(idx) = true;
                scatter(river_head_x(idx), river_head_y(idx), 70, 'r', 'filled', 'MarkerEdgeColor', 'k');
                selected_river_head = [selected_river_head; river_head_x(idx), river_head_y(idx)];
            else
                % 已选择的点被再次点击，提示用户
                disp('该点已被选择！');
            end
        else
            disp('点击位置未匹配到任何源头点，请重新点击。');
        end
    end
    hold off

    % 导出所选择的河流源头点到文件
    if isempty(selected_river_head)
        warning('未选择任何河流源头点。');
    else
        selected_river_head_path = fullfile(workspace, 'selected_river_head.txt');
        writematrix(selected_river_head, selected_river_head_path, 'Delimiter', '\t');
        fprintf('已导出所选河流源头点到文件：%s\n', selected_river_head_path);
    end

    close
end

