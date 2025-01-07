function [outlet_x, outlet_y] = ttb_search_outlet_coordinate(dem, flow_accumulation, threshold, info)
    %% 通过TopoToolBox获取流域DEM的出水口点坐标
    % 作者：Zhaolong Dai
    %
    % 核心思想：
    %       出水口是指水流离开汇水区的点，这个点是汇水区边界上的最低点。 
    %       通常一条河流的汇水区没有其他的地表径流流入且只有惟一的一个出水点。
    %       因此，我们可以通过流域的DEM进行水文分析得到河网，提取河网栅格的高程数据，
    %       找到这部分高程中的最低点，即为我们所需的出水口点。
    %
    % 基本逻辑：
    %       1. 根据阈值提取出河网 stream_grid；
    %       2. 提取河网栅格的高程数据 stream_elevation ，并将河网之外的部分设为NaN；
    %       3. 找到 Stream_elevation 中高程最低的栅格点下标 outlet_linear_index；
    %       4. 根据 outlet_linear_index 可得到出水口的二维矩阵下标以及地理坐标
    %
    % 输入参数：
    %       dem - 流域DEM（GRIDobj）
    %       flow_accumulation - 流域流量数据（GRIDobj）
    %       threshold - 启动河流所需的最小上游区域（单位为像素）
    %       info - 是否显示具体信息，输入1则为显示
    %
    % 输出数据：
    %       outlet_x 出水口点横坐标
    %       outlet_y 出水口点纵坐标
    %
    % 出水口点是后续计算河流源头点、河流剖面信息、构造隆升历史等过程的基础。

    % 提取河网二值栅格图
    stream_grid = flow_accumulation > threshold;  
    
    % 提取河网栅格的高程数据，并将河网之外的部分设为 NaN
    stream_elevation = dem.Z;
    stream_elevation(~(stream_grid.Z)) = NaN;  
    
    % 找到高程最低点
    % 'omitan' Ignores all NaN values and returns the minimum of the non-NaN elements.
    [min_elevation, outlet_linear_index] = min(stream_elevation(:), [], 'omitnan');
    
    % 找到 linearIndex 对应的二维矩阵下标
    [outlet_row, outlet_col] = ind2sub(size(dem.Z), outlet_linear_index);
    
    % 获取出水口的地理坐标
    [outlet_x, outlet_y] = dem.ind2coord(outlet_linear_index);
    
    % 需要观察这里的流量是否正常
    if info == 1
        fprintf("出水口的二维矩阵下标为：(%d, %d) \n" + ...
            "出水口的地理坐标为：(%f, %f) \n" + ...
            "出水口的高程为：%f \n" + ...
            "出水口的流量为：%f", ...
            outlet_row, outlet_col, ...
            outlet_x, outlet_y, ...
            min_elevation, ...
            flow_accumulation.Z(outlet_row, outlet_col));
    end
end