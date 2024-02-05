function ff = binarySearch(f, d_min,bx,bz,r,x,z,n)
    % f: 目标函数
    % d: 阈值
    % 初始化界限
    lo = 0;
    hi = 1e4;
    % 允许的误差范围
    tolerance = 1e-3;

    % 二分查找
    while lo < hi
        mid = (lo + hi) / 2;
        if f(mid,bx,bz,r,x,z,n) >= d_min
            % 如果 f(mid) 满足条件，搜索左半部分
            hi = mid;
        else
            % 否则，搜索右半部分
            lo = mid;
        end

        % 检查是否达到了足够的精度
        if abs(hi - lo) < tolerance
            break;
        end
    end

    ff = lo;
end