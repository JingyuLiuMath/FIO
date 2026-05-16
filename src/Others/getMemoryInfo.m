function mem = getMemoryInfo()
    %GETMEMORYINFO 跨平台获取 MATLAB 内存使用信息
    %   mem = getMemoryInfo() 返回结构体 mem，包含内存使用信息。
    %   Windows 下调用 MATLAB 内置 memory 函数；
    %   Linux/macOS 下通过 /proc 文件系统获取进程内存信息。
    %
    %   输出字段（尽量与 MATLAB memory 函数兼容）：
    %       MemUsedMATLAB   - 当前 MATLAB 进程占用的物理内存（字节）
    %       MemUsedMATLAB_GB - 同上（GB，方便打印）
    %       MemAvailableAllArrays - 系统可用内存（字节，Linux 近似值）
    %
    %   示例：
    %       mem = getMemoryInfo();
    %       fprintf("    total used memory: %.1e GB\n", mem.MemUsedMATLAB_GB);
    
    if ispc
        % Windows: 直接使用 MATLAB 内置函数
        mem = memory;
        mem.MemUsedMATLAB_GB = mem.MemUsedMATLAB / 1e9;
        
    else
        % Linux / macOS: 通过系统命令获取
        mem = struct();
        
        [status_code, status_out] = system('cat /proc/self/status | grep VmRSS');
        if status_code == 0
            tokens = regexp(status_out, '\d+', 'match');
            if ~isempty(tokens)
                mem.MemUsedMATLAB = str2double(tokens{1}) * 1024;  % kB -> bytes
                mem.MemUsedMATLAB_GB = mem.MemUsedMATLAB / 1e9;
            else
                mem.MemUsedMATLAB = NaN;
                mem.MemUsedMATLAB_GB = NaN;
            end
        else
            mem.MemUsedMATLAB = NaN;
            mem.MemUsedMATLAB_GB = NaN;
        end
        
        [status_code, meminfo] = system('cat /proc/meminfo | grep MemAvailable');
        if status_code == 0
            tokens = regexp(meminfo, '\d+', 'match');
            if ~isempty(tokens)
                mem.MemAvailableAllArrays = str2double(tokens{1}) * 1024;  % kB -> bytes
            else
                mem.MemAvailableAllArrays = NaN;
            end
        else
            mem.MemAvailableAllArrays = NaN;
        end

        [status_code, hwm_out] = system('cat /proc/self/status | grep VmHWM');
        if status_code == 0
            tokens = regexp(hwm_out, '\d+', 'match');
            if ~isempty(tokens)
                mem.PeakMemUsedMATLAB = str2double(tokens{1}) * 1024;  % kB -> bytes
                mem.PeakMemUsedMATLAB_GB = mem.PeakMemUsedMATLAB / 1e9;
            end
        end
        
    end
end