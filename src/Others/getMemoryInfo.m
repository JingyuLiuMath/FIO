function memGB = getMemoryInfo()
if ispc
    % Windows
    mem = memory;
    memGB = mem.MemUsedMATLAB / 1024^3;
else
    [status, result] = system('cat /proc/self/status | grep VmRSS');
    if status == 0
        tokens = regexp(result, '(\d+)\s*kB', 'tokens');
        if ~isempty(tokens)
            kb = str2double(tokens{1}{1});
            memGB = kb / 1024^2;
        else
            memGB = NaN;
        end
    else
        memGB = NaN;
    end
end

end