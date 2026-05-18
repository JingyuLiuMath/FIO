function memGB = getMaxAvailableMemory()

if ispc
    mem = memory;
    memGB = mem.MemAvailableAllArrays / 1024^3;    
elseif isunix
    [status, result] = system('cat /proc/meminfo | grep MemAvailable');
    if status == 0
        tokens = regexp(result, '(\d+)\s*kB', 'tokens');
        if ~isempty(tokens)
            kb = str2double(tokens{1}{1});
            memGB = kb / 1024^2;  % kB -> GB
        else
            memGB = NaN;
        end
    else
        memGB = NaN;
    end
else
    memGB = NaN;
end

memGB = memGB * 0.7;

end