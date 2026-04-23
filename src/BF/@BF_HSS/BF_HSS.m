classdef BF_HSS < HSS_Herm
    
    properties
        % *****************************************************************
        % PROPERTY: Col Points.
        freq_start_ (1, 1) double;
        freq_end_ (1, 1) double;
        % -----------------------------------------------------------------
    end
    
     methods
         function A = BF_HSS(...
                 N, freq_start, freq_end, ...
                 level, offset)
             arguments (Input)
                 N (1, 1) double;
                 freq_start (1, 1) double = -N / 2;
                 freq_end (1, 1) double = N / 2 - 1;
                 level (1, 1) double = 0;
                 offset (1, 1) double = 0;
             end

             arguments (Output)
                 A BF_HSS;
             end

             A.global_size_ = N;
             A.freq_start_ = freq_start;
             A.freq_end_ = freq_end;
             A.size_ = freq_end - freq_start + 1;
             A.level_ = level;
             A.offset_ = offset;

         end
         
     end

end