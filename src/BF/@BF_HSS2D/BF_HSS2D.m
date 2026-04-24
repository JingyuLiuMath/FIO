classdef BF_HSS2D < HSS_Herm
    
    properties
        nx_ (1, 1) double;
        ny_ (1, 1) double;
        x_freq_start_ (1, 1) double;
        x_freq_end_ (1, 1) double;
        x_size_ (1, 1) double;
        y_freq_start_ (1, 1) double;
        y_freq_end_ (1, 1) double;
        y_size_ (1, 1) double;

        perm_ (:, 1) double;
        perm_inv_ (:, 1) double;
        
        ind_ (:, 1) double;
        re_ (:, 1) double;
        sk_ (:, 1) double;
    end
    
     methods
         function A = BF_HSS2D(...
                 nx, ny, ...
                 x_freq_start, x_freq_end, ...
                 y_freq_start, y_freq_end, ...
                 level, offset)
             arguments (Input)
                 nx (1, 1) double;
                 ny (1, 1) double;
                 x_freq_start (1, 1) double = -nx / 2;
                 x_freq_end (1, 1) double = nx / 2 - 1;
                 y_freq_start (1, 1) double = -nx / 2;
                 y_freq_end (1, 1) double = nx / 2 - 1;
                 level (1, 1) double = 0;
                 offset (1, 1) double = 0;
             end

             arguments (Output)
                 A BF_HSS2D;
             end

             N = nx * ny;
             A.nx_ = nx;
             A.ny_ = ny;
             A.global_size_ = N;
             A.x_freq_start_ = x_freq_start;
             A.x_freq_end_ = x_freq_end;
             A.x_size_ = x_freq_end - x_freq_start + 1;
             A.y_freq_start_ = y_freq_start;
             A.y_freq_end_ = y_freq_end;
             A.y_size_ = y_freq_end - y_freq_start + 1;
             A.size_ = A.x_size_ * A.y_size_;
             A.level_ = level;
             A.offset_ = offset;

         end
         
     end

end