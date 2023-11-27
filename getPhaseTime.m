function [t_array] = getPhaseTime(t, y, id)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k = 1;
for i = 1:length(y)
    if y(i) == id
        prev_i = i;
        exit = 0;
        while exit == 0
            if i == length(y)
                post_id = i;
                break;
            end
            i = i + 1;
            if y(i) ~= id
                post_i = i;
                exit = 1;
            end
        end
        t_array(k) = t(post_i-1) - t(prev_i);
        k = k + 1;
    end
end

end

