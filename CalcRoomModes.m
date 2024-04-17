function [modes_length,modes_width,modes_depth,SchroederFreq] = CalcRoomModes(room_dim, RT60 ,num_modes)
    c = 343; % m/s
    f_length = c/room_dim(1);
    f_width = c/room_dim(2);
    f_depth = c/room_dim(3);
    modes_length = zeros(1,num_modes);
    modes_width = zeros(1,num_modes);
    modes_depth = zeros(1,num_modes);

    for i = 1:num_modes
        modes_length(i) = f_length*(0.5*i); 
        modes_width(i) = f_width*(0.5*i);
        modes_depth(i) = f_depth*(0.5*i);
    end

    SchroederFreq = 2000*sqrt(RT60/prod(room_dim));
end