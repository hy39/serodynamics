function [ init_prev ] = initial_prev_Wu_expon( )
%INITIAL_PREV_XU_UNIFORM Summary of this function goes here
%   Detailed explanation goes here
        init_prev = []; 
        
        for a=1:4
        if a<4
            Ys = [0.0379	0.0178	0.0083	0.0039	0.0019	0.0009	0.0004	0.0002	0.0001];
            seroconvert_obs = sum(Ys);
            naive = 1 - seroconvert_obs;
            init_prev(a,:) = [naive Ys];
        end
        if a==4
            Ys = [0.0238	0.0160	0.0108	0.0073	0.0050	0.0034	0.0023	0.0016	0.0011];
            seroconvert_obs = sum(Ys);
            naive = 1 - seroconvert_obs;
            init_prev(a,:) = [naive Ys];
        end
        end
end

