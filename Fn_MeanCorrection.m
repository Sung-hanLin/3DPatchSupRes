function [input_2,e] = Fn_MeanCorrection(input_2, input_1, c_low, c_input, TargDim, TargSize, SS_expan_L, SS_expan_Lz)
%MeanCorrection, created by Hank on 07/18/2022

input_2_temp = input_2(1+SS_expan_L:TargDim(1)+SS_expan_L,...
    1+SS_expan_L:TargDim(2)+SS_expan_L,...
    1+SS_expan_Lz:TargDim(3)+SS_expan_Lz);
input_2_sum = sum(input_2_temp,3);

%Downsampling
factor = size(c_low,1) / size(input_2_sum,1);
down_input_2_sum = imresize(input_2_sum, factor, 'bilinear');

error = c_low - down_input_2_sum; %errors in low resolution
error_up = imresize(error, 1/factor, 'bilinear');
error_up_3D =repmat(error_up/size(input_2_temp,3), 1,1,size(input_2_temp,3));

input_2_temp_Err = input_2_temp+error_up_3D;
input_2(1+SS_expan_L:TargDim(1)+SS_expan_L,...
    1+SS_expan_L:TargDim(2)+SS_expan_L,...
    1+SS_expan_Lz:TargDim(3)+SS_expan_Lz) = input_2_temp_Err;

e = max(abs(input_2 - input_1), [], 'all');

end
