function expanded_data = Fn_expand_space(data, K, N, Nz, SS_expan_L, SS_expan_Lz)
% expand the data size for the patch moving
% wrote by Hank 07/05/2022
DimNum = size(size(data),2);
%SS_expan_L  = N+1;
expan_L = (N+N)*2;
expan_Lz = (Nz+Nz)*2;

switch DimNum
    case 3
        expanded_data = zeros(size(data, 1)+expan_L, size(data, 2)+expan_L, size(data, 3)+expan_Lz);
        expanded_data(1+SS_expan_L: size(expanded_data,1)-SS_expan_L, ...
            1+SS_expan_L: size(expanded_data,2)-SS_expan_L, ...
            1+SS_expan_Lz: size(expanded_data,3)-SS_expan_Lz) = data;

    case 4
        expanded_data = zeros(K, size(data,2)+expan_L, size(data,3)+expan_L, size(data,4)+expan_Lz);
        expanded_data(:,1+SS_expan_L: size(expanded_data,2)-SS_expan_L, ...
            1+SS_expan_L: size(expanded_data,3)-SS_expan_L, ...
            1+SS_expan_Lz: size(expanded_data,4)-SS_expan_Lz) = data(1:K, :,:,:);
end
end