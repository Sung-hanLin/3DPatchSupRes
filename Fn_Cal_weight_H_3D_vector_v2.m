function [w_p] = Fn_Cal_weight_H_3D_vector_v2(loc,K,p,N,Nz,data, L_Vi, L_Viz)
% loc: [x,y] location of the center point, loc(1)=i, loc(2)=j
% K: number of tissue type
% p: probability maps;
% N: patch size (size of serch volume)
% data: under-sampled proton image
% Patch center = loc(1), loc(2), loc(3), dimension=NxNxN
% L_Vi: length of search volume (Vi)
%% Parameters &%% Search volume (Vi)
w_p = zeros(K, L_Vi, L_Vi, L_Viz, 'single');
X = loc(1); Y = loc(2); Z = loc(3);
patch_size = L_Vi^2*(L_Viz);
Cen_Ncub_idx = floor(patch_size /2)+1;

%% Patch Center (Ni)
Ni = reshape(data(X-N:X+N, Y-N:Y+N, Z-Nz:Z+Nz), [patch_size,1]); %stretch the patch into a vector.
Ni_neig = Ni([1:Cen_Ncub_idx-1, Cen_Ncub_idx+1:patch_size]); % remove the patch center

%% Neighborhood voxels (Nj)
Nj_all = zeros(patch_size, L_Vi, L_Vi, L_Viz, 'single');
Xx = X-N: X+N; Yy = Y-N: Y+N; Zz = Z-Nz: Z+Nz;
for z =1 : L_Viz
    for y = 1 : L_Vi
        for x = 1 : L_Vi
            Nj_all(:,x,y,z) = reshape(data(Xx(x)-N:Xx(x)+N, Yy(y)-N:Yy(y)+N, Zz(z)-Nz:Zz(z)+Nz), [patch_size,1]);
        end
    end
end
Nj_all_neig = Nj_all([1:Cen_Ncub_idx-1, Cen_Ncub_idx+1:patch_size], :,:,:);

%% Calulation the weighting factors
x = 1 : L_Vi;
y = 1 : L_Vi;
z = 1 : L_Viz;
Diff_Nij = squeeze(sum(Ni_neig - Nj_all_neig (:,x,y,z)));
EXP_Diff_Nij = exp (-0.5 * (1/patch_size) * sqrt(Diff_Nij.^2) ./ std(Ni_neig).^2);

for k = 1:K
    w_p(k,:,:,:) = squeeze(p(k,X,Y,Z).* p(k, Xx, Yy, Zz)).*EXP_Diff_Nij / K;
end

w_p = squeeze(sum(w_p,1));
%NoNaN = ~isnan(w_p);
%Idx_NoNaN = find(NoNaN);
Val_NoNaN = w_p(~isnan(w_p));

w_p(N+1, N+1, Nz+1) = 0; %the patch center should not be calculated
w_p(w_p==inf)=1*1e5; %eliminamte INFs
w_p(isnan(w_p))=1*1e-5; %eliminamte NaNs

if sum(Val_NoNaN) ~=0
    %w_p(Idx_NoNaN) = w_p(Idx_NoNaN)./ (sum(Val_NoNaN)*K); %Normalization for the Vi range
    w_p = w_p / sum(w_p,'all');
end
end