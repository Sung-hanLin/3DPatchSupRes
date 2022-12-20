% Made by Sung-Han Lin (internal version: V10)
clear; clc; close all;
load('Digi_Phantom_Cardi.mat');
h_HR = h_mri; %structural images
segment_map = p_maps; %map of segmented structural images
c_input = c_mri_bic_DS_SS;   %single slice C13 image (low-res)

%% Setting up parameters
TargSize = 192;  % target reconstruction matrix size in X-Y plane
TargDim = [TargSize TargSize size(h_HR, 3)];
K = 3;  %number of segmetation tissue
N = 2;  %patch size (actual patch  = (N*2)+1)
L_Vi = N*2+1; %length of search volume (Vi)
SS_expan_L  = N+N; % length of expandsion in X, Y directions
max_patch_z = floor((size(h_HR,3)-1)/2); %maximum patch size in Z direction
if N>=max_patch_z
    Nz = max_patch_z;
else 
    Nz = N; 
end
L_Viz = Nz*2+1;
SS_expan_Lz  = Nz+Nz; % length of expandsion in Z direction

iteration_target = 10; %initial iteration number for estimating total iteration number
iter_thresh = 0.001; %termination threshold slope of the tangent line of the fitting function
filename = ['Recon_Img_K',num2str(K),'N', num2str(N),'Nz', num2str(Nz),'.mat'];

%% Upsampling low-res 13C / Downdampling high-res 1H
%Upsampling low-res 13C
c_low = c_input;
[up_x,up_y] = meshgrid(1:(length(c_input)-1)/(TargSize-1):length(c_input));
c_init_1 = interp2(c_low,up_x,up_y,'spline');   % bilinear interpolation
%Image display
figure; subplot(1,2,1); imagesc(c_low); title('Original 13C'); colorbar; hold on; 
subplot(1,2,2); imagesc(c_init_1); colorbar; title('upsampled 13C'); colormap jet

% Downsampling High-res 1H
down_H = zeros(size(c_init_1,1), size(c_init_1,2), size(h_HR,3));
down_p = zeros(size(segment_map,1), size(c_init_1,1),size(c_init_1,2), size(segment_map,4));
downFactor = size(c_init_1,1) / size(segment_map,2);

for s = 1:size(h_HR,3)
    down_H(:,:,s) = imresize(h_HR(:,:,s), downFactor,'nearest');
    for k = 1:size(segment_map,1)
        down_p(k,:,:,s) = imresize(squeeze(segment_map(k,:,:,s)), downFactor,'nearest');
    end
end

%% Manual image alignment between 1H and 13C images
h_low_sum = sum(down_H,3);
shift_y = 0; %int, shifting pixels #  in y axis
shift_x = 0; %int, shifting pixels #  in x axis
p_shift = zeros(size(down_p)); down_H_shift = zeros(size(down_H));

for k = 1:size(down_p,1)
    for s = 1:size(down_H,3)
        down_H_shift(:,:,s) = immove(down_H(:,:,s),shift_y,shift_x);
        p_shift(k,:,:,s) = immove(squeeze(down_p(k,:,:,s)),shift_y,shift_x);
    end
end
down_p = p_shift;
down_H = down_H_shift;

% Image display (Overlay aligned 13C functional and structural images) 
figure; 
H1_base = down_H_shift(:,:, round(size(down_H_shift,3)/2,0));
imshow(H1_base, []); title('After overlay C13 and H1');hold on;
color = cat(3, ones(size(H1_base))*0.929, ones(size(H1_base))*0.694, ones(size(H1_base))*0.125);
h_HR = imshow(color); set(h_HR, 'AlphaData', c_init_1*2e-4); colorbar;

%% Create initial guess of multiple-slice 13C images
c_init_multi = c_init_1 / size(down_H,3); %divided the signal intensity to #slices
c_init =repmat(c_init_multi, 1,1,size(down_H,3));

% Image display
disrange_mx = max(c_init, [] , 'all'); disrange_mi = min (c_init,[], 'all');
figure; montage(c_init, 'size', [1,size(c_init,3)], 'DisplayRange', ...
    [disrange_mi, disrange_mx]); colorbar; colormap jet; title('initial multiple-slice 13C');

%% Preparation for Iterative super-res algorithm (expand space)
p_1 = expand_space(down_p, K, N, Nz, SS_expan_L, SS_expan_Lz); %Probability maps (segmented structural images)
p_1 = single(p_1);
h_low_1 = expand_space(down_H, K, N, Nz, SS_expan_L, SS_expan_Lz); %Down-sampled structural images
h_low_1 = single(h_low_1);
input_1 = expand_space(c_init, K, N, Nz, SS_expan_L, SS_expan_Lz); %Initial guess of 13C images
input_1 = single(input_1);

%% initial step of reconstruction: calculation the weight matrix & Mean Correction
input_2_1H = zeros(size(input_1), 'single');
xxi = 1+(N+N); xxj=size(input_1, 1)-(N+N);
yyi = 1+(N+N); yyj=size(input_1, 2)-(N+N);
zzi = 1+(Nz+Nz); zzj=size(input_1, 3)-(Nz+Nz);

poolobj = gcp('nocreate'); delete(poolobj); parpool(6);
w = zeros(TargDim(1)+SS_expan_L, TargDim(2)+SS_expan_L,...
    TargDim(3)+SS_expan_Lz, L_Vi, L_Vi, (Nz*2+1), 'single');

tic
parfor s = zzi:zzj
    for j = yyi:yyj
        for i = xxi:xxj
            w(i,j,s,:,:,:) = Cal_weight_H_3D_vector_v2([i,j,s],K,p_1,N,Nz,h_low_1, L_Vi, L_Viz);        
        end
    end
end
t = toc;
disp(['weighting patch cal t =',num2str(t) , 'sec']);
weighting_patch_cal_time = t; 
[input_2_1H_re, t] = SR_recon_s2(w,input_1, input_2_1H, N, Nz, xxi,xxj,yyi,yyj,zzi,zzj);
input_2 = input_2_1H_re*1;

%Mean Correction
iteration = 1;
error = zeros(1,iteration_target, 'single');
[input_2,error(1)] = MeanCorrection(input_2, input_1,c_low, c_input, TargDim,TargSize, SS_expan_L, SS_expan_Lz);
fprintf('#iter. = %d; Recon. t. =  %d sec; error = %d\n', iteration, t, error(iteration));

while iteration < iteration_target
    iteration = iteration+1;
    input_1 = input_2;
    input_2_1H = zeros(size(input_1), 'single');
    [input_2_1H, t] = SR_recon_s2(w,input_1, input_2_1H, N, Nz, xxi,xxj,yyi,yyj,zzi,zzj);

    input_2 = input_2_1H*1;
    [input_2 ,error(iteration)] = MeanCorrection(input_2, input_1,c_low, c_input, TargDim,TargSize, SS_expan_L, SS_expan_Lz);
    fprintf('#iter. = %d; Recon. t. =  %d sec; error = %d\n', iteration, t, error(iteration));
end

%% Calculation the Termination number of iteration
[xData, yData] = prepareCurveData([], error);
ft = fittype( 'power2' ); opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); opts.Display = 'Off';
[fitresult, gof] = fit( xData, yData, ft, opts);
coeffs = coeffvalues(fitresult); 
syms x
a = coeffs(1); b = coeffs(2); c = coeffs(3);
f = a*(x.^b)+c; df = diff(f);

predi_iteration = 1;
while abs(eval(subs(df, x, predi_iteration)))> iter_thresh
    slope = eval(subs(df, x, predi_iteration));
    predi_iteration = predi_iteration+1;

    if predi_iteration>1000
        predi_iteration = 1000;
        disp(['Predicted_iteration could not reach convergence, iter# set to', num2str(predi_iteration)]);
        break
    end
end
iteration_target = predi_iteration;
disp(['Predicted_iteration # = ', num2str(predi_iteration)]);

%Plotting errors of mean correction
figure; 
while iteration < iteration_target
    iteration = iteration+1;
    input_1 = input_2;
    input_2_1H = zeros(size(input_1), 'single');
    [input_2_1H, t] = SR_recon_s2(w,input_1, input_2_1H, N, Nz, xxi,xxj,yyi,yyj,zzi,zzj);

    input_2 = input_2_1H*1;
    [input_2 ,error(iteration)] = MeanCorrection(input_2, input_1,c_low, c_input, TargDim,TargSize, SS_expan_L, SS_expan_Lz);
     fprintf('#iter. = %d; Recon. t. =  %d sec; error = %d\n', iteration, t, error(iteration));
     plot(1:iteration,error); xlabel('Number of Iterations');ylabel('max(x(i+1)-x(i))'); hold on
     drawnow
end
hold off
e_final = max(abs(input_2 - input_1),[],'all');
disp(['# iteration_final = ',num2str(iteration)]);

% Recon results
c_super_res =  input_2(1+SS_expan_L:TargDim(1)+SS_expan_L,...
    1+SS_expan_L:TargDim(2)+SS_expan_L,...
    1+SS_expan_Lz:TargDim(3)+SS_expan_Lz);

%% Display and Save results 
save( filename,'c_super_res','error', 'c_init', 'iteration', 'K', 'N', 'Nz', 'input_1', 'input_2', ...
    'weighting_patch_cal_time', 'shift_x', 'shift_y');
poolobj = gcp('nocreate'); delete(poolobj);

DisRag_mx_SS = max(c_low, [], 'all'); DisRag_mi_SS = min(c_low, [], 'all');
DisRag_mx_3D = max(c_init, [], 'all'); DisRag_mi_3D = min(c_init, [], 'all');
DisRag_mx_H = max(down_H, [], 'all'); DisRag_mi_H = min(down_H, [], 'all');

%Display Input and Recon 13C
figure; subplot(1,2,1); imshow(c_init_1,[DisRag_mi_SS DisRag_mx_SS]); title('Input 13C'); colorbar; colormap jet; 
hold on;
subplot(1,2,2); imshow(sum(c_super_res,3),[DisRag_mi_SS DisRag_mx_SS]); title('Recon 13C'); colorbar; colormap jet; 
hold off;

%Display Recon 13C and 1H
figure; subplot(2,1,1); montage(c_super_res, 'size', [1,size(c_super_res,3)], 'DisplayRange', ...
    [DisRag_mi_3D, DisRag_mx_3D]); colorbar; hold on;
subp2 = subplot(2,1,2); montage(down_H, 'size', [1,size(c_super_res,3)], 'DisplayRange', ...
    [DisRag_mi_H, DisRag_mx_H]); colorbar; colormap jet; hold off; 