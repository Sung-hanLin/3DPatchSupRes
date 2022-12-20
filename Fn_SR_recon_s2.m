function [input_2_1H, t] = Fn_SR_recon_s2(w,input_1, input_2_1H, N, Nz, xxi,xxj,yyi,yyj,zzi,zzj)
% Super resolution reconstruction, writen by Hank on 07/13/2022

tic
for s = zzi:zzj
    for j = yyi:yyj
        for i = xxi:xxj
            input_2_1H(i,j,s) = sum(squeeze(w(i,j,s,:,:,:)).*input_1(i-N:i+N, j-N:j+N, s-Nz:s+Nz),'all');

        end
    end
end
t = toc;
end