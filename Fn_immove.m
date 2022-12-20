function outImage=Fn_immove(inImage,mx,my)
%%
% mx:   int, # of pixels moved along up-down direction (positive --> down)
%                                                 (negtive --> up)
% my:   int, # of pixels moved along left-right direction (positive --> right)
%                                                    (negtive --> left)
%%
[m,n] = size(inImage);
outImage = zeros(size(inImage));

if mx >=0
    if my >= 0
        outImage(1+mx:m,1+my:n) = inImage(1:m-mx,1:n-my);
        outImage(1:mx,:) = inImage(1:mx,:);
        outImage(:,1:my) = inImage(:,1:my);
    else
        outImage(1+mx:m,1:n+my) = inImage(1:m-mx,1-my:n);
        outImage(1:mx,:) = inImage(1:mx,:);
        outImage(:,1+n+my:n) = inImage(:,1+n+my:n);
    end
else
    if my >= 0
        outImage(1:m+mx,1+my:n) = inImage(1-mx:m,1:n-my);
        outImage(1+m+mx:m,:) = inImage(1+m+mx:m,:);
        outImage(:,1:my) = inImage(:,1:my);
    else
        outImage(1:m+mx,1:n+my) = inImage(1-mx:m,1-my:n);
        outImage(1+m+mx:m,:) = inImage(1+m+mx:m,:);
        outImage(:,1+n+my:n) = inImage(:,1+n+my:n);
    end
end
end

