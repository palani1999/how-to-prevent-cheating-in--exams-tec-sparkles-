function [hog,samples_x,samples_y]=compute_hog(img)

global config;

img=single(img);

if ~config.USE_PHOG
    [hog,samples_x,samples_y]=compute_hog_internal(img);
else
    
    
    bandwidth = config.HOG_CELL_DIMS ./ config.NUM_HOG_BINS;
    [H,W,D]=size(img);
    if ~all([W H]>=(bandwidth(1:2).*config.NUM_HOG_BINS(1:2)+bandwidth(1:2)+2)*2*2)
        hog=[];
        samples_x=[];
        samples_y=[];
        return;
    end
    
   
    leftover = mod([W H]-8,bandwidth(1:2).*config.NUM_HOG_BINS(1:2)*2);
    top_left = floor(leftover/2);
    if any(leftover>0)
        bot_right = leftover-top_left;
        img = img((top_left(2)+1):(end-bot_right(2)),(top_left(1)+1):(end-bot_right(1)),:);
        [H,W,D]=size(img);
        assert((mod(H-8,bandwidth(2)*4)==0) && (mod(W-8,bandwidth(1)*4)==0));
    end
       
   
    [hg,samples_x,samples_y]=compute_hog_internal(img(4:(end-3),4:(end-3),:));
    assert(samples_x(1)==1 && samples_y(1)==1);
    
    samples_x=samples_x+top_left(1)+3; % add back the margin
    samples_y=samples_y+top_left(2)+3;
    hog.hog1=hg;
    
    im2=(img(1:2:end,1:2:end,:)+img(2:2:end,1:2:end,:)+img(1:2:end,2:2:end,:)+img(2:2:end,2:2:end,:))/4;        
    band2=bandwidth(1:2)/2;
    for x=0:1
        for y=0:1
           
            [hg,sx,sy] = compute_hog_internal(im2((2+y*band2(2)):(end-1-(2-y)*band2(2)),...
                                                  (2+x*band2(1)):(end-1-(2-x)*band2(1)),:));
            assert(sx(1)==1 && sy(1)==1);
            hog.hog2(y+1,x+1,:,:,:) =hg;
        end
    end

    im4=(im2(1:2:end,1:2:end,:)+im2(2:2:end,1:2:end,:)+im2(1:2:end,2:2:end,:)+im2(2:2:end,2:2:end,:))/4;
    band4=bandwidth(1:2)/4;
    for x=0:3
        for y=0:3
            [hg,sx,sy] = compute_hog_internal(im4((1+y*band4(2)):(end-(4-y)*band4(2)),...
                                                  (1+x*band4(1)):(end-(4-x)*band4(1)),:));
            assert(sx(1)==1 && sy(1)==1);
            hog.hog4(y+1,x+1,:,:,:) =hg;                    
        end
    end
end

end

