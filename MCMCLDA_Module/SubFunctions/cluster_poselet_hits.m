function cluster_id=cluster_poselet_hits(hits_for_img,hyps_for_img,kl_dists)


global config;
cluster_thresh= config.GREEDY_CLUSTER_THRESH;
cluster_id=saliency_cluster_poselet_hits(hits_for_img,hyps_for_img,kl_dists,cluster_thresh);

end


function dist=torso_bounds_distance(hyps,hits)
    global K;
    bounds1 = torso_bounds_from_keypoints(hyps(1).mu([K.L_Shoulder K.R_Shoulder K.L_Hip K.R_Hip],:));
    bounds2 = torso_bounds_from_keypoints(hyps(2).mu([K.L_Shoulder K.R_Shoulder K.L_Hip K.R_Hip],:));
    dist=1-bounds_match(bounds1,0,bounds2,0);
end

function dist=bounds_distance(hyps,hits)
    global gmodel;
    bounds1 = predict_bounds(hits.bounds(1,:),0,gmodel.hough_votes(hits.poselet_id(1)));
    bounds2 = predict_bounds(hits.bounds(2,:),0,gmodel.hough_votes(hits.poselet_id(2)));
    dist=1-bounds_match(bounds1,0,bounds2,0);
end


function dist=kl_distance(hyps,hits)
    dist=hyps(1).distance(hyps(2));
end

function cluster_id=saliency_cluster_poselet_hits(hits_for_img,hyps_for_img,dists,dist_thresh)
global config;
global K;




if config.DEBUG>0
    DISPLAY_STEPS=10;
end

if hits_for_img.isempty
    cluster_id=[];
    return;
end


hyp_clusters{1,1}.hit_id(1)=1;
hyp_clusters{1,1}.poselet_id(1)=hits_for_img.poselet_id(1);

if config.DEBUG>2
    imshow(config.DEBUG_IMG);
    draw_hyp_onecolor(hyps_for_img(1),hits_for_img.select(1)  ,[0 1 0]);

    keyboard;
end

for j=2:hits_for_img.size
   
    dst=nan(length(hyp_clusters),1);
    for k=1:length(hyp_clusters)
       
            hits_to_check=1:length(hyp_clusters{k}.hit_id);
    end

        if 0
            dst_m=inf(1,length(hits_to_check));
            for m=1:length(hits_to_check)
                hit_m=hyp_clusters{k}.hit_id(hits_to_check(m));
                
                dst_m(m) = hyps_for_img(hyp_clusters{k}.hit_id(hits_to_check(m))).distance(hyps_for_img(j),config);
                
                
                if isinf(dst_m(m))
                    break;
                end
            end
        else
           
            hits_to_check_idx = hyp_clusters{k}.hit_id(hits_to_check);
            dst_m2 = dists(hits_to_check_idx,j);
            dst_m2 = dst_m2(:)';

            dst_m = dst_m2;

        end
      
        
        if any(isinf(dst_m))
            dst(k)=inf;
        else
            scores = hits_for_img.score(hyp_clusters{k}.hit_id(hits_to_check));
            dst(k)=(dst_m*scores)/sum(scores);
        end
        
        0;
    end
    if config.DEBUG>2
        mrg=hyp_clusters{find(dst==min(dst),1)}.hit_id(1);
        imshow(config.DEBUG_IMG);
        hold on;

        draw_hyp_onecolor(hyps_for_img(mrg),hits_for_img.select(mrg),[1 0 0]);
        draw_hyp_onecolor(hyps_for_img(j),hits_for_img.select(j)  ,[0 1 0]);
        title(sprintf('min distance=%f',min(dst)));

        if ~exist('poselets','var')
            poselets=load_var([config.TMP_DIR '/poselet_lib.mat'],'poselets');
            a=load_var([config.POSELET_DIR '/init.mat'],'a');
        end
          thumbs = combine_thumbs_of_poselets(extract_thumbs_of_poselet(a,set_poselets_max_examples(poselets(hits_for_img.poselet_id([j mrg])),5)));
        figure(2);
        display_thumbs(thumbs,[],[2 5]);
        title('Red (first row) and green (second row)');
        figure(1);
        keyboard;
    end
    if min(dst)<dist_thresh
        mrg=find(dst==min(dst),1);
        poselet_id = hits_for_img.poselet_id(j);
        
        if ~ismember(poselet_id,hyp_clusters{mrg}.poselet_id)
            hyp_clusters{mrg}.poselet_id(end+1)=poselet_id;
            hyp_clusters{mrg}.hit_id(end+1)=j;
        end
    elseif length(hyp_clusters)<config.HYP_CLUSTER_MAXNUM
       
        hyp_clusters{end+1,1}.hit_id(1)=j;
        hyp_clusters{end,1}.poselet_id(1)=hits_for_img.poselet_id(j);
    end
    if config.DEBUG>0 && mod(j,DISPLAY_STEPS)==0
        disp(sprintf('%d of %d',j,hits_for_img.size));
        imshow(config.DEBUG_IMG);
        hold on;
        for k=1:length(hyp_clusters)
            draw_cluster(hits_for_img.select(hyp_clusters{k}.hit_id),hyps_for_img(hyp_clusters{k}.hit_id),[K.L_Shoulder K.R_Shoulder]);
        end
        keyboard;
    end
end


cluster_id = zeros(hits_for_img.size,1);
for k=1:length(hyp_clusters)
    cluster_id(hyp_clusters{k}.hit_id) = k;
end

end


function draw_hyp_onecolor(merge_hyps,hit,color)
global K;
keypts_range=[K.L_Hip K.R_Hip K.L_Shoulder K.R_Shoulder];
merge_hyps.draw(keypts_range,repmat(color,length(keypts_range),1),'-',1);
rectangle('position',[merge_hyps.rect(1:2) merge_hyps.rect(3:4)-merge_hyps.rect(1:2)],'edgecolor',color,'linestyle',':');
hit.draw_bounds(color);
end

function draw_cluster(hits,hyps, keypts)
global K;
MAX_HITS = 100;
if hits.size>MAX_HITS
    hits = hits.select(1:MAX_HITS);
    hyps = hyps(1:MAX_HITS);
end

kp_mu=reshape([hyps(:).mu],K.NumPrimaryKeypoints,2,[]);

colors = jet(length(keypts));
for k=1:length(keypts)
    num_samples=min(10,size(kp_mu,3));
    scatter(kp_mu(keypts(k),1,1:num_samples),kp_mu(keypts(k),2,1:num_samples),'.','MarkerEdgeColor',colors(k,:));
end


end






function [bounds,angle] = predict_bounds(poselet_bounds, poselet_angle, poselet2bounds)
    

    scale = min(poselet_bounds(3:4)); 
    image2poselet_ctr = poselet_bounds(1:2)+poselet_bounds(3:4)/2;
    rad_angle = poselet_angle*pi/180;
    poselet_rot = [cos(rad_angle) sin(rad_angle); -sin(rad_angle) cos(rad_angle)];

    scaled_bounds = poselet2bounds.obj_bounds * scale;
    poselet2bounds_ctr = scaled_bounds(1:2) + scaled_bounds(3:4)/2;
    bounds_dims = scaled_bounds(3:4);

    image2bounds_ctr = image2poselet_ctr + poselet2bounds_ctr*poselet_rot;        
    bounds = [image2bounds_ctr - bounds_dims/2 bounds_dims];
    angle = poselet_angle;
end


function [torso_bounds,torso_angle,torso_score]=compute_torso_bounds(hits_for_torso,hyps_for_torso)
    global K;
    
    torso_kpts = [K.L_Shoulder K.R_Shoulder K.L_Hip K.R_Hip];
    
    kp_mu=reshape([hyps_for_torso(:).mu],size(hyps_for_torso(1).mu,1),2,[]);     
    torso_score=sum(hits_for_torso.score);

    if 1
        for kp=1:length(torso_kpts)
            coords = reshape([kp_mu(torso_kpts(kp),1,:) kp_mu(torso_kpts(kp),2,:)],2,[]);
            mean_coords(kp,:) = sum([hits_for_torso.score hits_for_torso.score].*coords',1)/torso_score;
        end    
    else
        kp_var=reshape([hyps_for_torso(:).sigma],size(hyps_for_torso(1).sigma,1),2,[]);     
        for kp=1:length(torso_kpts)
              coords = shiftdim([kp_mu(torso_kpts(kp),1,:) kp_mu(torso_kpts(kp),2,:)],1)';
              var    = shiftdim([kp_var(torso_kpts(kp),1,:) kp_var(torso_kpts(kp),2,:)],1)';
              mean_coords(kp,:) =get_mode(coords,var,hits_for_torso.score);
        end
    end

    [torso_bounds,torso_angle] = torso_bounds_from_keypoints(mean_coords);
end


function md = get_mode(x,sigma,w)

   w = w./sigma;
   md = (x'*w)/sum(w);
   return;


  [modes,mode_w] = meanshift(x, sigma.^2, w);
  md=modes(find(mode_w==max(mode_w),1),:);
end



function [modes,mode_w,mode_of_x] = meanshift(at, sigma2, w)

x = at;
[N,D] = size(x);

if N==1
   modes = at;
   mode_w = w;
   mode_of_x = 1;
   return;
end

THRESH = 1e-5;
MODE_EPS = 20;
MAX_ITERS = 100;

for iter=1:MAX_ITERS
    for i=1:N
       xdiff = (x - repmat(x(i,:),N,1))./sigma2;
       wt = w.*exp(-sum(xdiff.^2,2) / 2)./sqrt(prod(sigma2,2));
       
       sumw = sum(wt);
       if sumw>0       
           nvalue(i,:) = sum(repmat(wt,1,D).*x,1)./sumw;
       else
           nvalue(i,:) = x(i,:);
       end
    end
    
    shift_sqrd_dist = sum((x - nvalue).^2,2);

    if mean(shift_sqrd_dist)<THRESH
        break;
    end

    x = nvalue;
end


[md,foo,mode_of_x] = unique(round(x*MODE_EPS),'rows');
modes = md./MODE_EPS;


for i=1:size(modes,1)
    xdiff = (at - repmat(modes(i,:),N,1))./sigma2;
    dist = exp(-sum(xdiff.^2,2) / 2)./sqrt(prod(sigma2,2));
    mode_w(i,1) = w' * dist;
end
end


