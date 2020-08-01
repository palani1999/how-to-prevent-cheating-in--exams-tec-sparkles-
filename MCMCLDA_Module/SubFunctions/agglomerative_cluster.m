function [x,w,src_idx]=agglomerative_cluster(x,w,thresh,dist_fn,merge_fn,min_clusters)


global config;

if ~exist('min_clusters','var')
   min_clusters=1;
end

if length(w)>config.MAX_AGGLOMERATIVE_CLUSTER_ELEMS
   [srt,srtd]=sort(w);
   weak_samples = srtd(1:(length(w)-config.MAX_AGGLOMERATIVE_CLUSTER_ELEMS));
   orig_idx=setdiff(1:length(w),weak_samples);
   w(weak_samples) = [];
   x(weak_samples,:) = [];
else
   orig_idx=1:length(w);
end


src_idx = cell(size(x,1),1);
for i=1:length(src_idx)
   src_idx{i}=uint32(orig_idx(i));
end


distance_matrix = ones(size(x,1),size(x,1));
for i=1:size(x,1)
    d = dist_fn(x(i,:), x(i+1:size(x, 1),:));
    distance_matrix(i,:) = [ones(i,1)', d' ];
end

while 1
    
    [best_d,min_dist_ind] = min(distance_matrix(:)); 
    [best_i, best_j] = ind2sub(size(distance_matrix),min_dist_ind);
    best_pair = [best_i, best_j];

    
    if best_d>thresh
       break;
    end
    
   
    if size(x, 1) <= min_clusters 
        break;
    end
    
    
    i=best_pair(1);
    j=best_pair(2);

    
    [x_ij,w_ij] = merge_fn(x(i,:),w(i), x(j,:),w(j),  src_idx{i},src_idx{j});
    x(i,:) = x_ij;
    w(i) = w_ij;
    src_idx{i} = [src_idx{i}; src_idx{j}];
    
    
    d = dist_fn(x(i,:), x(i+1:size(x,1),:));
    distance_matrix(i,:) = [ones(i,1)' , d'];
    d = dist_fn(x(i,:), x(1:i-1,:));
    distance_matrix(:, i) = [d ; ones(size(x,1) - i + 1, 1)];
    
    
    m = 1 : size(x, 1);
    m_bin = (m ~= best_j);
    distance_matrix = distance_matrix(m_bin, m_bin);
    
   
    x(j,:) = [];
    w(j) = [];
    src_idx(j) = [];
 end

