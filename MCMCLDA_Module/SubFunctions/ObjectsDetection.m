function [bounds_predictions,poselet_hits,torso_predictions]=ObjectsDetection(img, model)
global config;



fprintf('Computing pyramid HOG ... ');
total_start_time=clock;
phog=image2phog(img);
fprintf('Done in %4.2f secs.\n',etime(clock,total_start_time));

fprintf('Detecting poselets... ');
start_time=clock;
poselet_hits = nonmax_suppress_hits(detect_poselets(phog,model.svms));
poselet_hits.score = -poselet_hits.score.*model.logit_coef(poselet_hits.poselet_id,1)-model.logit_coef(poselet_hits.poselet_id,2);
poselet_hits.score = 1./(1+exp(-poselet_hits.score));
[srt,srtd]=sort(poselet_hits.score,'descend');
poselet_hits = poselet_hits.select(srtd);
fprintf('Done in %4.2f secs.\n',etime(clock,start_time));
           
if isfield(model,'bigq_weights')
   fprintf('Big Q...');
    start_time=clock;
   hyps=[model.hough_votes.hyp];
   [features,~,kl_dists] = get_context_features_in_image(hyps,poselet_hits);
   poselet_hits.score = sum(features.*model.bigq_weights(poselet_hits.poselet_id,1:(end-1)),2) + model.bigq_weights(poselet_hits.poselet_id,end);
   poselet_hits.score = -poselet_hits.score.*model.bigq_logit_coef(poselet_hits.poselet_id,1)-model.bigq_logit_coef(poselet_hits.poselet_id,2);
   poselet_hits.score = 1./(1+exp(-poselet_hits.score));
   fprintf('Done in %4.2f secs.\n',etime(clock,start_time));
end

fprintf('Clustering poselets... ');
start_time=clock;
hyps = instantiate_hypotheses([model.hough_votes.hyp],poselet_hits);
cluster_labels = cluster_poselet_hits(poselet_hits,hyps,kl_dists);

fprintf('Done in %4.2f secs.\n',etime(clock,start_time));

fprintf('Predicting bounding boxes... ');
start_time=clock;
[H W D] = size(img);
[torso_predictions,bounds_predictions] = cluster2bounds(poselet_hits, hyps, cluster_labels, model, [W H]);
fprintf('Done in %4.2f secs.\n',etime(clock,start_time));

disp(sprintf('Total time: %f secs',(etime(clock,total_start_time))));

