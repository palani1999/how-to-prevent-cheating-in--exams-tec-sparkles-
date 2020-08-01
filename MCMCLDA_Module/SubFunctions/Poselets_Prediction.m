function [bounds_predictions,poselet_hits,torso_predictions] = Poselets_Prediction(img)
global config;

time=clock;
init;

config.USE_MEX_HOG = true;

faster_detection = true;  
interactive_visualization = false;
enable_bigq = true; 

if faster_detection
    config.DETECTION_IMG_MIN_NUM_PIX = 240^2; 
    config.DETECTION_IMG_MAX_NUM_PIX = 640^2;  
    config.PYRAMID_SCALE_RATIO = 2;
end


load('person-model.mat'); 
if ~enable_bigq
   model =rmfield(model,'bigq_weights');
   model =rmfield(model,'bigq_logit_coef');
end

[bounds_predictions,poselet_hits,torso_predictions] = ObjectsDetection(img,model);

warning('off','MATLAB:structOnObject')
bounds_predictions = struct(bounds_predictions);
torso_predictions = struct(torso_predictions);
poselet_hits = struct(poselet_hits);
warning('on','MATLAB:structOnObject')
0;