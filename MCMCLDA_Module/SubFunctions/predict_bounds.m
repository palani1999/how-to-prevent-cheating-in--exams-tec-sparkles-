

function bounds = predict_bounds(poselet_bounds, poselet2bounds)
   

    scale = min(poselet_bounds(3:4));
    image2poselet_ctr = poselet_bounds(1:2)+poselet_bounds(3:4)/2;

    scaled_bounds = poselet2bounds.obj_bounds * scale;
    poselet2bounds_ctr = scaled_bounds(1:2) + scaled_bounds(3:4)/2;
    bounds_dims = scaled_bounds(3:4);

    image2bounds_ctr = image2poselet_ctr + poselet2bounds_ctr;        
    bounds = [image2bounds_ctr - bounds_dims/2 bounds_dims];
end