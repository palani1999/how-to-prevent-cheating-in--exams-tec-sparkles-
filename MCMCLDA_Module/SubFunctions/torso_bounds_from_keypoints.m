function [torso_bounds,torso_angle]=torso_bounds_from_keypoints(lrshoulder_lrhip_coords)



    global config;
    mShoulder = shiftdim(mean(lrshoulder_lrhip_coords(1:2,1:2,:)))';
    mHip      = shiftdim(mean(lrshoulder_lrhip_coords(3:4,1:2,:)))';
    torsoCtr = (mShoulder+mHip)/2;

    spine = mShoulder - mHip;
    torso_length = sqrt(sum(spine.^2,2));
    torso_dir = spine./[torso_length torso_length];

    torso_angle = atan2(torso_dir(:,2),torso_dir(:,1))+pi/2;
    torso_angle(torso_angle>pi) = torso_angle(torso_angle>pi)-2*pi;
    torso_dims = [torso_length/config.TORSO_ASPECT_RATIO torso_length];   
    torso_bounds = [torsoCtr-torso_dims/2 torso_dims]';
end