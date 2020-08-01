clc

FDetect = vision.CascadeObjectDetector;

video = videoinput('winvideo', 1);
set(video, 'ReturnedColorSpace', 'RGB');
capturedImage = getsnapshot(video);

Image = capturedImage;

BB = step(FDetect, Image);

figure, imshow(Image); hold on
for i = 1: size(BB, 1)
   rectangle('Position', BB(i, :), 'LineWidth', 2, 'LineStyle', '-');
end    

title('Face Recognition');
hold off;