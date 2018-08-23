% Cancer Detection
% Pre-processing: Median and Gabor filter 
% Segmentation: Thresholding
% Feature Extraction: Classification using threshold
% Written by Shu Leong & Moana Yamanouchi

% Read image
I = imread('IM-0001-0055.jpg');

I_resize=imresize(I, [500,500]);
I_resize=im2double(I_resize);
figure, imshow(I_resize);

% Apply median filter 
I_resize = medfilt2(I_resize);

%Initialize the paramters for Gabor filter
gamma=0.3;
psi=0;
theta=0;
bw=2.8;
lambda=3.0;
pi=0;

% Apply Gabor filter
for x =1:500
    for y=1:500
   
        x_theta=I_resize(x,y)*cos(theta)+I_resize(x,y)*sin(theta);
        y_theta=-I_resize(x,y)*sin(theta)+I_resize(x,y)*cos(theta);
        
        gb(x,y)=exp(-(x_theta.^2/2*bw^2+gamma^2*y_theta.^2/2*bw^2))*cos(2*pi/lambda*x_theta+psi);
    end
end
gb = imcomplement(gb);

% Set up variables
roundnessThreshold = 0.25; % tumors generally have some semblance of roundness
areaThreshold = 100; % blood vessels are generally smaller than this size
sizeOfObject = 0; % to record the largest object's size
traced = 0; % tumors found or not found

% Segment lung images using thresholding method (Otsu's)
BW = imbinarize(I);
figure, imshow(BW);

% Clear unnecessary part of the image using imclearborder
objectsOnly = imclearborder(BW);

% Find largest 5 objects (tumors are generally greater in size than blood vessels)
largest5Objects = bwareafilt(objectsOnly, 5);
figure, imshow(largest5Objects); title('largest 5 objects');

% Traces the exterior boundaries of objects
[B,L] = bwboundaries(largest5Objects,'noholes');

% Returns measurements for area and center of mass of region for each labeled region in the label matrix L
stats = regionprops(L,'Area','Centroid');

% Loop over the boundaries
for k = 1:length(B)

  % Obtain (X,Y) boundary coordinates corresponding to label 'k'
  boundary = B{k};

  % Compute a simple estimate of the object's perimeter
  delta_sq = diff(boundary).^2;    
  perimeter = sum(sqrt(sum(delta_sq,2)));
  
  % Obtain the area calculation corresponding to label 'k'
  area = stats(k).Area;
  
  % Compute the roundness metric
  metric = 4*pi*area/perimeter^2;
  metric_string = sprintf('%2.2f',metric);

  % Compare against pre-defined roundness and area threshold
  % Cancer nodules should have some semblance of roundness and greater in
  % size than blood vessels (which normally are around 30 to 60 pixels in
  % size)
  if metric > roundnessThreshold
      if area > sizeOfObject && area > areaThreshold
          sizeOfObject = area;
          traced = B{k};
      end
  end
  
  % For each potential cancer nodules, list out their roundness metrics
  text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','r',...
       'FontSize',14,'FontWeight','bold');
  
end

% If there is a detected cancer nodules, draw the boundary of said nodules
% in the original CT scan image in red
if traced ~= 0
    figure, imshow(I)
    hold on;
    plot(traced(:,2), traced(:,1), 'r', 'LineWidth', 2)
% Else, notify user that the CT scan belongs to healthy patient
else
    disp('CT scan belongs to normal patient');
end