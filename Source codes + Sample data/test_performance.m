% Test performance for cancer detection
% Written by Shu Leong & Moana Yamanouchi

%Calculate TAR, TRR, FAR, FRR
%setting up directory for cancer CT --> measure TAR & FRR
cancer_dir = '/Users/shuling/Desktop/Cancer CT';
%build full file name from parts
%returns a character vector containing the full path to the file
%filePattern = '/Users/shuling/Desktop/CT/*.jpg'
filePattern = fullfile(cancer_dir, '*.jpg');
%folder content
jpegFiles = dir(filePattern);
cancer_nFiles = length(jpegFiles);
TAR = 0;
TRR = 0;
FAR = 0;
FRR = 0;

for i=1:cancer_nFiles
   baseFilename = jpegFiles(i).name;
   fullFilename = fullfile(cancer_dir, baseFilename);
   I = imread(fullFilename);
   % Set up variables
    roundnessThreshold = 0.25; % tumors generally have some semblance of roundness
    areaThreshold = 100; % blood vessels are generally smaller than this size
    sizeOfObject = 0; % to record the largest object's size
    traced = 0; % tumors found or not found

    % Segment lung images using thresholding method (Otsu's)
    BW = imbinarize(I);

    % Clear unnecessary part of the image using imclearborder
    objectsOnly = imclearborder(BW);

    % Find largest 5 objects (tumors are generally greater in size than blood vessels)
    largest5Objects = bwareafilt(objectsOnly, 5);

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
    end

    % Calculate TAR and FRR
    if traced ~= 0
        TAR = TAR + 1;
    else
        FRR = FRR + 1;
    end
end

%setting up directory for normal CT --> measure TRR & FAR
normal_dir = '/Users/shuling/Desktop/Normal CT';
%build full file name from parts
%returns a character vector containing the full path to the file
%filePattern = '/Users/shuling/Desktop/CT/*.jpg'
filePattern = fullfile(normal_dir, '*.jpg');
%folder content
jpegFiles = dir(filePattern);
normal_nFiles = length(jpegFiles);

for i=1:normal_nFiles
   baseFilename = jpegFiles(i).name;
   fullFilename = fullfile(normal_dir, baseFilename);
   I = imread(fullFilename);
   % Set up variables
    roundnessThreshold = 0.25; % tumors generally have some semblance of roundness
    areaThreshold = 100; % blood vessels are generally smaller than this size
    sizeOfObject = 0; % to record the largest object's size
    traced = 0; % tumors found or not found

    % Segment lung images using thresholding method (Otsu's)
    BW = imbinarize(I);

    % Clear unnecessary part of the image using imclearborder
    objectsOnly = imclearborder(BW);

    % Find largest 5 objects (tumors are generally greater in size than blood vessels)
    largest5Objects = bwareafilt(objectsOnly, 5);

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
    end

    % Calculate TRR & FAR
    if traced ~= 0
        FAR = FAR + 1;
    else
        TRR = TRR + 1;
    end
end
