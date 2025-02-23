%% Edge decay stuff

filePath = '//data/visitor/ihma587/id03/20241202/RAW_DATA/Al_CR50_2_REX/Al_CR50_2_REX_Int_comp_mono_pink/Al_CR50_2_REX_Int_comp_mono_pink.h5';
tf = '/1.1/measurement/pco_ff/';
mono = '/2.1/measurement/pco_ff/';
pink = '/3.1/measurement/pco_ff/';

mono  = h5read(filePath, mono);
mono = medfilt2(mono(:,:,1)-115);

tf  = h5read(filePath, tf);
tf = medfilt2(tf(:,:,1)-120);

pink  = h5read(filePath, pink);
pink = medfilt2(pink(:,:,1)-160);

figure; 
    subplot(1,3,1) ; imagesc(mono(900:1500,800:1400));axis image; title 'Monochromatic Beam';set(gca, 'ColorScale', 'log');colormap ('gray')
    subplot(1,3,2) ; imagesc(tf(900:1500,800:1400));axis image; title ' Monochromatic Beam with Focusing';set(gca, 'ColorScale', 'log')
    subplot(1,3,3) ; imagesc(pink(750:1350,800:1400));axis image; title ' Pink Beam';set(gca, 'ColorScale', 'log')
%%

% File path
filePath = '//data/visitor/ihma587/id03/20241202/RAW_DATA/Al_CR50_2_REX/Al_CR50_2_REX_Int_comp_mono_pink/Al_CR50_2_REX_Int_comp_mono_pink.h5';

% Dataset paths
tf_path = '/1.1/measurement/pco_ff/';
mono_path = '/2.1/measurement/pco_ff/';
pink_path = '/3.1/measurement/pco_ff/';

% Read and process images
mono = h5read(filePath, mono_path);
mono = medfilt2(mono(:,:,1) - 115);

tf = h5read(filePath, tf_path);
tf = medfilt2(tf(:,:,1) - 125);

pink = h5read(filePath, pink_path);
pink = medfilt2(pink(:,:,1) - 160);

% Get maximum intensity from monochromatic image
mono_max = mean(mono(:));

% Normalize images relative to mono_max
mono_norm = mono / mono_max;  % Monochromatic max is now 1
tf_norm = tf / mono_max;      % Normalize TF image
pink_norm = pink / mono_max;  % Normalize Pink Beam image

% Display images
figure;
subplot(1, 3, 1);
imagesc(mono_norm(900:1500, 800:1400)); axis image;
title('Monochromatic Beam', 'FontSize', 15, 'FontName', 'Helvetica');
set(gca, 'ColorScale', 'log'); colormap('gray');

subplot(1, 3, 2);
imagesc(tf_norm(900:1500, 800:1400)); axis image;
title('Monochromatic Beam with Focusing', 'FontSize', 15, 'FontName', 'Helvetica');
set(gca, 'ColorScale', 'log');

subplot(1, 3, 3);
imagesc(pink_norm(750:1350, 800:1400)); axis image;
title('Pink Beam', 'FontSize', 15, 'FontName', 'Helvetica');
set(gca, 'ColorScale', 'log');

%% Pink Beam (3.1/pco_ff)

%%% Pink Beam (3.1/pco_ff)

% File and dataset paths
filePath = '//data/visitor/ihma587/id03/20241202/RAW_DATA/Al_CR50_2_REX/Al_CR50_2_REX_Int_comp_mono_pink/Al_CR50_2_REX_Int_comp_mono_pink.h5';
pinkPath = '/1.1/measurement/pco_ff/';

% Effective pixel size in microns
pixelSize = 0.19;

% Read and preprocess the Pink dataset
pink = h5read(filePath, pinkPath);
pink = medfilt2(pink(:,:,1) - 100);

% Replace NaN values and set negative values to 0
pink(isnan(pink)) = 0; % Replace NaNs with 0
pink(pink < 0) = 0;    % Set negative values to 0

% Display the Pink Beam image
figure;
imagesc(pink(750:1350,800:1400)); set(gca, 'ColorScale', 'log');
colormap('gray');
colorbar;
axis image;
title('Pink Beam Image (3.1/pco_ff)');
xlabel('X (pixels)');
ylabel('Y (pixels)');
hold on;

% Let user draw the line and store coordinates
disp('Draw a line on the image for analysis.');
%h = imline(gca);
%position = wait(h); % Wait for user to finish drawing
%x1 = position(1, 1); y1 = position(1, 2);
%x2 = position(2, 1); y2 = position(2, 2);
x1 = 200;
y2 = 300;
x2 = 220;
y2 = 320;

% Draw the final line on the image
line([x1, x2], [y1, y2], 'Color', 'red', 'LineWidth', 1.5);

% Extract intensity profile along the line
pinkProfile = improfile(pink, [x1, x2], [y1, y2]);

% Replace NaN values and set negative values to 0 in the profile
pinkProfile(isnan(pinkProfile)) = 0; % Replace NaNs with 0
pinkProfile(pinkProfile < 0) = 0;    % Set negative values to 0

% Convert distance to microns
lineLength = sqrt((x2 - x1)^2 + (y2 - y1)^2); % Line length in pixels
distanceMicrons = linspace(0, lineLength * pixelSize, length(pinkProfile)); % Distance in microns

% Remove NaN or infinite values from both x and pinkProfile
validIndices = ~isnan(pinkProfile) & ~isinf(pinkProfile);
x = distanceMicrons(validIndices);
pinkProfile = pinkProfile(validIndices);

% Normalize the intensity profile
pinkProfile = pinkProfile - min(pinkProfile);
pinkProfile = pinkProfile / max(pinkProfile);

% Fit an error function to the Edge Spread Function (ESF)
erfFit = fit(x', pinkProfile, 'a*erf((x-b)/c) + d', 'StartPoint', [1, mean(x), 1, 0]);

% Derive the Line Spread Function (LSF)
lsf = diff(erfFit(x)) / pixelSize;

% Calculate FWHM of the LSF
halfMax = max(lsf) / 2;
fwhmIndices = find(lsf >= halfMax);
FWHM = x(fwhmIndices(end)) - x(fwhmIndices(1));

% Display results
figure;
subplot(2,1,1);
plot(x, pinkProfile, 'b-', 'LineWidth', 1.5);
hold on;
plot(x, erfFit(x), 'r--', 'LineWidth', 1.5);
title('Pink Beam Edge Spread Function (ESF)');
xlabel('Distance (microns)');
ylabel('Normalized Intensity');
legend('Profile', 'Error Function Fit');
grid on;

subplot(2,1,2);
plot(x(1:end-1), lsf, 'k-', 'LineWidth', 1.5);
title('Pink Beam Line Spread Function (LSF)');
xlabel('Distance (microns)');
ylabel('LSF Intensity');
grid on;

% Print the FWHM result
fprintf('Pink Beam FWHM of the Line Spread Function (LSF): %.2f microns\n', FWHM);

%%

% Define input and output directories
input_directory = '/data/projects/drex/PINKBEAM/ForgedAl/Rocking_live_2x/';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/Rocking_live_2x_videos/';

% File range
num_files = 150; % Number of files
file_prefix = 'rci_'; % File prefix
file_extension = '.h5'; % File extension

% Group to read from HDF5
group_path = '/entry/FWHM/FWHM';

% Conversion factor for axes
pixel_to_microns = 0.18;

% Colormap and limits
colormap_used = 'jet';
caxis_limits = [0.013, 0.04]; % Updated lower threshold to 0.013

% Create a video writer object
video_file = fullfile(output_directory, 'Rocking_live_video.avi');
v = VideoWriter(video_file, 'Motion JPEG AVI');
v.FrameRate = 1 / 5; % Frame rate = 1 frame every 5 seconds
open(v);

% Process each file
for i = 1:num_files
    % Read the HDF5 file
    file_name = sprintf('%s%d%s', file_prefix, i, file_extension);
    file_path = fullfile(input_directory, file_name);
    if ~isfile(file_path)
        fprintf('File %s not found. Skipping.\n', file_path);
        continue;
    end
    data = h5read(file_path, group_path);
    
    % Convert pixel dimensions to microns
    [num_rows, num_cols] = size(data);
    x_microns = (0:num_cols-1) * pixel_to_microns;
    y_microns = (0:num_rows-1) * pixel_to_microns;
    
    % Plot the data
    figure('Visible', 'off');
    imagesc(x_microns, y_microns, rot90(data)); % Rotate image by 90 degrees
    colormap(colormap_used);
    caxis(caxis_limits);
    cbar = colorbar; % Create colorbar object
    ylabel(cbar, 'FWHM (°)'); % Set colorbar label
    axis image;
    xlabel('X (microns)');
    ylabel('Y (microns)');
    title(sprintf('Time: %d s', i * 5), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    
    % Save the image
    output_image_file = fullfile(output_directory, sprintf('frame_%03d.png', i));
    saveas(gcf, output_image_file);
    
    % Write the frame to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Close the figure
    close(gcf);
end

% Close the video writer object
close(v);

disp('Processing complete. Video and images saved.');

%%

% Define input and output directories
input_directory = '/data/projects/drex/PINKBEAM/ForgedAl/Rocking_live_2x/';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/Rocking_live_2x_videos/';

% File range
num_files = 150; % Number of files
file_prefix = 'rci_'; % File prefix
file_extension = '.h5'; % File extension

% Group to read from HDF5
group_path = '/entry/Peak position/Peak position';

% Conversion factor for axes
pixel_to_microns = 0.18;

% Normalization midpoint
midpoint = (8.265 + 8.29) / 2;

% Colormap and limits
colormap_used = 'jet';
caxis_limits = [8.265 - midpoint, 8.29 - midpoint]; % Normalize around 0

% Create a video writer object
video_file = fullfile(output_directory, 'Rocking_live_video_peak_position.avi');
v = VideoWriter(video_file, 'Motion JPEG AVI');
v.FrameRate = 1 / 5; % Frame rate = 1 frame every 5 seconds
open(v);

% Process each file
for i = 1:num_files
    % Read the HDF5 file
    file_name = sprintf('%s%d%s', file_prefix, i, file_extension);
    file_path = fullfile(input_directory, file_name);
    if ~isfile(file_path)
        fprintf('File %s not found. Skipping.\n', file_path);
        continue;
    end
    data = h5read(file_path, group_path);
    
    % Normalize data around 0
    data_normalized = data - midpoint;
    
    % Convert pixel dimensions to microns
    [num_rows, num_cols] = size(data_normalized);
    x_microns = (0:num_cols-1) * pixel_to_microns;
    y_microns = (0:num_rows-1) * pixel_to_microns;
    
    % Plot the data
    figure('Visible', 'off');
    imagesc(x_microns, y_microns, rot90(data_normalized)); % Rotate image by 90 degrees
    colormap(colormap_used);
    caxis(caxis_limits); % Set normalized limits
    cbar = colorbar; % Create colorbar object
    ylabel(cbar, 'Local orientation in degrees'); % Set colorbar label
    axis image;
    xlabel('X (microns)');
    ylabel('Y (microns)');
    title(sprintf('Time: %d s', i * 5), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    
    % Save the image
    output_image_file = fullfile(output_directory, sprintf('frame_%03d.png', i));
    saveas(gcf, output_image_file);
    
    % Write the frame to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Close the figure
    close(gcf);
end

% Close the video writer object
close(v);

disp('Processing complete. Video and images saved.');


%%
%% Define the directory and file path
data_directory = '/data/visitor/ihma587/id03/20241202/RAW_DATA/Forged_Al/Forged_Al_heating3_ff2x_setpoint700now';

% Output directory to save the images
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Define the range of datasets (4.1 to 13.1)
start_scan = 5.1;
end_scan = 8.1;

% Generate the scan list
scan_list = start_scan:1:end_scan; % Create a list of scan numbers (4.1, 5.1, ..., 13.1)

% Define the ROI (800:1400 in X and Y)
roi_rows = 1000:1600; % Y-direction
roi_cols = 900:1500; % X-direction

% Initialize file numbering
file_counter = 4; % Start numbering from 04 (change as needed)

% Process each dataset
for scan_num = scan_list
    % Define the dataset path
    dataset_path = sprintf('/%.1f/measurement/pco_ff', scan_num);
    disp(['Reading dataset: ', dataset_path]);
    
    % Read the HDF5 file
    h5_file = fullfile(data_directory, 'Forged_Al_heating3_ff2x_setpoint700now.h5');
    try
        image_data = h5read(h5_file, dataset_path);
    catch ME
        fprintf('Failed to read dataset %s. Skipping.\n', dataset_path);
        continue;
    end
    
    % Get the dimensions of the dataset
    [N, M, num_frames] = size(image_data);
    disp(['Data Dimensions for Scan ', num2str(scan_num), ': ', num2str(N), ' x ', num2str(M), ' x ', num2str(num_frames)]);

    % Save every 5th image in the dataset
    for frame_idx = 1:5:num_frames % Process every 5th image
        % Extract the current image
        img = squeeze(image_data(:, :, frame_idx));
        
 % Apply median filtering and scale the intensity to the desired range (80 to 1600)
    img_scaled = uint8(255 * mat2gray(medfilt2(img_roi), [80, 1900]));

        % Extract the ROI (800:1400 in X and Y)
        img_roi = img(roi_rows, roi_cols);
        colormap_hot = parula(256); % Define the 'hot' colormap with 256 levels
        img_colored = ind2rgb(img_scaled, colormap_hot); % Convert the grayscale image to RGB using the colormap

    % Construct the filename (e.g., 04.png, 05.png, etc.)
    filename = sprintf('%02d.png', file_counter); % File numbering starts at 04
    output_path = fullfile(output_directory, filename);

    % Save the color-coded image as a PNG
    imwrite(img_colored, output_path);
     
        
        % Increment the file counter
        file_counter = file_counter + 1;
    end
    
    disp(['Saved every 5th image for scan ', num2str(scan_num)]);
end

disp('Processing complete. All images saved.');


%%
%% Define Input and Output Directories
input_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/'; % Directory where PNGs are saved
output_video_path = fullfile(input_directory, 'stabilized_video.avi'); % Output video file name

% List all PNG files in the directory
image_files = dir(fullfile(input_directory, '*.png'));
num_images = length(image_files);

% Check if there are images to process
if num_images == 0
    error('No PNG images found in the directory: %s', input_directory);
end

% Define Parameters
background_threshold = 80; % Background intensity threshold
roi_size = [601, 601]; % ROI size: [height, width]
frame_rate = 10; % Video frame rate

% Initialize the video writer object
save_video = true; % Set to true to save the stabilized video
if save_video
    v = VideoWriter(output_video_path, 'Motion JPEG AVI');
    v.FrameRate = frame_rate;
    open(v);
end

% Placeholder for the video frames
stabilized_frames = [];

% Process Each Image
disp('Processing images...');
for i = 1:num_images
    % Read the current image
    img_path = fullfile(input_directory, image_files(i).name);
    img = imread(img_path);
    
    % Convert to grayscale if the image is RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    % Find the brightest region (above the background threshold)
    binary_mask = img > background_threshold;
    props = regionprops(binary_mask, img, 'WeightedCentroid'); % Find intensity-weighted centroid
    
    if isempty(props)
        warning('No region found above the background threshold in frame %d. Skipping.', i);
        continue;
    end
    
    % Get the weighted centroid (stabilized center)
    centroid = round(props(1).WeightedCentroid);
    center_y = centroid(2); % y-coordinate (row)
    center_x = centroid(1); % x-coordinate (column)
    
    % Define ROI bounds
    half_height = floor(roi_size(1) / 2);
    half_width = floor(roi_size(2) / 2);
    row_start = max(center_y - half_height, 1);
    row_end = min(center_y + half_height, size(img, 1));
    col_start = max(center_x - half_width, 1);
    col_end = min(center_x + half_width, size(img, 2));
    
    % Extract ROI
    stabilized_roi = img(row_start:row_end, col_start:col_end);
    
    % Pad ROI if necessary to maintain consistent size
    padding_rows = max(0, roi_size(1) - size(stabilized_roi, 1));
    padding_cols = max(0, roi_size(2) - size(stabilized_roi, 2));
    stabilized_roi = padarray(stabilized_roi, ...
        [padding_rows, padding_cols], ...
        90, 'post'); % Pad with background value
    
    % Save the stabilized frame to the list
    stabilized_frames{i} = stabilized_roi; %#ok<AGROW>
    
    % Write the frame to the video
    if save_video
        writeVideo(v, mat2gray(stabilized_roi));
    end
end

% Close the video writer
if save_video
    close(v);
    disp(['Stabilized video saved to: ', output_video_path]);
end

% Play the Stabilized Video in MATLAB
disp('Playing stabilized video...');
figure;
for i = 1:length(stabilized_frames)
    imshow(stabilized_frames{i}, []);
    title(sprintf('Frame %d of %d', i, length(stabilized_frames)));
    pause(1 / frame_rate);
end
disp('Video playback complete.');

%%
% Define Directories
data_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';

% List all PNG files for scans 5.1 to 8.1
image_files = dir(fullfile(data_directory, '*.png'));
num_images = length(image_files);

if num_images == 0
    error('No PNG images found in the directory: %s', data_directory);
end

% Parameters
pixel_size = 0.19; % Effective pixel size in micrometers (190 nm)
time_per_image = 0.5; % Time per image in seconds
times = (1:num_images) * time_per_image; % Time vector

% Initialize Arrays to Store Object Properties
object_areas = zeros(1, num_images); % Object areas in micrometers squared
boundary_lengths = zeros(1, num_images); % Object boundary lengths in micrometers

% Process Images
disp('Analyzing images...');
for i = 1:num_images
    % Read the current image
    img_path = fullfile(data_directory, image_files(i).name);
    img = imread(img_path);
    
    % Convert to grayscale if the image is RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    % Threshold the image to isolate the object
    binary_img = imbinarize(img, 'adaptive', 'Sensitivity', 0.5); % Adaptive thresholding
    
    % Remove small noise and fill holes
    binary_img = imopen(binary_img, strel('disk', 2)); % Morphological opening to remove noise
    binary_img = imfill(binary_img, 'holes'); % Fill holes
    
    % Measure object properties
    stats = regionprops(binary_img, 'Area', 'Perimeter'); % Measure area and perimeter
    
    % Assume the largest object is the one of interest
    [~, largest_idx] = max([stats.Area]);
    object_area_pixels = stats(largest_idx).Area; % Area in pixels
    boundary_length_pixels = stats(largest_idx).Perimeter; % Perimeter in pixels
    
    % Convert to physical units (micrometers)
    object_areas(i) = object_area_pixels * (pixel_size^2); % Area in micrometers squared
    boundary_lengths(i) = boundary_length_pixels * pixel_size; % Perimeter in micrometers
end

disp('Image analysis complete.');

% Compute Boundary Velocity
% Boundary velocity is the change in boundary position over time
boundary_velocities = [0, diff(boundary_lengths) ./ diff(times)]; % micrometers per second

% Plot Results
figure;

% Plot object area vs. time
subplot(3, 1, 1);
plot(times, object_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Object Area (\mum^2)');
title('Object Area as a Function of Time');
grid on;

% Plot boundary length vs. time
subplot(3, 1, 2);
plot(times, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Length (\mum)');
title('Boundary Length as a Function of Time');
grid on;

% Plot boundary velocity vs. time
subplot(3, 1, 3);
plot(times, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Velocity (\mum/s)');
title('Boundary Velocity as a Function of Time');
grid on;

disp('Plots generated successfully.');

%%
% Define Input and Output Directories
data_directory = '/data/visitor/ihma587/id03/20241202/RAW_DATA/Forged_Al/Forged_Al_heating3_ff2x_setpoint700now';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Define the Range of Datasets and ROI
start_scan = 7.1;
end_scan = 8.1;
scan_list = start_scan:1:end_scan; % Generate a list of scan numbers (5.1, 6.1, ..., 8.1)

roi_rows = 1000:1600; % Y-direction
roi_cols = 900:1500; % X-direction

pixel_size = 0.19; % Effective pixel size in micrometers (190 nm)
time_per_image = 5; % Time per image in seconds

% Initialize Arrays for Object Properties
object_areas = []; % Object areas in micrometers squared
boundary_lengths = []; % Boundary lengths in micrometers
time_vector = []; % Time vector

% Process Each Dataset
disp('Processing datasets...');
for scan_num = scan_list
    % Define the dataset path
    dataset_path = sprintf('/%.1f/measurement/pco_ff', scan_num);
    disp(['Reading dataset: ', dataset_path]);

    % Read the HDF5 file
    h5_file = fullfile(data_directory, 'Forged_Al_heating3_ff2x_setpoint700now.h5');
    try
        image_data = h5read(h5_file, dataset_path);
    catch ME
        fprintf('Failed to read dataset %s. Skipping.\n', dataset_path);
        continue;
    end

    % Get the dimensions of the dataset
    [N, M, num_frames] = size(image_data);
    disp(['Data Dimensions for Scan ', num2str(scan_num), ': ', num2str(N), ' x ', num2str(M), ' x ', num2str(num_frames)]);

    % Process Every 5th Image in the Dataset
    for frame_idx = 1:10:num_frames
        % Extract the current image
        img = squeeze(image_data(:, :, frame_idx));

        % Extract the ROI (apply the defined rows and columns)
        img_roi = img(roi_rows, roi_cols);

        % Normalize and filter the ROI
        img_filtered = medfilt2(img_roi); % Apply median filtering
        img_normalized = mat2gray(img_filtered, [80, 1900]); % Normalize intensity to [0, 1]

        % Detect the object using intensity thresholding or edge detection
        mask = img_normalized > 0.1; % Apply an intensity threshold to create a binary mask
        mask = imopen(mask, strel('disk', 1)); % Morphological opening to remove noise
        mask = imfill(mask, 'holes'); % Fill holes in the mask

        % Optionally, use edge detection to enhance the boundary
        edges = edge(img_normalized, 'sobel');

        % Measure object properties
        stats = regionprops(mask, 'Area', 'Perimeter'); % Get object area and boundary length
        [~, largest_idx] = max([stats.Area]); % Assume the largest object is of interest
        object_area_pixels = stats(largest_idx).Area; % Object area in pixels
        boundary_length_pixels = stats(largest_idx).Perimeter; % Boundary length in pixels

        % Convert to physical units (micrometers)
        object_areas = [object_areas, object_area_pixels * (pixel_size^2)]; % Area in micrometers squared
        boundary_lengths = [boundary_lengths, boundary_length_pixels * pixel_size]; % Length in micrometers
        time_vector = [time_vector, frame_idx * time_per_image]; % Time in seconds

        % Visualize the object detection on one of the frames
        if frame_idx == 1 && scan_num == start_scan
            figure;
            subplot(1, 3, 1);
            imshow(img_normalized, []);
            title('Normalized Image');

            subplot(1, 3, 2);
            imshow(mask, []);
            title('Binary Mask');

            subplot(1, 3, 3);
            imshow(img_normalized, []);
            hold on;
            boundaries = bwboundaries(mask);
            plot(boundaries{1}(:, 2), boundaries{1}(:, 1), 'r', 'LineWidth', 2); % Overlay boundary
            title('Object Boundary Overlay');
            hold off;
        end
    end
end

% Compute Boundary Velocity
boundary_velocities = [0, diff(boundary_lengths) ./ diff(time_vector)]; % Boundary velocity in micrometers per second

% Plot Results
figure;

% Plot object area vs. time
subplot(3, 1, 1);
plot(time_vector, object_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Object Area (\mum^2)');
title('Object Area as a Function of Time');
grid on;

% Plot boundary length vs. time
subplot(3, 1, 2);
plot(time_vector, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Length (\mum)');
title('Boundary Length as a Function of Time');
grid on;

% Plot boundary velocity vs. time
subplot(3, 1, 3);
plot(time_vector, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Velocity (\mum/s)');
title('Boundary Velocity as a Function of Time');
grid on;

disp('Analysis complete.');

%%
data_directory = '/data/visitor/ihma587/id03/20241202/RAW_DATA/Forged_Al/Forged_Al_heating3_ff2x_setpoint700now';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

start_scan = 5.1;
end_scan = 8.1;
scan_list = start_scan:1:end_scan;

roi_rows = 1000:1600; % Y-direction
roi_cols = 900:1500; % X-direction

pixel_size = 0.19; % Effective pixel size in micrometers (190 nm)
time_per_image = 0.5; % Time per image in seconds

object_areas = [];
boundary_lengths = [];
time_vector = [];
boundary_overlay_images = [];
boundary_overlay_times = [];

disp('Processing datasets...');
for scan_num = scan_list
    dataset_path = sprintf('/%.1f/measurement/pco_ff', scan_num);
    disp(['Reading dataset: ', dataset_path]);

    h5_file = fullfile(data_directory, 'Forged_Al_heating3_ff2x_setpoint700now.h5');
    try
        image_data = h5read(h5_file, dataset_path);
    catch ME
        fprintf('Failed to read dataset %s. Skipping.\n', dataset_path);
        continue;
    end

    [N, M, num_frames] = size(image_data);
    disp(['Data Dimensions for Scan ', num2str(scan_num), ': ', num2str(N), ' x ', num2str(M), ' x ', num2str(num_frames)]);

    for frame_idx = 1:10:num_frames
        img = squeeze(image_data(:, :, frame_idx));
        img_roi = img(roi_rows, roi_cols);
        img_filtered = medfilt2(img_roi);
        img_normalized = mat2gray(img_filtered, [80, 1900]);

        mask = img_normalized > 0.1;
        mask = imopen(mask, strel('disk', 1));
        mask = imfill(mask, 'holes');

        edges = edge(img_normalized, 'sobel');

        stats = regionprops(mask, 'Area', 'Perimeter');
        [~, largest_idx] = max([stats.Area]);
        object_area_pixels = stats(largest_idx).Area;
        boundary_length_pixels = stats(largest_idx).Perimeter;

        object_areas = [object_areas, object_area_pixels * (pixel_size^2)];
        boundary_lengths = [boundary_lengths, boundary_length_pixels * pixel_size];
        time_vector = [time_vector, frame_idx * time_per_image];

        if length(boundary_overlay_images) < 4 && (frame_idx == 1 || mod(length(boundary_overlay_images), floor(num_frames / 4)) == 0)
            boundary_overlay_images{end + 1} = {img_normalized, mask};
            boundary_overlay_times(end + 1) = frame_idx * time_per_image;
        end
    end
end

boundary_velocities = [0, diff(boundary_lengths) ./ diff(time_vector)];

figure;
subplot(3, 1, 1);
plot(time_vector, object_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Object Area (\mum^2)');
title('Object Area as a Function of Time');
grid on;

subplot(3, 1, 2);
plot(time_vector, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Length (\mum)');
title('Boundary Length as a Function of Time');
grid on;

subplot(3, 1, 3);
plot(time_vector, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Velocity (\mum/s)');
title('Boundary Velocity as a Function of Time');
grid on;

disp('Analysis complete.');

figure;
for i = 1:length(boundary_overlay_images)
    subplot(2, 2, i);
    imshow(boundary_overlay_images{i}{1}, []);
    hold on;
    boundaries = bwboundaries(boundary_overlay_images{i}{2});
    plot(boundaries{1}(:, 2), boundaries{1}(:, 1), 'r', 'LineWidth', 2);
    title(sprintf('Boundary Overlay at %.1f s', boundary_overlay_times(i)));
    hold off;
end

%% grain statistics

data_directory = '/data/visitor/ihma587/id03/20241202/RAW_DATA/Forged_Al/Forged_Al_heating3_ff2x_setpoint700now';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

start_scan = 5.1;
end_scan = 8.1;
scan_list = start_scan:1:end_scan;

roi_rows = 1000:1600; % Y-direction
roi_cols = 900:1500; % X-direction

pixel_size = 0.19; % Effective pixel size in micrometers (190 nm)
time_per_image = 0.5; % Time per image in seconds

grain_areas = [];
boundary_lengths = [];
boundary_smoothness = []; % Measure of grain boundary smoothness (Boundary Length / Area)
boundary_curvatures = []; % Measure of average boundary curvature
time_vector = [];
boundary_overlay_images = {}; % Store boundary overlays for all scans
boundary_overlay_times = []; % Store corresponding times for all scans

disp('Processing datasets...');

% Initialize cumulative time counter
cumulative_time = 0;

for scan_num = scan_list
    dataset_path = sprintf('/%.1f/measurement/pco_ff', scan_num);
    disp(['Reading dataset: ', dataset_path]);

    h5_file = fullfile(data_directory, 'Forged_Al_heating3_ff2x_setpoint700now.h5');
    try
        image_data = h5read(h5_file, dataset_path);
    catch ME
        fprintf('Failed to read dataset %s. Skipping.\n', dataset_path);
        continue;
    end

    [N, M, num_frames] = size(image_data);
    disp(['Data Dimensions for Scan ', num2str(scan_num), ': ', num2str(N), ' x ', num2str(M), ' x ', num2str(num_frames)]);

    for frame_idx = 1:10:num_frames
        img = squeeze(image_data(:, :, frame_idx));
        img_roi = img(roi_rows, roi_cols);
        img_filtered = medfilt2(img_roi);
        img_normalized = mat2gray(img_filtered, [80, 1900]);

        mask = img_normalized > 0.1;
        mask = imopen(mask, strel('disk', 1));
        mask = imfill(mask, 'holes');

        stats = regionprops(mask, 'Area', 'Perimeter');
        [~, largest_idx] = max([stats.Area]);
        grain_area_pixels = stats(largest_idx).Area;
        boundary_length_pixels = stats(largest_idx).Perimeter;

        % Curvature Calculation
        boundaries = bwboundaries(mask);
        boundary_coords = boundaries{largest_idx};
        x_coords = boundary_coords(:, 2);
        y_coords = boundary_coords(:, 1);

        % Use a smoothed spline to fit the boundary and calculate curvature
        pp_x = csaps(1:length(x_coords), x_coords, 0.9); % Spline fit for x-coordinates
        pp_y = csaps(1:length(y_coords), y_coords, 0.9); % Spline fit for y-coordinates
        dx = fnval(fnder(pp_x, 1), 1:length(x_coords)); % First derivative
        dy = fnval(fnder(pp_y, 1), 1:length(y_coords));
        ddx = fnval(fnder(pp_x, 2), 1:length(x_coords)); % Second derivative
        ddy = fnval(fnder(pp_y, 2), 1:length(x_coords));
        curvature = abs(dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2)); % Curvature formula
        mean_curvature = mean(curvature(~isnan(curvature))); % Average curvature

        grain_areas = [grain_areas, grain_area_pixels * (pixel_size^2)];
        boundary_lengths = [boundary_lengths, boundary_length_pixels * pixel_size];
        boundary_smoothness = [boundary_smoothness, (boundary_length_pixels * pixel_size) / (grain_area_pixels * (pixel_size^2))];
        boundary_curvatures = [boundary_curvatures, mean_curvature];

        % Add cumulative time for the image
        time_vector = [time_vector, cumulative_time + frame_idx * time_per_image];

        % Store boundary overlay for visualization
        if frame_idx == 1 || frame_idx == num_frames
            boundary_overlay_images{end + 1} = {img_normalized, mask};
            boundary_overlay_times = [boundary_overlay_times, cumulative_time + frame_idx * time_per_image];
        end
    end

    % Update cumulative time for the next scan
    cumulative_time = cumulative_time + num_frames * time_per_image;
end

boundary_velocities = [0, diff(boundary_lengths) ./ diff(time_vector)];

% Plot grain area, boundary length, velocity, smoothness, and curvature
figure;
subplot(5, 1, 1);
plot(time_vector, grain_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Grain Area (\mum^2)');
title('Grain Area as a Function of Time');
grid on;

subplot(5, 1, 2);
plot(time_vector, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Length (\mum)');
title('Boundary Length as a Function of Time');
grid on;

subplot(5, 1, 3);
plot(time_vector, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Velocity (\mum/s)');
title('Boundary Velocity as a Function of Time');
grid on;

subplot(5, 1, 4);
plot(time_vector, boundary_smoothness, 'm-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Smoothness (Length / Area)');
title('Boundary Smoothness as a Function of Time');
grid on;

subplot(5, 1, 5);
plot(time_vector, boundary_curvatures, 'c-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Average Curvature (\mum^{-1})');
title('Boundary Curvature as a Function of Time');
grid on;

disp('Analysis complete.');

%% Visualize the first and last images of each scan
figure;
num_panels = length(boundary_overlay_images); % Number of images stored (2 images per scan)
for panel_idx = 1:num_panels
    subplot(ceil(num_panels / 2), 4, panel_idx); % Create a grid with enough rows for all panels
    imshow(boundary_overlay_images{panel_idx}{1}, []);
    hold on;
    boundaries = bwboundaries(boundary_overlay_images{panel_idx}{2});
    plot(boundaries{1}(:, 2), boundaries{1}(:, 1), 'r', 'LineWidth', 2);
    title(sprintf('Cumulative Time: %.1f s', boundary_overlay_times(panel_idx)));
    hold off;
end

%%

data_directory = '/data/visitor/ihma587/id03/20241202/RAW_DATA/Forged_Al/Forged_Al_heating3_ff2x_setpoint700now';
output_directory = '/data/projects/drex/PINKBEAM/ForgedAl/timescans_700setpoint/';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

start_scan = 5.1;
end_scan = 8.1;
scan_list = start_scan:1:end_scan; % Ensure scans use floating-point numbers

roi_rows = 1000:1600; % Y-direction
roi_cols = 900:1500; % X-direction

pixel_size = 0.19; % Effective pixel size in micrometers (190 nm)
time_per_image = 0.5; % Time per image in seconds

grain_areas = [];
boundary_lengths = [];
boundary_smoothness = [];
boundary_curvatures = [];
time_vector = [];
boundary_overlay_images = {};
boundary_overlay_times = [];

disp('Processing datasets...');

% Initialize cumulative time counter
cumulative_time = 0;

for scan_num = scan_list
    % Correct dataset path with floating-point scan number formatting
    dataset_path = sprintf('/%.1f/measurement/pco_ff', scan_num);
    disp(['Reading dataset: ', dataset_path]);

    h5_file = fullfile(data_directory, 'Forged_Al_heating3_ff2x_setpoint700now.h5');
    try
        % Read HDF5 file for the specific dataset
        image_data = h5read(h5_file, dataset_path);
    catch ME
        fprintf('Failed to read dataset %s. Skipping.\n', dataset_path);
        continue;
    end

    [N, M, num_frames] = size(image_data);
    disp(['Data Dimensions for Scan ', num2str(scan_num), ': ', num2str(N), ' x ', num2str(M), ' x ', num2str(num_frames)]);

    for frame_idx = 1:10:num_frames
        img = squeeze(image_data(:, :, frame_idx));
        img_roi = img(roi_rows, roi_cols);
        img_filtered = medfilt2(img_roi);
        img_normalized = mat2gray(img_filtered, [80, 1900]);

        mask = img_normalized > 0.1;
        mask = imopen(mask, strel('disk', 1));
        mask = imfill(mask, 'holes');

        stats = regionprops(mask, 'Area', 'Perimeter');
        if isempty(stats)
            continue; % Skip if no valid regions are detected
        end

        [~, largest_idx] = max([stats.Area]);
        grain_area_pixels = stats(largest_idx).Area;
        boundary_length_pixels = stats(largest_idx).Perimeter;

        % Curvature Calculation
        boundaries = bwboundaries(mask);
        boundary_coords = boundaries{largest_idx};
        x_coords = boundary_coords(:, 2);
        y_coords = boundary_coords(:, 1);

        % Use a smoothed spline to fit the boundary and calculate curvature
        pp_x = csaps(1:length(x_coords), x_coords, 0.9); % Spline fit for x-coordinates
        pp_y = csaps(1:length(y_coords), y_coords, 0.9); % Spline fit for y-coordinates
        dx = fnval(fnder(pp_x, 1), 1:length(x_coords)); % First derivative
        dy = fnval(fnder(pp_y, 1), 1:length(y_coords));
        ddx = fnval(fnder(pp_x, 2), 1:length(x_coords)); % Second derivative
        ddy = fnval(fnder(pp_y, 2), 1:length(x_coords));
        curvature = abs(dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2)); % Curvature formula
        mean_curvature = mean(curvature(~isnan(curvature))); % Average curvature

        grain_areas = [grain_areas, grain_area_pixels * (pixel_size^2)];
        boundary_lengths = [boundary_lengths, boundary_length_pixels * pixel_size];
        boundary_smoothness = [boundary_smoothness, (boundary_length_pixels * pixel_size) / (grain_area_pixels * (pixel_size^2))];
        boundary_curvatures = [boundary_curvatures, mean_curvature];

        % Add cumulative time for the image
        time_vector = [time_vector, cumulative_time + frame_idx * time_per_image];

        % Store boundary overlay for visualization
        if frame_idx == 1 || frame_idx == num_frames
            boundary_overlay_images{end + 1} = {img_normalized, mask};
            boundary_overlay_times = [boundary_overlay_times, cumulative_time + frame_idx * time_per_image];
        end
    end

    % Update cumulative time for the next scan
    cumulative_time = cumulative_time + num_frames * time_per_image;
end

boundary_velocities = [0, diff(boundary_lengths) ./ diff(time_vector)];

% Plot grain area, boundary length, velocity, smoothness, and curvature
figure;
subplot(4, 1, 1);
plot(time_vector, grain_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Grain Area (\mum^2)');
title('Grain Area as a Function of Time');
grid on;

subplot(4, 1, 2);
plot(time_vector, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Length (\mum)');
title('Boundary Length as a Function of Time');
grid on;

subplot(4, 1, 3);
plot(time_vector, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Boundary Velocity (\mum/s)');
title('Boundary Velocity as a Function of Time');
grid on;

%subplot(4, 1, 3);
%plot(time_vector, boundary_smoothness, 'm-', 'LineWidth', 2);
%xlabel('Time (s)');
%ylabel('Boundary Smoothness (Length / Area)');
%title('Boundary Smoothness as a Function of Time');
%grid on;

subplot(4, 1, 4);
plot(time_vector, boundary_curvatures, 'c-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Av. Curvature (\mum^{-1})');
title('Boundary Curvature as a Function of Time');
grid on;

disp('Analysis complete.');

% Visualize the first and last images of each scan
figure;
num_panels = length(boundary_overlay_images);
for panel_idx = 1:num_panels
    subplot(ceil(num_panels / 4), 4, panel_idx);
    imshow(boundary_overlay_images{panel_idx}{1}, []);
    hold on;
    boundaries = bwboundaries(boundary_overlay_images{panel_idx}{2});
    plot(boundaries{1}(:, 2), boundaries{1}(:, 1), 'r', 'LineWidth', 2);
    title(sprintf('Cumulative Time: %.1f s', boundary_overlay_times(panel_idx)));
    hold off;
end

%%
% Define aspect ratio for single-column A4 format
fig_width = 85 / 25.4;  % Convert mm to inches (~3.35 inches)
fig_height = 110 / 25.4; % Adjusted height (~4.33 inches)

% First Figure: Grain Area & Boundary Length
figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);

subplot(2, 1, 1);
plot(time_vector, grain_areas, 'b-', 'LineWidth', 2);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 15);
ylabel('Grain Area (\mum^2)', 'FontName', 'Helvetica', 'FontSize', 15);
title('Grain Area as a Function of Time', 'FontName', 'Helvetica', 'FontSize', 15);
set(gca, 'FontName', 'Helvetica', 'FontSize', 15);
grid on;

subplot(2, 1, 2);
plot(time_vector, boundary_lengths, 'r-', 'LineWidth', 2);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 15);
ylabel('Boundary Length (\mum)', 'FontName', 'Helvetica', 'FontSize', 15);
title('Boundary Length as a Function of Time', 'FontName', 'Helvetica', 'FontSize', 15);
set(gca, 'FontName', 'Helvetica', 'FontSize', 15);
grid on;

% Second Figure: Boundary Velocity & Boundary Curvature
figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);

subplot(2, 1, 1);
plot(time_vector, boundary_velocities, 'g-', 'LineWidth', 2);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 15);
ylabel('Boundary Velocity (\mum/s)', 'FontName', 'Helvetica', 'FontSize', 15);
title('Boundary Velocity as a Function of Time', 'FontName', 'Helvetica', 'FontSize', 15);
set(gca, 'FontName', 'Helvetica', 'FontSize', 15);
grid on;

subplot(2, 1, 2);
plot(time_vector, boundary_curvatures, 'c-', 'LineWidth', 2);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 15);
ylabel('Av. Curvature (\mum^{-1})', 'FontName', 'Helvetica', 'FontSize', 15);
title('Boundary Curvature as a Function of Time', 'FontName', 'Helvetica', 'FontSize', 15);
set(gca, 'FontName', 'Helvetica', 'FontSize', 15);
grid on;


