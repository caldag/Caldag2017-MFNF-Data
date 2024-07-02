%% About image_processing.m

% This function takes an experiment recording, crops it, applies image processing and extracts swimmer position and orientation.
% Data from x-y plane image and x-z plane image are saved separately.
%
% This code can be transformed into a function by uncommenting lines 25 and 301 and commenting lines 29-37.
%
%% --- INPUTS TO THE FUNCTION --- %%
%%%% filepath=Destination of the folder the experiment recording is placed in (string)
%%%% freq=Rotation rate of tha magnetic field which is also the name of the recording (integer).
%%%% rectxy = [xmin ymin width height] --> Channel boundaries of x-y plane image of the swimmer, in pixels (1x4 vector)
%%%% rectxz = [xmin ymin width height] --> Channel boundaries of x-z olane image of the swimmer, in pixels (1x4 vector)
%%%% upl --> unit pixel length in mm (double)
%%%% orlimitxy --> Binary image threshold to extract tail profile from x-y plane image (double, between 0 and 1)
%%%% orlimitxz --> Binary image threshold to extract tail profile from x-z plane image (double, between 0 and 1)
% Lower orlimit values mean the intensity threshold for selecting the pixels are high (fewer pixels are selected)
%%
%
%% --- OUTPUTS OF THE FUNCTION --- %%
% The function as is does not return anything.
% Position and orientation information is exported as .dat files.
% Position and orientation data from x-y plane images are exported into .dat file
% with "xy" suffix.
% Position and orientation data from x-z plane images are exported into
% .dat file with "xz" suffix.

%The number before "xy" or "xz" suffixes denote rotation rate.
%%

%% Initialization

%function [] = image_processing(filepath, freq, rectxy, rectxz, upl, orlimitxy, orlimitxz) %Uncomment this line to call image_processing as a function

%Part to modify the code to call it directly starts here. Default values are for D1.6-L4-Q0 at 15 Hz rotation rate.

filepath=['.\D1.6-L4-Q0\']; %Set the destination folder the movie is placed in.
%The default filepath can be used if this code is placed in Experiment
%Videos folder available to download from microswimmer.sabanciuniv.edu
freq=15; %Set the frequency, i.e., name of the video file (without extension).
rectxy=[364 217 137 22];
rectxz=[364 279 134 23];
upl=20/(501-363); %20 is the length of the mm paper at the background. Pixels 501 and 363 correspond to the x- coordinate of the sides of the paper.
orlimitxy=0.75;
orlimitxz=0.78;

%Part to modify the code to call it directly ends here.

fpref = strcat(filepath, num2str(freq))  %Full filepath is formed

clearvars -except filepath fpref freq rectxy rectxz upl orlimitxy orlimitxz; %clear everything except the inputs.

READ = 1; 													  % READ = 1 if you want to read file
											
if READ
    swimmerObj = VideoReader(strcat(fpref,'.mp4')); % Video name
    frames = read(swimmerObj);	%Read video and save the frames in variable "frames"
end

OUT = 1;	% OUT = 1 if one desires to display data extraction process, as shown in the article.

if OUT
    out = VideoWriter('test.avi'); %prepere data extraction process video file
    open(out);
end

k1 = 1;  % Starting frame
a=size(frames);
k2 = a(4); % Ending frame (last frame)

rect = rectxy; % The channel boundaries
H = rect(4); % Height of the channel in pixels

X1 = 5; %Parameters that define the search region around the swimmer body
X2 = 30; %The parameters here are for D1.6-L4 configuration.
Y1 = 10;
Y2 = 10;

orlimit = orlimitxy; % Image black & white limit to extract tail profile
poslimit = 0.05; % Image black & white limit to detect the position of the head

WRITE = 1;  %Boolean to write extracted position and orientation data as .dat file

if k1==k2; WRITE = 0; end %If one frame is to be processed, do not save extracted data (for testing)

fref = imcrop(frames(:,:,:,k1),rect); % Cropping the reference (first) frame
fref1 = rgb2gray(fref); % Black & white version of reference frame
fref2 = imadjust(fref1); % Increase contrast
figure(1)
clf
clear cent cent2 theta dir xsw ysw
xsw = zeros(k2-k1,1); % Space allocation for x-, y- (or z-) position and orientation data
ysw = zeros(k2-k1,1);
theta = zeros(k2-k1,1);
dir =  zeros(k2-k1,1);

j = 0; % Loop counter

%% Frame by frame processing loop for x-y plane images

for k = k1:k2
    a = imcrop(frames(:,:,:,k),rect); % Crop the image down to channel boundaries
    subplot(411);imshow(a); % Place the cropped image onto top of figure
    b1 = rgb2gray(a); % Obtain grayscale image
    b2 = imadjust(b1); % Increase contrast
    b5 = imhistmatch(b2,fref2); % Histogram matching to reference frame
    subplot(412);imshow(b5) % Display the modified image
    c = im2bw(b2,orlimit);  % Binary image generation (for extraction of orientation)
    c2 = im2bw(b2,poslimit);  % Binary image generation (for extraction of orientation)
    subplot(413);imshow(c); % Show the black and white image for orientation extraction
    j = j+1; % Increase loop counter by one
    
    if rect == rectxy % If x-y plane is being extracted
        cc = flipud(~c);
        cc2 = flipud(~c2); % Flip "not" of the images in transformation from image coordinate system to Cartesian
    elseif rect == rectxz
        cc = ~c;
        cc2 = ~c2; % Take not of the images only as the mirrorred image is upside down.
    else
        disp(['rectangle is not specified!!']);
    end
    
    s = regionprops(cc,'Centroid','Orientation','FilledArea'); % Find the centroid, orientation and total filled area of black and white image
    s2 = regionprops(cc2,'Centroid','FilledArea'); % Find the centroid and total filled area of black and white image
    [famax,kmax] = max(cell2mat({s.FilledArea})); % Find the largest area, which should be the head
    cents = cell2mat({s.Centroid}); % Find the center of that area
    cent(j,:) = cents(2*kmax-1:2*kmax); % Record the position values
    cents2 = cell2mat({s2.Centroid}); % Find the centroid from the other image as well
    [famax2,kmax2] = max(cell2mat({s2.FilledArea}));
    cent2(j,:) = cents2(2*kmax2-1:2*kmax2);
    
    [ii,jj,ss] = find(cc); %Collect the points
    [m,n] = size(cc);
    jsw = find(jj > cent2(j,1)-X1 & jj < cent2(j,1)+X2 & ...  % Mask the region where the swimmer is
        ii < cent2(j,2)+Y2 & ii > cent2(j,2)-Y1 );         % approximately
    p = polyfit(jj(jsw),ii(jsw),1); % Fit a line to the points
    
    
    if p(1) > 0					% Bounding box is modified to find the tail pixels for extreme tilts
        Y2b = Y2+p(1)*X2; Y1b = Y1;
    else
        Y1b = Y1-p(1)*X2; Y2b = Y2;
    end
    
    clear jsw p
    
    jsw = find(jj > cent2(j,1)-X1 & jj < cent2(j,1)+X2 & ...  % Mask the region where the swimmer is
        ii < cent2(j,2)+Y2b & ii > cent2(j,2)-Y1b );         % approximately (recalculation)
    p = polyfit(jj(jsw),ii(jsw),1);	%The recalculation is carried out in case the swimmer is more tilted than anticipated by the masking region
    
    theta(j) = atan(p(1));	% Find orientation angle (either theta_xy or theta_xz)
    subplot(414)
    iic = ii(jsw);
    
    plot([0 rect(3)+1 rect(3)+1 0 0],[0 0 rect(4)+1 rect(4)+1 0],'k-',... % Plot the data extraction process
        jj(jsw),iic,'.',...
        jj(jsw),p(2)+jj(jsw)*p(1),cent(j,1),cent(j,2),'o',cent2(j,1),cent2(j,2),'*');
    axis('equal')
    
    dir(j) = p(1); % Slope of the fit line
    
    if OUT; %If output video is to be recorded, add the figure
        frm = getframe(gcf);
        writeVideo(out,frm);
    end
end

if OUT; close(out); end % Close data extraction video

xsw = cent2(:,1)*upl; % Transformation from pixel coordinates to Cartesian
ysw = (cent2(:,2)-(H+1)/2)*upl; % y- axis is along the center of the bounding box, define Cartesian coordinates accordingly
thdg = theta*180/pi; % Orientation angle in degrees

if (k2>k1)
    figure(3); %Plot extracted data
    subplot(211);
    plot(xsw,ysw,'.');axis('equal');xlabel('xsw [mm]'); ylabel('{ysw,zsw} [mm]');
    subplot(212)
    plot(theta,'.');xlabel('index'); ylabel('orientation angle [deg]')
    tmp= [xsw';ysw';thdg';dir']';
    if rect == rectxy & WRITE %Record extracted data
        save(strcat(fpref,'xy.dat'),'tmp','-ascii');
    elseif rect == rectxz & WRITE
        save(strcat(fpref,'xz.dat'),'tmp','-ascii');
    end
end

%Extraction for x-y plane image is complete.

%% Initialization for data extraction in x-z plane images.

rect = rectxz; % The channel boundaries
H = rect(4); % Height of the channel in pixels

X1 = 5; %Parameters that define the search region around the swimmer body
X2 = 30; %The parameters here are for D1.6-L4 configuration.
Y1 = 10;
Y2 = 10;

orlimit = orlimitxy; % Image black & white limit to extract tail profile
poslimit = 0.05; % Image black & white limit to detect the position of the head

WRITE = 1;  %Boolean to write extracted position and orientation data as .dat file

if k1==k2; WRITE = 0; end %If one frame is to be processed, do not save extracted data (for testing)

fref = imcrop(frames(:,:,:,k1),rect); % Cropping the reference (first) frame
fref1 = rgb2gray(fref); % Black & white version of reference frame
fref2 = imadjust(fref1); % Increase contrast
figure(1)
clf
clear cent cent2 theta dir xsw ysw
xsw = zeros(k2-k1,1); % Space allocation for x-, y- (or z-) position and orientation data
ysw = zeros(k2-k1,1);
theta = zeros(k2-k1,1);
dir =  zeros(k2-k1,1);

j = 0; % Loop counter

%% Frame by frame processing loop for x-z plane images

for k = k1:k2
    a = imcrop(frames(:,:,:,k),rect); % Crop the image down to channel boundaries
    subplot(411);imshow(a); % Place the cropped image onto top of figure
    b1 = rgb2gray(a); % Obtain grayscale image
    b2 = imadjust(b1); % Increase contrast
    b5 = imhistmatch(b2,fref2); % Histogram matching to reference frame
    subplot(412);imshow(b5) % Display the modified image
    c = im2bw(b2,orlimit);  % Binary image generation (for extraction of orientation)
    c2 = im2bw(b2,poslimit);  % Binary image generation (for extraction of orientation)
    subplot(413);imshow(c); % Show the black and white image for orientation extraction
    j = j+1; % Increase loop counter by one
    
    if rect == rectxy % If x-y plane is being extracted
        cc = flipud(~c);
        cc2 = flipud(~c2); % Flip "not" of the images in transformation from image coordinate system to Cartesian
    elseif rect == rectxz
        cc = ~c;
        cc2 = ~c2; % Take not of the images only as the mirrorred image is upside down.
    else
        disp(['rectangle is not specified!!']);
    end
    
    s = regionprops(cc,'Centroid','Orientation','FilledArea'); % Find the centroid, orientation and total filled area of black and white image
    s2 = regionprops(cc2,'Centroid','FilledArea'); % Find the centroid and total filled area of black and white image
    [famax,kmax] = max(cell2mat({s.FilledArea})); % Find the largest area, which should be the head
    cents = cell2mat({s.Centroid}); % Find the center of that area
    cent(j,:) = cents(2*kmax-1:2*kmax); % Record the position values
    cents2 = cell2mat({s2.Centroid}); % Find the centroid from the other image as well
    [famax2,kmax2] = max(cell2mat({s2.FilledArea}));
    cent2(j,:) = cents2(2*kmax2-1:2*kmax2);
    
    [ii,jj,ss] = find(cc); %Collect the points
    [m,n] = size(cc);
    jsw = find(jj > cent2(j,1)-X1 & jj < cent2(j,1)+X2 & ...  % Mask the region where the swimmer is
        ii < cent2(j,2)+Y2 & ii > cent2(j,2)-Y1 );         % approximately
    p = polyfit(jj(jsw),ii(jsw),1); % Fit a line to the points
    
    
    if p(1) > 0					% Bounding box is modified to find the tail pixels for extreme tilts
        Y2b = Y2+p(1)*X2; Y1b = Y1;
    else
        Y1b = Y1-p(1)*X2; Y2b = Y2;
    end
    
    clear jsw p
    
    jsw = find(jj > cent2(j,1)-X1 & jj < cent2(j,1)+X2 & ...  % Mask the region where the swimmer is
        ii < cent2(j,2)+Y2b & ii > cent2(j,2)-Y1b );         % approximately (recalculation)
    p = polyfit(jj(jsw),ii(jsw),1);	%The recalculation is carried out in case the swimmer is more tilted than anticipated by the masking region
    
    theta(j) = atan(p(1));	% Find orientation angle (either theta_xy or theta_xz)
    subplot(414)
    iic = ii(jsw);
    
    plot([0 rect(3)+1 rect(3)+1 0 0],[0 0 rect(4)+1 rect(4)+1 0],'k-',... % Plot the data extraction process
        jj(jsw),iic,'.',...
        jj(jsw),p(2)+jj(jsw)*p(1),cent(j,1),cent(j,2),'o',cent2(j,1),cent2(j,2),'*');
    axis('equal')
    
    dir(j) = p(1); % Slope of the fit line
    
    if OUT %If output video is to be recorded, add the figure
        frm = getframe(gcf);
        writeVideo(out,frm);
    end
end

if OUT; close(out); end % Close data extraction video

xsw = cent2(:,1)*upl; % Transformation from pixel coordinates to Cartesian
ysw = (cent2(:,2)-(H+1)/2)*upl; % y- axis is along the center of the bounding box, define Cartesian coordinates accordingly
thdg = theta*180/pi; % Orientation angle in degrees

if (k2>k1)
    figure(3); %Plot extracted data
    subplot(211);
    plot(xsw,ysw,'.');axis('equal');xlabel('xsw [mm]'); ylabel('{ysw,zsw} [mm]');
    subplot(212)
    plot(theta,'.');xlabel('index'); ylabel('orientation angle [deg]')
    tmp= [xsw';ysw';thdg';dir']';
    if rect == rectxy & WRITE %Record extracted data
        save(strcat(fpref,'xy.dat'),'tmp','-ascii');
    elseif rect == rectxz & WRITE
        save(strcat(fpref,'xz.dat'),'tmp','-ascii');
    end
end

% end