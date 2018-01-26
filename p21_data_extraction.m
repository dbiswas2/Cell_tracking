function p21_data_extraction



 format short
 warning('off','all');  
 
 % loading the data
 full_path1 = 'C:\Users\Debojyoti\Desktop\New folder\Video_Result_p21_RPE\mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)\H2_image.mat';
 data1  = load(full_path1);
 H2_image = data1.H2_image;
 
 full_path2 = 'C:\Users\Debojyoti\Desktop\New folder\Video_Result_p21_RPE\mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)\P21_image.mat';
 data2  = load(full_path2);
 P21_image = data2.P21_image;
 
 full_path3 = 'C:\Users\Debojyoti\Desktop\New folder\Video_Result_p21_RPE\mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)\DIC_image.mat';
 data3  = load(full_path3);
 DIC_image = data3.DIC_image;
 
 full_path = 'C:\Users\Debojyoti\Desktop\New folder\DATA.mat';
 data = load(full_path);
 DATA = data.Data;
 
 
 
 for i = 1:1:50 %size(H2_image,3)
 
     fprintf('frame = %d\n',i)  
    
   nucleus = H2_image(:,:,i);
   p21 = P21_image(:,:,i);
   DIC = DIC_image(:,:,i);
   
   bg_nucleus = imopen(nucleus,strel('disk',150)); %  morphological opening on the binary image  with the structuring element: disk type
   nucleus1 = nucleus-bg_nucleus; % background subtracted image
   
   bg_p21 = imopen(p21,strel('disk',50));
   p21_1 = p21-bg_p21; % background subtracted image
   
   

% user specific background selection
    Background = roipoly(imadjust(DIC));
    Background_masked_nucleus  = (nucleus1.*uint16(Background));
    Background_masked_p21  = (p21_1.*uint16(Background));
    
    s_background = sum(Background_masked_p21(:)); 
    
     % computation of the area cropped
        T_background = adaptthresh(Background_masked_nucleus,0.5);
        background_binarized = imbinarize(Background_masked_nucleus,T_background);
        background_tmp = bwareaopen(background_binarized,500);
        background_tmp = imclearborder(background_tmp);
        background_tmp = imfill(background_tmp,'holes');
 
        stats_bg = regionprops(background_tmp,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
        area_bg = stats_bg.Area;
        length(stats_bg);
        Result.I_bg_avg(i) = s_background/area_bg;
%-----------------------------------------------------------------------------------------------------------
   figure('color','white')
   imshow(imadjust(p21_1))
   hold on

   Boundary = DATA.boundary{i};
   Bounding_box = DATA.bounding_box{i};
   Tracking_index = DATA.Tracking_index{i};
   Area = DATA.AREA_tracking{i};
   
   clear Avg_I
   Avg_I = zeros(1,length(Boundary));
   for j = 1:length(Boundary)
       boundary = Boundary{j};
       
        mask = uint16(poly2mask(boundary(:,2),boundary(:,1),size(nucleus1,1),size(nucleus1,2)));
        image_masked = (p21_1.*mask);
        s = sum(image_masked(:));
        Avg_I(j) = s/Area(j);
       
       
       plot(boundary(:,2) + Bounding_box(j,1),boundary(:,1)+ Bounding_box(j,2), 'r', 'LineWidth', 2)
       text(boundary(1,2) + Bounding_box(j,1)-30,boundary(1,1) + Bounding_box(j,2)+30,num2str(Tracking_index(j)),'Color','y','fontweight','b','fontsize',12,'BackgroundColor','none')
   end
   
   Result.AVG_I{i} = Avg_I; % average intensity in the nuclei
   Result.AVG_I_bg{i} = Avg_I/Result.I_bg_avg(i); % fold change compare to background
   
 

    
    
   pause(0.5)
 end
   
 save('RESULT','Result');  

