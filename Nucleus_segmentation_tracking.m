function Nucleus_segmentation_tracking

clear 
close all
clc

 format short
 warning('off','all');  
 
 % loading the data
 full_path = 'C:\Users\Debojyoti\Desktop\New folder\Video_Result_p21_RPE\mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)\H2_image.mat';
 D3  = load(full_path);

 H2_image = D3.H2_image;
 
 % creating video file
 video_name = strcat(full_path,'.mp4');
 vidObj = VideoWriter('H2','MPEG-4');
 open(vidObj);
 
for i = 1:1:50 %size(H2_image,3)
fprintf('frame = %d\n',i)  
    
   nucleus = H2_image(:,:,i);
   
   [BW, nucleus1, stats] = nuclei(nucleus); 
   

   
   imshow(imadjust(nucleus1))
   hold on
   title(['frame = ', num2str(i)])
   

   
   % Initialization of storage vectors/matrices
    area = zeros(length(stats),1);
    TMP_centroid = zeros(length(stats),2);
    TMP_cov_r = zeros(length(stats),1);
    
    vi = 1;
    
    clear valid_index BB_valid Boundary_Valid Area_Valid Centroid_valid Cov_r_valid  
    
    
    
    for k = 1:length(stats)

        area(k) = stats(k).Area;
        perimeter(k) = stats(k).Perimeter;
        major_axis_length(k) = stats(k).MajorAxisLength;
        
        TMP_BB = stats(k).BoundingBox + 3*[-1 -1 2 2];
        TMP_centroid(k,:) = stats(k).Centroid;
      
        TMP_BW = imcrop(BW,TMP_BB);
        [TMP_B,~] = bwboundaries(TMP_BW,'noholes');

        TMP_enclosed_area = zeros(1,length(TMP_B));
        
        for j = 1:length(TMP_B)
             TMP_boundary = TMP_B{j};
             TMP_enclosed_area(j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
        end
        
        [~,TMP_I_max] = max(TMP_enclosed_area);
        TMP_boundary_valid = TMP_B{TMP_I_max};
        
        TMP_r = zeros(1,length(TMP_boundary_valid));
            for b_idx = 1:length(TMP_boundary_valid)
                TMP_r(b_idx) = sqrt((TMP_boundary_valid(b_idx,2)+ TMP_BB(1)-TMP_centroid (k,1))^2 + (TMP_boundary_valid(b_idx,1)+ TMP_BB(2)-TMP_centroid (k,2))^2);
            end
        TMP_cov_r(k) = cov(TMP_r);
        
% the small areas below a certain threshold need to be ignored        
        if area(k)>7500   
            valid_index(vi) = k;
            BB_valid(vi,:) = TMP_BB;
            Boundary_Valid{vi} = TMP_boundary_valid;
            
            
            Area_Valid(vi) = area(k);
            Perimeter_Valid(vi) = perimeter(k);
            Major_Axis_Length_Valid(vi) = major_axis_length(k);
            

            Centroid_valid(vi,:) = TMP_centroid(k,:);
            Cov_r_valid(vi) = TMP_cov_r(k);
            
            vi = vi+1;
        end   
    end
    
    
  % storing all the the information together
  Data.CENTROID_tracking{i} = Centroid_valid;
  Data.AREA_tracking{i} = Area_Valid;
  Data.PERIMETER_tracking{i} = Perimeter_Valid;
  Data.Major_Axis_Length_tracking{i} = Major_Axis_Length_Valid;
  

  Data.index{i} = (1:1:vi-1);
  Data.bounding_box{i} = BB_valid;
  Data.boundary{i} = Boundary_Valid;
  

  



    if i == 1
        
        for index = 1:length(valid_index)
        
    
            boundary = Boundary_Valid{index};
            plot(boundary(:,2) + BB_valid(index,1),boundary(:,1)+ BB_valid(index,2), 'r', 'LineWidth', 2)
            plot(Centroid_valid(index,1),Centroid_valid(index,2), 'm+', 'LineWidth', 2,'Markersize',8)
    
            text(boundary(1,2) + BB_valid(index,1),boundary(1,1) + BB_valid(index,2),num2str((index)),'Color','r','fontweight','b','fontsize',10,'BackgroundColor','w')
            print(['frame',num2str(i)],'-djpeg')
            
            v1 = Data.CENTROID_tracking{i};
            old_indices = 1:1:length(v1);
            
        Data.Tracking_index{i} = old_indices;    
        end
           
           
    else


            
        v1 = Data.CENTROID_tracking{i-1};
        v2 = Data.CENTROID_tracking{i};
        a1 = Data.AREA_tracking{i-1};
        a2 = Data.AREA_tracking{i};
        p1 = Data.PERIMETER_tracking{i-1};
        p2 = Data.PERIMETER_tracking{i};
        MAL_1 = Data.Major_Axis_Length_tracking{i-1};
        MAL_2 = Data.Major_Axis_Length_tracking{i};


        
        [new_indices] = tracking_cell(v1,v2,a1,a2,old_indices, p1, p2, MAL_1,MAL_2);
        
        
        
        for index = 1:length(valid_index)
        
    
            boundary = Boundary_Valid{index};
            plot(boundary(:,2) + BB_valid(index,1),boundary(:,1)+ BB_valid(index,2), 'r', 'LineWidth', 2)
            plot(Centroid_valid(index,1),Centroid_valid(index,2), 'm+', 'LineWidth', 2,'Markersize',8)
    
            text(boundary(1,2) + BB_valid(index,1),boundary(1,1) + BB_valid(index,2),num2str(new_indices(index)),'Color','r','fontweight','b','fontsize',10,'BackgroundColor','w')
            print(['frame',num2str(i)],'-djpeg')


                    
        end
        
  
           
        old_indices = new_indices;
      Data.Tracking_index{i} = old_indices;    
    end
writeVideo(vidObj, getframe(gcf));      
end
close(vidObj);
save('DATA','Data');

end

%% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [tmpBW_er, nucleus1, stats] = nuclei(nucleus) 
% separtes the connected regions 
bg = imopen(nucleus,strel('disk',150));
nucleus1 = nucleus-bg;



T = adaptthresh(nucleus,0.5); % this is an extension of Wellner's Method (Wellner 1993)
BW = imbinarize(nucleus,T);
BW = bwareaopen(BW,500); % removes all connected components  that have fewer than 500 pixels
BW = imclearborder(BW); % removes structures that are lighter than their surroundings and that are connected to the image border
BW = imfill(BW,'holes'); % fills holes in the input binary image

tmpBW_er = BW;

  for i = 1:1:40 % image  errosion to separate wealy connected components 
   tmpBW_er = imerode(BW,strel('sphere',5));
   tmpBW_er = bwareaopen(tmpBW_er,10);
  end


%  measurements for the set of properties : 'Area','BoundingBox'....   for each 8-connected component in the binary image
stats = regionprops(tmpBW_er,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter','MajorAxisLength');

end

%%
function [new_indices] = tracking_cell(v1,v2,a1,a2,old_indices, p1, p2, MAL_1,MAL_2)
% tracking centroids of the cells 
% the comparison of the  two frames based on closest centroid and other
% conditions for example 1. relative change in the area of the cell between
% consequitive frame above certain threshold will be considered as the new
% entry/merging event; 2. proximity of the correlated centroids etc.


new_indices = 1:1:length(v2);


% initialization
 dist = zeros(length(v2),length(v1));
 m = zeros(1,length(v2));
 next_m = zeros(1,length(v2));
 index = zeros(1,length(v2));
 next_index = zeros(1,length(v2));
 
for i = 1:length(v2)
    for j = 1:length(v1)
        
       dist(i,j) = sqrt((v2(i,1)-v1(j,1))^2+(v2(i,2)-v1(j,2))^2);
       
   end
    m(i) = min(dist(i,:));
    index(i) = find(dist(i,:) == m(i));
    
    next_m(i) = min(setdiff(dist(i,:),min(dist(i,:))));
    next_index(i) = find(dist(i,:) == next_m(i));
end


unique_sorted_index = unique(index);
unique_sorted_index_m = [unique_sorted_index,zeros(1,length(old_indices)-length(unique_sorted_index))];

%     missing_cell_no = nnz(old_indices-unique_sorted_index_m);
%     missing_cell_index = find((old_indices-unique_sorted_index_m)>0);

tmp = max(old_indices);

for i = 1:1:length(unique_sorted_index)
    tmp_index = find(index == unique_sorted_index(i));
    

    if length(tmp_index)== 1 || any(abs(diff(m(tmp_index)))>150) == 1  
    
    [~,min_index] = min(m(tmp_index));

    new_indices(tmp_index(min_index)) = old_indices(unique_sorted_index(i));
    
    if length(tmp_index)>1
        tmp_index_m = setdiff(tmp_index,tmp_index(min_index));
        for j = 1:1:length(tmp_index_m)
        new_indices(tmp_index_m(j))= tmp+1;
        tmp = tmp+1;
        end
    end
    
    else
        
    if any(m(tmp_index)>110)==1 

        new_indices(tmp_index) = tmp+(1:1:length(tmp_index));
        tmp = tmp+length(tmp_index);
    else    
        [~,min_index] = min(m(tmp_index));
        new_indices(tmp_index(min_index)) = old_indices(unique_sorted_index(i));
        
        tmp_index_m = setdiff(tmp_index,tmp_index(min_index));
        new_indices(tmp_index_m) = old_indices(next_index(i));
        
    end
    end
end

tmp_index = max(union(old_indices,new_indices));

for i = 1:length(new_indices)
   if m(i)>200 
       if isempty(find(old_indices == new_indices(i), 1))~=1
            new_indices(i)= tmp_index+1;
            tmp_index = tmp_index+1;
       end
   end    
end

%% checking the change of area

C = intersect(old_indices,new_indices);

tmp_index = max(union(old_indices,new_indices))+1;
for i = 1:1:length(C)
    old_area_index = find(old_indices ==  C(i));
    new_area_index = find(new_indices ==  C(i));
    
    rel_change_area = 100*abs(a2(new_area_index)-a1(old_area_index))/a1(old_area_index);
    
    if rel_change_area > 40
        
      new_indices(find(new_indices ==  C(i))) = tmp_index;
      tmp_index = tmp_index+1;
    end
end




% to break any tie
unique_sorted_newindex = unique(new_indices);
tmp_index = max(union(old_indices,new_indices))+1;
for i = 1:1:length(unique_sorted_newindex)
    tie_index = find(new_indices == unique_sorted_newindex(i));
    if length(tie_index)>1

        [~,min_index] = min(m(tie_index));
   
    
    tie_index_m = setdiff(tie_index,tie_index(min_index));
        for j = 1:1:length(tie_index_m)
        new_indices(tie_index_m(j))= tmp_index;
        tmp_index = tmp_index + 1;
        end
        
    end
end    
end
