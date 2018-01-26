function p21_tracking_plot

clear
close all
clc

path1 = 'C:\Users\Debojyoti\Desktop\New folder\Video_Result_p21_RPE\';
path2 = 'mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)';

full_path1 = [path1, path2, '\DATA.mat']; % contains tracking indices
data1 = load(full_path1);
DATA = data1.Data;
 
full_path2 = [path1, path2, '\RESULT.mat']; % contains p21 data
data2 = load(full_path2);
RESULT = data2.Result;
 
 
 Tracking_index = DATA.Tracking_index{1};
 Average_intensity = RESULT.AVG_I{1}; % average intensity 
 Average_intensity_bg = RESULT.AVG_I_bg{1}; % fold increase in the average intensity compare to background
 
 for j = 1:length(Tracking_index)
      Average_intensity_tracking{j} = Average_intensity(j);
      Average_intensity_tracking_bg{j} = Average_intensity_bg(j);
      
end
 



 old_tracking_index = Tracking_index;
 frame_newentries{1} = old_tracking_index;
 for i = 2:1:50
     
     Tracking_index = DATA.Tracking_index{i};
     Average_intensity = RESULT.AVG_I{i};
     Average_intensity_bg = RESULT.AVG_I_bg{i};
     
     for j = 1:length(old_tracking_index)
         index = find(Tracking_index == old_tracking_index(j));
         tmp_content1 = [Average_intensity_tracking_bg{old_tracking_index(j)}  Average_intensity_bg(index)];
         Average_intensity_tracking_bg{old_tracking_index(j)} = tmp_content1;   
         
         tmp_content2 = [Average_intensity_tracking{old_tracking_index(j)}  Average_intensity(index)];
         Average_intensity_tracking{old_tracking_index(j)} = tmp_content2; 
     end
     
     new_entries = setdiff(Tracking_index,intersect(old_tracking_index,Tracking_index));

     frame_newentries{i} = setdiff(new_entries,frame_newentries{i-1});
     
     
     for k = 1:length(new_entries)
         index = find(Tracking_index == new_entries(k));
         Average_intensity_tracking_bg{new_entries(k)} = Average_intensity_bg(index);
         Average_intensity_tracking{new_entries(k)} = Average_intensity(index);
     end
     old_tracking_index = Tracking_index;
 end
 
 
 colors = jet(length(Average_intensity_tracking_bg));
 figure('color','white')
 set(gca,'linewidth',1,'fontweight','b')
 hold on
for i = 1:1:length(Average_intensity_tracking_bg)

   
    for k = 1:length(frame_newentries)
     if any(find(frame_newentries{1,k}== i)) == 1
         frame_index = k;
     end
    end
    
    y1 = Average_intensity_tracking_bg{i};
    y2 = Average_intensity_tracking{i};
    x = frame_index:1:frame_index+length(y1)-1;

%     h1 = plot(x,y1,'.--','MarkerSize',10);
%     set(h1,'Color',colors(i,:))
    hold on
    if length(y1) > 1
    xq = frame_index:0.5:frame_index+length(y1)-1;
    yq = interp1(x,y1,xq,'spline');
    h2 = plot(xq,yq,'-','linewidth',1);
    set(h2,'Color',colors(i,:))
    end
    xlim([0 50])
    pause(0.5)

end
ylabel('p21')    
ylim([0 1.2])   

end
