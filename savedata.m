function savedata
% To save the channel data separately 
 
 
 file_name = 'C:\Users\Debojyoti\Desktop\New folder\RPEp21GFP_movies\mpspt548hrs then meduia change add h2 start movie .nd2 - mpspt548hrs then meduia change add h2 start movie .nd2 (series 01)';
 Path_to_Save = 'Video_Result_p21_RPE'; 

  index_backslash = strfind(file_name,'\');
  name = file_name(index_backslash(end)+1:end);

  
  Path_to_Save_complete = [Path_to_Save,'\',name];
  if exist('fname','dir') == 0
     mkdir(Path_to_Save_complete);
  end


     tiffile = strcat(file_name,'.tif');
     info = imfinfo(tiffile);
     nimages = length(info);
     cols = info(1).Width;
     rows = info(1).Height;
     I_tmp = uint16(zeros(rows,cols,nimages));
     
     
     %% Saving  channel data   separately

      for i = 1:1:nimages
         tmp = imread(tiffile,i);    
         I_tmp(:,:,i) = tmp;
      end

 DIC_image = I_tmp(:,:,1:4:end);
 PCNA_image = I_tmp(:,:,2:4:end);
 H2_image = I_tmp(:,:,3:4:end);
 P21_image = I_tmp(:,:,4:4:end);

 save([Path_to_Save_complete, '\DIC_image.mat'],'DIC_image')
 save([Path_to_Save_complete, '\PCNA_image.mat'],'PCNA_image')
 save([Path_to_Save_complete, '\H2_image.mat'],'H2_image')
 save([Path_to_Save_complete, '\P21_image.mat'],'P21_image')
 
end
