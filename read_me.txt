1.savedata.m: saves the different channel data separately so that once they are saved they can be used separately decreasing the time to read them for future use.

2.Nucleus_segmentation_tracking.m:segment out the cell nuclei and assign unique number to track them in the latter frames. It outputs Data.mat containing all the tracking indices of the cells. It also creates a video H2.mp4 to visualize the cell tracking.

3.p21_data_extraction.m:extract the p21 average intensity information using the nucleui masks based on the data contained in Data.mat . It saves the p21 data in Result.mat .

4.p21_tracking_plot.m : Generates plot of p21 levels in cells over time.


