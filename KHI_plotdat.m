% load in .mat in the /EECS587/
% myFolder = 'C:\Users\sophmaca\Desktop\Fall 2018\EECS587\timeseries_KHI\';
% load matrix3d.mat

%% watch movie of KHI
for i=1:100 %lenght of t
   data_slice = reshape(matrix3d(:,:,i),[128,128]);
   contourf(data_slice')
   colormap(jet)
   M(i) = getframe;
end

figure
movie(M,1)

%% Make gif of KHI
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimatedKHI.gif';
for n = 1:100
    data_slice = reshape(matrix3d(:,:,n),[128,128]);
    contourf(data_slice','edgecolor','none')  
    colormap(jet)
    title('Kelvin-Helmholtz Instability')
    drawnow 
    
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif','Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
      end 
  end