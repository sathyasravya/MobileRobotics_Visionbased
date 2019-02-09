S = load('C:\Users\sathya sravya\Desktop\MR_20161121_Assignment_4\MatFilesQues1\cube_imgs.mat');
%disp(S);
S2= load('C:\Users\sathya sravya\Desktop\MR_20161121_Assignment_4\MatFilesQues1\projMatrices.mat');
for j=1:56 %for different images 
  A = [];
for i=1:8 %takes care of number of projection matrices under consideration
pts2D_view3 = squeeze(S.image_pts(i,:,j));
use_x = pts2D_view3(1);
use_y = pts2D_view3(2);
%disp(use_x);disp(use_y);
projMat_view3 =  S2.projMatrices{i};  %each projection matrix
%disp(projMat_view3);
A1 = [use_y*projMat_view3(3,1)-projMat_view3(2,1) use_y*projMat_view3(3,2)-projMat_view3(2,2) use_y*projMat_view3(3,3)-projMat_view3(2,3) use_y*projMat_view3(3,4)-projMat_view3(2,4)];
A2 = [use_x*projMat_view3(3,1)-projMat_view3(1,1) use_x*projMat_view3(3,2)-projMat_view3(1,2) use_x*projMat_view3(3,3)-projMat_view3(1,3) use_x*projMat_view3(3,4)-projMat_view3(1,4)];
A=vertcat(A,A1,A2);
end
%disp(A);
[U,S11,V] = svd(A);
%disp(V);
    Coordinates_3D = V(1:4,end);
   % disp(Coordinates_3D);
    Coordinates_3D = Coordinates_3D / Coordinates_3D(4,1);
    X(j) = Coordinates_3D(1,1);
    Y(j) = Coordinates_3D(2,1);
    Z(j) = Coordinates_3D(3,1);
    L(j) = Coordinates_3D(4,1);
    %plo(X(j),Y(j),Z(j)); hold on;
end
figure;
pts3D = [X;Y;Z;L];
scatter3(X,Y,Z,'filled');
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
title('Reconstructed cube from all views given - 3D world points ');



