function [pts3D] =algebraicTriangulation(x1,x2,p1,p2)
% Loading img Frame Points from 1st and 5th camera.
%load('image15.mat')
%img1_x = (image15(1,1,:)-50);
%mg1_y = (image15(1,2,:)-100);
%img2_x = (image15(2,1,:)-50);
%img2_y = (image15(2,2,:)-100);
img1_x = x1(1,:);
img1_y = x1(2,:);
img2_x = x2(1,:);
img2_y = x2(2,:);
% Generating Projection Matrix for the 1st and 5th frame respectively.
X = [];
Y = [];
Z = [];
L = [];
% Calculating the 3D coordinates using Triangulation.
for i=1:562
    i1_x = img1_x(i);
    i1_y = img1_y(i);
    i2_x = img2_x(i);
    i2_y = img2_y(i);
    A = [i1_y*p1(3,1)-p1(2,1),i1_y*p1(3,2)-p1(2,2),i1_y*p1(3,3)-p1(2,3),i1_y*p1(3,4)-p1(2,4);
        i1_x*p1(3,1)-p1(1,1),i1_x*p1(3,2)-p1(1,2),i1_x*p1(3,3)-p1(1,3),i1_x*p1(3,4)-p1(1,4);
        i2_y*p2(3,1)-p2(2,1),i2_y*p2(3,2)-p2(2,2),i2_y*p2(3,3)-p2(2,3),i2_y*p2(3,4)-p2(2,4);
        i2_x*p2(3,1)-p2(1,1),i2_x*p2(3,2)-p2(1,2),i2_x*p2(3,3)-p2(1,3),i2_x*p2(3,4)-p2(1,4)];
    [U,S,V] = svd(A);
    Coordinates_3D = V(1:4,end);
   % Coordinates_3D = Coordinates_3D / Coordinates_3D(4,1);
    X(i) = Coordinates_3D(1,1);
    Y(i) = Coordinates_3D(2,1);
    Z(i) = Coordinates_3D(3,1);
    L(i) = Coordinates_3D(4,1);
    
end
pts3D = [X;Y;Z;L];
%disp(pts3D);
% Display 3D coordinates.
figure;

scatter3(X,Y,Z,'filled');
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

title('Display 3D world points ');
