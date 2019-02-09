%SathyaSravya.Vallabhajyosyula 
%-------------------------------------------------------
close all;clc;
% Define images to process 
h =zeros(3,3);H=zeros(9,9);e = zeros(3,1);
Zz =zeros(3,3);Yy =zeros(3,3);Xx =zeros(3,3);
v11=zeros(3,6); v12=zeros(3,6);v22=zeros(3,6);  
imageFileNames = {'C:\Users\sathya sravya\Documents\MATLAB\mr\images\img1.png',...
    'C:\Users\sathya sravya\Documents\MATLAB\mr\images\img2.png',...
    'C:\Users\sathya sravya\Documents\MATLAB\mr\images\img3.png',...
    'C:\Users\sathya sravya\Documents\MATLAB\mr\images\img4.png',...
    'C:\Users\sathya sravya\Documents\MATLAB\mr\images\img5.png',...
    };
% Detect checkerboards in images
[imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints(imageFileNames);
imageFileNames = imageFileNames(imagesUsed);
% Read the first image to obtain image size
originalImage = imread(imageFileNames{1});
[mrows, ncols, ~] = size(originalImage);
% Generate world coordinates of the corners of the squares
squareSize = 24;  % in units of 'millimeters'
for i = 1:numel(imageFileNames)
  I = imread(imageFileNames{i});
  figure
  %subplot(5, 1, i);
  imshow(I);
  hold on;
  plot(imagePoints(:,1,i),imagePoints(:,2,i),'ro');
end
worldPoints = generateCheckerboardPoints(boardSize, squareSize);
%display(imagePoints);

%ax^T = (-Xi,-Yi,-Zi,-1,0,0,0,0,xiXi,xiYi,xiZi,xi);
%ay^T = (0,0,0,0,-Xi,-Yi,-Zi,-1,yiXi,yiYi,yiZi,yi);
%for i=1:1:4
X(1) =worldPoints(1,1);Y(1) = worldPoints(1,2);Z(1,1) =0;%worldPoints() ;
   X(2) = worldPoints(6,1);Y(2)= worldPoints(6,2);Z(2)=0;
   X(3) = worldPoints(24,1);Y(3)=worldPoints(24,2);Z(3)=0;
   X(4) = worldPoints(48,1);Y(4)=worldPoints(48,2);Z(4)=0;
for imnum=1:1:3
axT=zeros(4,9);ayT=zeros(4,9);x=zeros(4,1);y=zeros(4,1);
%X=zeros(4,1);Y=zeros(4,1);Z=zeros(4,1);
%end
%if(i==1)
x(1)=imagePoints(1+48*(1-1),1,imnum);
y(1)=imagePoints(1+48*(1-1),2,imnum);
%end
%if(i==2)
x(2)=imagePoints(6+48*(1-1),1,imnum);
y(2)=imagePoints(6+48*(1-1),2,imnum);
%end
%if(i==3)
x(3)=imagePoints(24+48*(1-1),1,imnum);
y(3)=imagePoints(24+48*(1-1),2,imnum);
%end
%if(i==4)
x(4)=imagePoints(48+48*(1-1),1,imnum);
y(4)=imagePoints(48+48*(1-1),2,imnum);
%end
for i=1:1:4
axT(i,:) =[-X(i) -Y(i) -1 0 0 0 x(i)*X(i) x(i)*Y(i) x(i)];
ayT(i,:) =[0 0 0 -X(i) -Y(i) -1 y(i)*X(i) y(i)*Y(i) y(i)];
end
a = [axT(1,:);ayT(1,:);axT(2,:);ayT(2,:);axT(3,:);ayT(3,:);axT(4,:);ayT(4,:)];
%disp('mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm');
%disp(a);

% See additional examples of how to use the calibration data.  At the prompt type:
% showdemo('MeasuringPlanarObjectsExample')
% showdemo('StructureFromMotionExample')
[U,S,V] = svd(a);
%disp(U);disp(S);
%disp("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
%disp(V);
%disp("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");

%homography
h= [V(1,9) V(2,9) V(3,9) ;V(4,9) V(5,9) V(6,9);V(7,9) V(8,9) V(9,9)];
disp('hhhhhhhhhhh');
disp(h);
if(imnum==1)Xx=h; end
if(imnum==2)Yy=h; end
if(imnum==3)Zz=h; end
%H(3*(k-1)+1:3*(k-1)+3,3*(k-1)+1:3*(k-1)+3)=h;
v12(imnum,:)=[h(1,1)*h(1,2) h(1,1)*h(2,2)+h(2,1)*h(1,2) h(3,1)*h(1,2)+h(1,1)*h(3,2) h(2,1)*h(2,2) h(3,1)*h(2,2)+h(2,1)*h(3,2) h(3,1)*h(3,2)];
v11(imnum,:)=[h(1,1)*h(1,1) h(1,1)*h(2,1)+h(2,1)*h(1,1) h(3,1)*h(1,1)+h(1,1)*h(3,1) h(2,1)*h(2,1) h(3,1)*h(2,1)+h(2,1)*h(3,1) h(3,1)*h(3,1)];
v22(imnum,:)=[h(1,2)*h(1,2) h(1,2)*h(2,2)+h(2,2)*h(1,2) h(3,2)*h(1,2)+h(1,2)*h(3,2) h(2,2)*h(2,2) h(3,2)*h(2,2)+h(2,2)*h(3,2) h(3,2)*h(3,2)];
end
VV = [v12(1,:);v11(1,:)-v22(1,:);v12(2,:);v11(2,:)-v22(2,:);v12(3,:);v11(3,:)-v22(3,:)];%v12(4,:);v11(4,:)-v22(4,:); v12(5,:);v11(5,:)-v22(5,:);];
%disp("yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy");
disp('The matrix "V" to be solved');
disp(VV);
[u1,s1,v1]=svd(VV);
%disp("finaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaal");
%disp(u1);disp(s1);
%disp('');
%disp(v1);
b= [v1(1,6) v1(2,6) v1(3,6) ;v1(2,6) v1(4,6) v1(5,6);v1(3,6) v1(5,6) v1(6,6)];
b = (b+b.')/2;
e=eig(b);
while(min(e)<=0)
    e=eig(b);
    if(min(e)==0)
       b = b +pow(10,-12)*[1 0 0;0 1 0;0 0 1];  
   
    else (min(e)<0)
            b = b + (min(e)*min(e))*[1 0 0;0 1 0;0 0 1];
    end
end
disp('The matrix B on which chol is to be applied');
disp(b);
k1=chol(b);
disp('The K matrix');
display(inv(k1));%This is the final K matrix obtained
disp('Rotation matrices  R1,2|t');
%disp(Xx);
disp(k1*Xx);%this gives Rotation matrices
disp(k1*Yy);
disp(k1*Zz);
%disp();
