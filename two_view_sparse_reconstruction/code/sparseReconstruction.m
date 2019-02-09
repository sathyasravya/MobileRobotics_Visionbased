clc;close all;
K = [  558.7087 ,   0.0000,  310.3210 ;
     0.0000,  558.2827 , 240.2395 ;
     0.0000,    0.0000   , 1.0000 ] ;

I1 = rgb2gray(imread('C:\Users\sathya sravya\Documents\MATLAB\mr\Assignment_3\Assignment_3\20161121\images\img1.png'));
I2 = rgb2gray(imread('C:\Users\sathya sravya\Documents\MATLAB\mr\Assignment_3\Assignment_3\20161121\images\img2.png'));
%points1 = detectSURFFeatures(I1);
%points2 = detectSURFFeatures(I2);
points1 = detectSURFFeatures(I1,'MetricThreshold',10);
points2 = detectSURFFeatures(I2,'MetricThreshold',10);
[f1,vpts1] = extractFeatures(I1,points1);
[f2,vpts2] = extractFeatures(I2,points2);
indexPairs = matchFeatures(f1,f2) ;
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));
figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
legend('matched points 1','matched points 2');
%disp(matchedPoints1.Location);
%disp('---------------------');
%disp(matchedPoints2.Location);
a= matchedPoints1.Location;b=matchedPoints2.Location;
mu1 = mean(a,'double');mu2=mean(b,'double');
mu = [mu1;mu2];
%mu = mean(use,'double');
%disp(mu);
d1=0;d2=0;
for i=1:length(a)
    d1 = d1+ sqrt(power(a(i,1),2)+power(a(i,2),2));
    d2 = d2+ sqrt(power(b(i,1),2)+power(b(i,2),2));
end
d1 = d1/length(a);d2 = d2/length(b);
%disp(d1);disp(d2);
%disp(d);
scale1 = sqrt(2)/d1; 
%disp(scale1);
scale2 = sqrt(2)/d2;% disp(scale2);
T1 = [scale1,0,-scale1*mu1(1);0,scale1,-scale1*mu1(2);0,0,1];
%disp(T1);
T2 = [scale2,0,-scale2*mu2(1);0,scale2,-scale2*mu2(2);0,0,1];
%disp(T2);
ap=ones(1,length(a));x1 =[a.';ap];
ap=ones(1,length(b));x2 =[b.';ap];
%disp('=====================================================================================================================================');

x1 = T1*x1;
%disp(x1);
%disp(x1);
x2 = T2*x2;
%disp(x2);
%[F,in]=ransacfitfundmatrix(x1, x2,0.0001);%,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',1e-4);
%disp(F);imin =1; imax=562;y2 = imin+rand(1,562)*(imax-imin);disp(y2);
tresh = ones(8,1)*(1e-7);
%zeros(30,1);
use11 = zeros(3,8);use22 = zeros(3,8);max_inli=0;
fin_list1 = [];fin_list2 = [];
for i=1:60
cnt_inliers =0;
r = randi(562,8,1,'uint32');
usemat =[];
for j=1:8
usemat=[usemat;x1(1,r(j))*x2(1,r(j)),x1(1,r(j))*x2(2,r(j)),x1(1,r(j)),x1(2,r(j))*x2(1,r(j)),x1(2,r(j))*x2(2,r(j)),x1(2,r(j)),x2(1,r(j)),x2(2,r(j)),1];
%use8_1 = [use8_1;x1(1,r(j)),x1(2,r(j))];use8_2=[use8_2;];
use11(:,j) =x1(:,r(j));use22(:,j) =x2(:,r(j)); 
end
[U,D,V] = svd(usemat);F  = reshape(V(:,9),3,3);
tresh1= usemat*V(:,9);
for k=1:8
    if abs(tresh1(k))<tresh(k)
    cnt_inliers=cnt_inliers+1;
    fin_list1=[fin_list1;r(k)];
    fin_list2=[fin_list2;r(k)];
    end
end
 if cnt_inliers >=max_inli
            tresh = tresh1;max_inli= cnt_inliers; 
            %index = i;mmat = r;mmatm = usemat;
 %           disp('===========================');
            x11 = use11;x22 = use22;
  %          disp(tresh1);
 %           disp('===============================');
            ans_F = F;ans_tresh =tresh1 ;
           % disp(cnt_inliers);disp(max_inli);

 end
%if cnt_inliners ==8
 %   fin_list(1,end)=use11();
end
%disp(F);
%end
%disp(ans_F);
disp('====================Fundamental Matrix=============');
%disp(ans_tresh);
de_F =(T2')*(ans_F)*(T1);
[u d v] = svd(de_F); new_d = diag([d(1,1) d(2,2) , 0]); de_F = u * new_d * v.' ;
disp(de_F);%disp(rank(de_F));
E = (K.')*de_F*(K);
disp('=========Essential Matrix========');
disp(E);
[u,d,v] =svd(E);
new_d =  diag([(d(1,1)+d(2,2))/2, (d(1,1)+d(2,2))/2 , 0]); 
E = u * new_d * v.' ;
%disp(E);disp(x11);disp(x22);
%disp(max_inli); disp(rank(E));
%disp(E);
%disp(rank(E));disp(rank(de_F));
%disp(length(x11(0)));
[R, t] = decomposeEssentialMatrix(E, x1, x2, K);%disp(x11);
disp('=========Rotation matrix=======================');%disp(x22);
disp(R);disp('=========Translation matrix=======================');
disp(t);
%disp(fin_list1);
%disp(fin_list2);
figure;
showMatchedFeatures(I1,I2,matchedPoints1(fin_list1,:),matchedPoints2(fin_list2,:),'montage','PlotOptions',{'ro','go','y-'});
title('montage');
TT=horzcat(R,t);TT= [TT;0 0 0 1];
plotCameraFrustum(TT,'k', 1.5);hold on;
