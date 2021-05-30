clear
close all

I = imread('Test1.png');
%I = rgb2gray(I);

[xy,K] = InterfaceCurvature( I , 0.5 , 10 , 10 );

figure;
hold on;
[Nc,~] = size(xy);  % Nc: number of closed contours

for i=1:Nc
    plot(xy{i}(:,1),xy{i}(:,2));
    quiver(xy{i}(:,1),xy{i}(:,2),K{i}(:,1),K{i}(:,2));
end

axis equal;
axis tight;

title('Negative curvature vector field');
xlabel('X');
ylabel('Y');