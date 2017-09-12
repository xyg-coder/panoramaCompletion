function [ y8 ] = GXYfft(gray_img)
% source_img=imread(image_path);
% gray_img=rgb2gray(gray_img);
gray_img=int64(gray_img);
img_size=size(gray_img);
gray_img2=gray_img(:,2:img_size(2));
y2_1_1=gray_img(:,1:img_size(2)-1)+gray_img2;
y2_1_2=gray_img(:,1:img_size(2)-1)-gray_img2;
y2=cat(3,y2_1_1,y2_1_2);
d8_0_1=y2_1_1(9:img_size(1)-1,:);%用于计算d8(j);
d8_0_2=y2_1_2(9:img_size(1)-1,:);%第二个维度计算d8(j);
d8_1_1=y2_1_1(10:img_size(1),:);
d8_1_2=y2_1_2(10:img_size(1),:);%计算d8(j+1);
d8_j_1=y2_1_1(1:img_size(1)-9,:)-d8_0_1;
d8_j1_1=y2_1_1(2:img_size(1)-8,:)-d8_1_1;
d8_j_2=y2_1_2(1:img_size(1)-9,:)-d8_0_2;
d8_j1_2=y2_1_2(2:img_size(1)-8,:)-d8_1_2;

t2_0{1}=d8_j_1+d8_j1_1;
t2_0{2}=d8_j_2+d8_j1_2;
t2_1{1}=d8_j_1-d8_j1_1;
t2_1{2}=d8_j_2-d8_j1_2;

y8_2=cell(2,8);%为patch尺寸为8*2时的cell;
for i=1:16
    y8_2{i}=int64(zeros(img_size(1)-7,img_size(2)-1));
end;

for i=1:2
    for j=1:2
         y8_2{i,1}(j,:)=y2(j,:,i)+y2(j+1,:,i)+y2(j+2,:,i)+y2(j+3,:,i)+y2(j+4,:,i)+y2(j+5,:,i)+y2(j+6,:,i)+y2(j+7,:,i);
         y8_2{i,2}(j,:)=y2(j,:,i)+y2(j+1,:,i)+y2(j+2,:,i)+y2(j+3,:,i)-y2(j+4,:,i)-y2(j+5,:,i)-y2(j+6,:,i)-y2(j+7,:,i);
         y8_2{i,3}(j,:)=y2(j,:,i)+y2(j+1,:,i)-y2(j+2,:,i)-y2(j+3,:,i)-y2(j+4,:,i)-y2(j+5,:,i)+y2(j+6,:,i)+y2(j+7,:,i);
         y8_2{i,4}(j,:)=y2(j,:,i)+y2(j+1,:,i)-y2(j+2,:,i)-y2(j+3,:,i)+y2(j+4,:,i)+y2(j+5,:,i)-y2(j+6,:,i)-y2(j+7,:,i);
         y8_2{i,5}(j,:)=y2(j,:,i)-y2(j+1,:,i)-y2(j+2,:,i)+y2(j+3,:,i)+y2(j+4,:,i)-y2(j+5,:,i)-y2(j+6,:,i)+y2(j+7,:,i);
         y8_2{i,6}(j,:)=y2(j,:,i)-y2(j+1,:,i)-y2(j+2,:,i)+y2(j+3,:,i)-y2(j+4,:,i)+y2(j+5,:,i)+y2(j+6,:,i)-y2(j+7,:,i);
         y8_2{i,7}(j,:)=y2(j,:,i)-y2(j+1,:,i)+y2(j+2,:,i)-y2(j+3,:,i)-y2(j+4,:,i)+y2(j+5,:,i)-y2(j+6,:,i)+y2(j+7,:,i);
         y8_2{i,8}(j,:)=y2(j,:,i)-y2(j+1,:,i)+y2(j+2,:,i)-y2(j+3,:,i)+y2(j+4,:,i)-y2(j+5,:,i)+y2(j+6,:,i)-y2(j+7,:,i); 
    end;
end;
for i=1:2
    for row=1:img_size(1)-9
        y8_2{i,1}(row+2,:)=y8_2{i,1}(row,:)-t2_0{i}(row,:);
        y8_2{i,2}(row+2,:)=-y8_2{i,3}(row,:)+t2_0{i}(row,:);
        y8_2{i,3}(row+2,:)=y8_2{i,2}(row,:)-t2_0{i}(row,:);
        y8_2{i,4}(row+2,:)=-y8_2{i,4}(row,:)+t2_0{i}(row,:);
        y8_2{i,5}(row+2,:)=-y8_2{i,5}(row,:)+t2_1{i}(row,:);
        y8_2{i,6}(row+2,:)=y8_2{i,7}(row,:)-t2_1{i}(row,:);
        y8_2{i,7}(row+2,:)=-y8_2{i,6}(row,:)+t2_1{i}(row,:);
        y8_2{i,8}(row+2,:)=y8_2{i,8}(row,:)-t2_1{i}(row,:);
    end;
end;
y8=cell(8,8);
for i=1:64
    y8{i}=zeros(img_size(1)-7,img_size(2)-7);
end;
% for i=1
for j=1:8
   y8{j,1}(:,1:img_size(1)-7)=y8_2{1,j}(:,1:img_size(1)-7)+y8_2{1,j}(:,3:img_size(1)-5)...
        +y8_2{1,j}(:,5:img_size(1)-3)+y8_2{1,j}(:,7:img_size(1)-1);
    y8{j,2}(:,1:img_size(1)-7)=y8_2{1,j}(:,1:img_size(1)-7)+y8_2{1,j}(:,3:img_size(1)-5)...
        -y8_2{1,j}(:,5:img_size(1)-3)-y8_2{1,j}(:,7:img_size(1)-1);
    y8{j,3}(:,1:img_size(1)-7)=y8_2{1,j}(:,1:img_size(1)-7)-y8_2{1,j}(:,3:img_size(1)-5)...
        -y8_2{1,j}(:,5:img_size(1)-3)+y8_2{1,j}(:,7:img_size(1)-1);
    y8{j,4}(:,1:img_size(1)-7)=y8_2{1,j}(:,1:img_size(1)-7)-y8_2{1,j}(:,3:img_size(1)-5)...
        +y8_2{1,j}(:,5:img_size(1)-3)-y8_2{1,j}(:,7:img_size(1)-1);
    y8{j,5}(:,1:img_size(1)-7)=y8_2{2,j}(:,1:img_size(1)-7)-y8_2{2,j}(:,3:img_size(1)-5)...
        +y8_2{2,j}(:,5:img_size(1)-3)-y8_2{2,j}(:,7:img_size(1)-1);
    y8{j,6}(:,1:img_size(1)-7)=y8_2{2,j}(:,1:img_size(1)-7)-y8_2{2,j}(:,3:img_size(1)-5)...
        -y8_2{2,j}(:,5:img_size(1)-3)+y8_2{2,j}(:,7:img_size(1)-1);
    y8{j,7}(:,1:img_size(1)-7)=y8_2{2,j}(:,1:img_size(1)-7)+y8_2{2,j}(:,3:img_size(1)-5)...
        -y8_2{2,j}(:,5:img_size(1)-3)-y8_2{2,j}(:,7:img_size(1)-1);
     y8{j,8}(:,1:img_size(1)-7)=y8_2{2,j}(:,1:img_size(1)-7)+y8_2{2,j}(:,3:img_size(1)-5)...
        +y8_2{2,j}(:,5:img_size(1)-3)+y8_2{2,j}(:,7:img_size(1)-1);
end;

%以下可供检查用
% end;
% m1=0;
% m2=0;
% q=zeros(8,8);
% for i=1:8
%     for j=1:8
%         if(rem(i,8)==1||rem(i,8)==4||rem(i,8)==6||rem(i,8)==7) heng=1;
%         else heng=-1;
%         end;
%         if(rem(j,8)==1||rem(j,8)==2||rem(j,8)==7||rem(j,8)==0) shu=1;
%         else shu=-1;
%         end;
%         q(j,i)=heng*shu;
%     end;
% end;
% q2=zeros(8,8);
% for i=1:8
%     for j=1:8
%         if(rem(i,8)==1||rem(i,8)==4||rem(i,8)==6||rem(i,8)==7) shu=1;
%         else shu=-1;
%         end;
%         if(rem(j,8)==1||rem(j,8)==2||rem(j,8)==7||rem(j,8)==0) heng=1;
%         else heng=-1;
%         end;
%         q2(i,j)=heng*shu;
%     end;
% end;
% q=q2';
% for i=1:8
%     for j=1:8
%         m1=m1+q(i,j)*gray_img(12+i,8+j);
%         m2=m2+q2(i,j)*gray_img(12+i,8+j);
%     end;
% end;