function LRI_A_Index = LRIA(input1,MaxK)

[x y]=size(input1);
L_d_1=zeros(x,y);
L_d_2=zeros(x,y);
L_d_3=zeros(x,y);
L_d_4=zeros(x,y);
L_d_5=zeros(x,y);
L_d_6=zeros(x,y);
L_d_7=zeros(x,y);
L_d_8=zeros(x,y);

Threshold=std(input1(:))/2;
for i=MaxK+2:x-(MaxK+2),
    for j=MaxK+2:y-(MaxK+2),
        % DIRECTION = LEFT TOP DIAGONAL
   
       in=input1(i-1:i+1,j-1:j+1);
       curr=in(5);
       diff=abs(curr-in(:))-Threshold;
       diff(find(diff<0))=0;
   % Condition 1
   if(abs(input1(i,j)-input1(i-1,j-1))<Threshold)
         L_d_1(i,j)=0;  
   end
   
   % Condition 2
   R=i-1:-1:i-(MaxK+1);
   C=j-1:-1:j-(MaxK+1);
   index=(C-1).*x + R;
   a=input1(index);
%    a=[input(i-1,j-1) input(i-2,j-2) input(i-3,j-3) input(i-4,j-4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
%            if(k~=1)
           L_d_1(i,j)=k-1;
%            end
           break;
       end
   end

   if(k==MaxK+1)
       L_d_1(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
           L_d_1(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_1(i,j)=-(k-1);
   end
   
   
   % DIRECTION = UP
        
   % Condition 1
   if(abs(input1(i,j)-input1(i-1,j))<Threshold)
         L_d_2(i,j)=0;  
   end
   
   % Condition 2
   
   R=i-1:-1:i-(MaxK+1);
   C=ones(1,MaxK+1)*j;
   index=(C-1).*x + R;
   a=input1(index);
    
   
%    a=[input(i-1,j) input(i-2,j) input(i-3,j) input(i-4,j)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_2(i,j)=k-1;
           end
           break;
       end
   end
   if(k==MaxK+1)
       L_d_2(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_2(i,j)=-(k-1);
           end
           break;
       end
   end
   if(k==MaxK+1)
       L_d_2(i,j)=-(k-1);
   end
   
   
   % DIRECTION = RIGHT TOP DIAGONAL
        
   % Condition 1
   if(abs(input1(i,j)-input1(i-1,j+1))<Threshold)
         L_d_3(i,j)=0;  
   end
   
   % Condition 2
   
   R=i-1:-1:i-(MaxK+1);
   C=j+1:1:j+(MaxK+1);
   ind = sub2ind(size(input1), R, C);
   a=input1(ind);
   
%    a=[input(i-1,j+1) input(i-2,j+2) input(i-3,j+3) input(i-4,j+4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_3(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_3(i,j)=k-1;
   end
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_3(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_3(i,j)=-(k-1);
   end
   
   % DIRECTION = RIGHT
        
   % Condition 1
   if(abs(input1(i,j)-input1(i,j+1))<Threshold)
         L_d_4(i,j)=0;  
   end
   
   % Condition 2
   R=ones(1,(MaxK+1))*i;
   C=j+1:1:j+(MaxK+1);
    index=(C-1).*x + R;
   a=input1(index);
   
%    a=[input(i,j+1) input(i,j+2) input(i,j+3) input(i,j+4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_4(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_4(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_4(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_4(i,j)=-(k-1);
   end
   
   
   % DIRECTION = BOTTOM RIGHT DIAGONAL
        
   % Condition 1
   if(abs(input1(i,j)-input1(i+1,j+1))<Threshold)
         L_d_5(i,j)=0;  
   end
   
   % Condition 2
   
   R=i+1:1:i+(MaxK+1);
   C=j+1:1:j+(MaxK+1);
    
   index=(C-1).*x + R;
   a=input1(index);
   
%    a=[input(i+1,j+1) input(i+2,j+2) input(i+3,j+3) input(i+4,j+4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_5(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_5(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_5(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_5(i,j)=-(k-1);
   end
   
   % DIRECTION = DOWN
        
   % Condition 1
   if(abs(input1(i,j)-input1(i+1,j))<Threshold)
         L_d_6(i,j)=0;  
   end
   
   % Condition 2
   
   R=i+1:1:i+(MaxK+1);
   C=ones(1,(MaxK+1))*j;
    index=(C-1).*x + R;
   a=input1(index);
   
%    a=[input(i+1,j) input(i+2,j) input(i+3,j) input(i+4,j)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_6(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_6(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_6(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_6(i,j)=-(k-1);
   end
   
   
   % DIRECTION = BOTTOM LEFT DIAGONAL
        
   % Condition 1
   if(abs(input1(i,j)-input1(i+1,j-1))<Threshold)
         L_d_7(i,j)=0;  
   end
   
   % Condition 2
   R=i+1:1:i+(MaxK+1);
   C=j-1:-1:j-(MaxK+1);
   index=(C-1).*x + R;
   a=input1(index);
   
%    a=[input(i+1,j-1) input(i+2,j-2) input(i+3,j-3) input(i+4,j-4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_7(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_7(i,j)=k-1;
   end
   
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_7(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_7(i,j)=-(k-1);
   end
   
   % DIRECTION = LEFT
        
   % Condition 1
   if(abs(input1(i,j)-input1(i,j-1))<Threshold)
         L_d_8(i,j)=0;  
   end
   
   % Condition 2
   
   R=ones(1,(MaxK+1))*i;
   C=j-1:-1:j-(MaxK+1);
   index=(C-1).*x + R;
   a=input1(index);
   
%    a=[input(i,j-1) input(i,j-2) input(i,j-3) input(i,j-4)];
   
   for k=1:MaxK+1,
       if(a(k)<=input1(i,j)+Threshold)
           if(k~=1)
               L_d_8(i,j)=k-1;
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_8(i,j)=k-1;
   end
   % Condition 3
   for k=1:MaxK+1,
       if(a(k)>=input1(i,j)-Threshold)
           if(k~=1)
               L_d_8(i,j)=-(k-1);
           end
           break;
       end
   end
   
   if(k==MaxK+1)
       L_d_8(i,j)=-(k-1);
   end
   
    end
end

% subplot(2,4,1)
% imshow(L_d_1);
% subplot(2,4,2)
% imshow(L_d_2);
% subplot(2,4,3)
% imshow(L_d_3);
% subplot(2,4,4)
% imshow(L_d_4);
% subplot(2,4,5)
% imshow(L_d_5);
% subplot(2,4,6)
% imshow(L_d_6);
% subplot(2,4,7)
% imshow(L_d_7);
% subplot(2,4,8)
% imshow(L_d_8);

H_1=hist(L_d_1(:),MaxK*2+1);
H_2=hist(L_d_2(:),MaxK*2+1);
H_3=hist(L_d_3(:),MaxK*2+1);
H_4=hist(L_d_4(:),MaxK*2+1);
H_5=hist(L_d_5(:),MaxK*2+1);
H_6=hist(L_d_6(:),MaxK*2+1);
H_7=hist(L_d_7(:),MaxK*2+1);
H_8=hist(L_d_8(:),MaxK*2+1);

LRI_A_Index=[H_1 H_2 H_3 H_4 H_5 H_6 H_7 H_8];

