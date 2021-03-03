%%%% last copied: 01/29 
%%%% this file works on windows for process sth.
%% not from json but draw myself for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start & global setting
clc;
clear all;
global test_show_im
test_show_im=1;
global xpcnt;global ypcnt;
xpcnt=0.33;ypcnt=0.33;
global skip;global opc_width;
skip=10;     %采样点间隔
opc_width=5;
global minEPE;
minEPE=0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%?
%% original pic draw
w=200;h=200;
img_target=zeros(w,h);
% draw rec L-1
w1=20;h1=50;
x1=70;y1=70;
x2=x1+w1;y2=y1+h1;
img_target=draw_rec(img_target,x1,y1,x2,y2,1);  %%%%%%%
% draw rec L-2
w2=50;h2=20;
x3=x2;y4=y2;
x4=x3+w2;y3=y4-h2;
x2=70;y2=60; %L横线的右上角
img_target=draw_rec(img_target,x3,y3,x4,y4,1);

img_target=cloned_part_img(img_target,140,20,img_target,x1,y1,40,40);
img_target=cloned_part_img(img_target,20,140,img_target,x1,y1,40,60);
img_target=cloned_part_img(img_target,1,1,img_target,x1,y1,20,70);
img_target=cloned_part_img(img_target,160,150,img_target,x1,y1,40,40);
% % OPC ear
% r=0;
% img_target=ear4(img_target,x1,y1,x3,y3,r,1);

% img_target=draw_rec(img_target,1,1,20,20,1);
% % 靠近边缘环形测试
% xxx=5;yyy=95;www=10;
% img_target=draw_rec(img_target,xxx,xxx,xxx+www,yyy,1);         %%%%%%%%%%%%%%%
% img_target=draw_rec(img_target,xxx,xxx,yyy,xxx+www,1);
% img_target=draw_rec(img_target,yyy-www,xxx,yyy,yyy,1);        
% img_target=draw_rec(img_target,xxx,yyy-www,yyy,yyy,1);

% % L形凹进去
% r=3;
% img_target=draw_rec(img_target,x3-r,y2-r,x3+r,y2+r,0);

img_target_wb=1-img_target;        
% img是黑0白1，现取反，为看起来方便画图区域应该是黑
% ？取消这步之后没有边框问题了,
% 所以实际处理时应该白1表示画了的掩膜图形,即处理过程中显示的是黑底
imshow(img_target_wb);
show_edge(img_target,'r');
%% 低通滤波
img_process=filtering(img_target);
% figure,
% imshow(img_process,[]);
%% cut center pic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_target_cut=cut_center_img(img_target);
img_process_cut=cut_center_img(img_process);
if test_show_im==1
    figure()
    subplot(1,2,1),imshow(img_target_cut);
    show_edge(img_target_cut,'r');
    subplot(1,2,2),imshow(img_process_cut,[]);
    show_edge(img_process_cut,'b');
end
%% OPC
OPC(img_target_cut,img_process_cut);

%%
% img_process only 0/1
bw_p=1-logical(img_process);      % 黑白显色
imwrite(bw_p,save_p_path,'png');

%% show result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('name','result','color','w'),
% show_edge(img_process,'b');
% % 原图形轮廓
% show_edge(img_target,'r');
% set(gca,'YDir','reverse');        %将x轴方向设置为反向(从上到下递增)。
% axis equal;
% axis off;

%%
%---------------function---------------%
% 从gds2ascii的结果文件中读取数据到cell中
function data_cell=read_data_from(fpath)
    % fpath='./json/2.json';
    fid=fopen(fpath,'r');
    
    %  如果读取不成功
    if fid==-1
        fprintf("ERROR: can't read from %s\n",fpath);
        data_cell={};
        return
    end
    
    %  如果读取成功
    xy_find=0;
    data_cell={};
    data_cnt=1;
    xy_array=[];
    xy_array_cnt=1;

    while feof(fid) ~= 1                %用于判断文件指针p在其所指的文件中的位置，如果到文件末，函数返回1，否则返回0  
        line = fgetl(fid);            %获取文档第一行    

        %    "XY",      [   100,    200,    300,    400     ]
        % 开始读取XY
        if contains(line,'"XY"')
            xy_find=1;
            line = fgetl(fid);   % 跳过该行
        end

        if xy_find==1 
            % 正式开始读取
            if contains(line,'[')
                line = fgetl(fid);   % 跳过该行
            end

            % 结束这次到XY读取
            if contains(line,']')
                xy_find=0;
                data_cell{data_cnt}=xy_array;
                data_cnt=data_cnt+1;
                xy_array=[];
                xy_array_cnt=1;
            else
                % 读入到cell
                x=str2double(line(isstrprop(line,'digit')));
                line = fgetl(fid);            %获取Y
                y=str2double(line(isstrprop(line,'digit')));

                xy_array(xy_array_cnt,:)=[x,y];
                xy_array_cnt=xy_array_cnt+1;
            end     % end of real start getting XY
        end     % end of start getting XY
    end     % end of while

    fclose(fid);
end         % end of function :  read_data_from

%轮廓提取
function show_edge(img,color)
    BW = imbinarize(img);
    % imshow(BW);
    [B,L] = bwboundaries(BW);
    hold on
    for k = 1:length(B)
       boundary = B{k};     %包含最外的矩形框
       plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1)    %plot坐标(0,0)左下角
       skip=10;
       scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1))    %采样点
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function img=draw_rec(img,x1,y1,x3,y3,z)
% (i,j)矩阵第i行第j列，等价于MATLAB坐标位置(j,MAMY-i)，设置左上角为(0,0)则(j,i)
% 现在必须输入坐左上右下，以后可以优化 ###TODO
function img=draw_rec(img,y1,x1,y3,x3,z)
    % [a,b] = size(img);
    for i=x1:1:x3
        for j=y1:1:y3
            img(i,j)=z;
        end
    end
end

function img=ear4(img,x1,y1,x3,y3,r,z)
    img=draw_rec(img,x1-r,y1-r,x1+r,y1+r,z);
    img=draw_rec(img,x3-r,y1-r,x3+r,y1+r,z);
    img=draw_rec(img,x1-r,y3-r,x1+r,y3+r,z);
    img=draw_rec(img,x3-r,y3-r,x3+r,y3+r,z);
end

% clone part from img_source(xs,ys) to img_dest(xd,yd) 
% 现在不能超出矩阵范围，以后可优化只粘贴重叠部分 ###TODO
function img=cloned_part_img(img_dest,xd,yd,img_source,xs,ys,w,h)
    img=img_dest;
    img(xd:xd+w-1,yd:yd+h-1)=img_source(xs:xs+w-1,ys:ys+h-1);
end

function img=cut_center_img(img_source)
    global xpcnt;global ypcnt;
    [a,b]=size(img_source);
    w=floor(a*xpcnt);
    h=floor(b*ypcnt);
    x1=floor((a-w)/2);
    y1=floor((b-h)/2);
    img=img_source(x1:x1+w,y1:y1+h);
end

function img_source=restore_center_img(img_source,img_filled)
    global xpcnt;global ypcnt;
    [a,b]=size(img_source);
    w=floor(a*xpcnt);
    h=floor(b*ypcnt);
    x1=floor((a-w)/2);
    y1=floor((b-h)/2);
    img_source(x1:x1+w,y1:y1+h)=img_filled;
end

function img_process=filtering(img_target)
    % 低通滤波
    % 调用lybo.m中的"好像比较靠谱 但是这个没有转灰度，"
    global test_show_im
    if test_show_im==1
        figure,
        subplot(2,2,1),imshow(img_target);
        title('原图像');
    end

    % 将低频移动到图像的中心，这个也很重要
    s=fftshift(fft2(img_target));

    if test_show_im==1
        subplot(2,2,3),imshow(log(abs(s)),[]);% imagesc(abs(s));
        title('图像傅里叶变换取对数所得频谱');
    end

    % 求解变换后的图像的中心，我们后续的处理是依据图像上的点距离中心点的距离来进行处理
    [a,b] = size(img_target);
    a0 = round(a/2);
    b0 = round(b/2);
    d = min(a0,b0)/12;      %12  此处决定距离中心多远的频率不要  
    % d = 5;
    d = d^2;
    low_filter=zeros(a,b);
    for i=1:a
        for j=1:b
            distance = (i-a0)^2+(j-b0)^2;
            if distance<d
                low_filter(i,j) = s(i,j);
            end
        end
    end

    if test_show_im==1
        subplot(2,2,4),imshow(log(abs(low_filter)),[]);% imagesc(abs(low_filter));
        title('低通滤波频谱');
    end

    img_process = uint8(real(ifft2(ifftshift(low_filter))));

    if test_show_im==1
        subplot(2,2,2),imshow(img_process,[]);
        title('低通滤波后的图像');
    end
end

function [img,bs]=cal_opc(img_source,bs)

end

function EPE=cal_EPE(img_source,img_process)

end

function [img_process,bs]=opc_process(img_process,bs,img_source,img_source_i,EPE)
    global minEPE
    % 终止条件：达到minEPE or 全部情况测试完
    function flag=types_opc_process(type)
        [img_process,bs]=cal_opc(img_process,bs,type);

        img_process_o=restore_center_img(img_source,img_process);       
        img_filtered_o=filtering(img_process_o);        
        img_filtered=cut_center_img(img_filtered_o);

        EPE=cal_EPE(img_source_i,img_filtered);
        
        if( EPE>minEPE || sum(sum(cell2mat(bs(:,3))==0))==0 )  
            flag=true
            return
        else
            [img_process,bs]=opc_process(img_process,bs,img_source_o,img_source_i,EPE);
        end
        flag=false
    end
    
    for type = [0,1,-1]
        if types_opc_process(type)
            return
        end
    end
    % if all failed TODO
    
    % cal_opc: choose a line, draw new img, 0
    % cal_opc: choose a line, draw new img, 1
    % cal_opc: choose a line, draw new img, -1
end

function OPC(img_source)
    img_process_o=img_source;
    img_source_i=cut_center_img(img_source);
    img_process=img_source_i;   % --> img_opc not filtered
    % _o后缀是原大小图片、未切割
    % img_source_i 原图片切割后结果
    % img_source_o 切割前，滤波前，一直在变化
    
    
    % to get initial work, sample & stack
    global skip;
    BW = imbinarize(img_process);
    [B,L] = bwboundaries(BW);
    bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
    %cell2mat(p(1,1)) to get data
    for k = 1:length(B)
       boundary = B{k};     %包含最外的矩形框
       xs=boundary(1:skip:length(boundary),2);
       ys=boundary(1:skip:length(boundary),1);
       bs(k,1)={xs};
       bs(k,2)={ys};
       bs(k,3)={zeros(1,len(xs))};
       bs(k,4)={len(xs)};
    end
     
    img_filtered_o=filtering(img_process_o);  % filtering need no cut image    
    img_filtered=cut_center_img(img_filtered_o);
    
    EPE=cal_EPE(img_source_i,img_filtered);
    
    %%%% cal_min_EPE and its img_process TODO
    opc_process(img_process,bs,img_source,img_source_i,EPE);  % get: img_process EPE
    
%     global minEPE
%     % 终止条件：达到minEPE or 全部情况测试完
%     while( EPE>minEPE || sum(sum(cell2mat(bs(:,3))==0))==0 )  
%         [img_process,bs]=cal_opc(img_process,bs);
%         
%         img_process_o=restore_center_img(img_source,img_process);       
%         img_filtered_o=filtering(img_process_o);        
%         img_filtered=cut_center_img(img_filtered_o);
% 
%         EPE=cal_EPE(img_source_i,img_filtered);
%     end
    
    % end
    img_process_o=restore_center_img(img_source,img_process);       
    img_filtered_o=filtering(img_process_o);        
    img_filtered=cut_center_img(img_filtered_o);
    
    % show result
    figure();
    subplot(1,2,1),imshow(img_source_i);
    show_edge(img_source_i,'r');    
    subplot(1,3,2),imshow(img_process);
    show_edge(img_process,'r');
    subplot(1,3,3),imshow(img_filtered,[]);
    show_edge(img_filtered,'r');
    EPE
    
end