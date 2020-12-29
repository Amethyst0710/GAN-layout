
%%  read data from json
% clc;
% clear all;

% % jsondir='./json/';
% fname='11.json';
% fpath=strcat(jsondir,fname);

data_cell=read_data_from(fpath);
if isempty(data_cell)
    return
end

fig=figure;
if test_show_im==0
    set(fig, 'visible', 'off');
end

hold on
for k = 1:length(data_cell)
   boundary = data_cell{k};    
% % plot坐标(0,0)左下角  仅描边
%    plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
   patch(boundary(:,1), boundary(:,2),'k');  % 填充
end

set(gcf,'unit','normalized','position',[0.1,0.1,0.25,0.454]);  % 1920*1080=16：9 手动设定输出为方形
set(gca,'position',[0 0 1 1]); % axes在figure中的左边界，下边界，宽度，高度，最小0，最大1

% % 截取版图对应的XY范围
% xx=1100;
% yy=1200;
% radius=1500;
axis([xx-radius,xx+radius,yy-radius,yy+radius]);
axis off;

% save_o_path=strcat(picdir,fname(1:end-5),'.png');
saveas(gcf, save_o_path);
% saveas(gcf, 'save.jpg');
close(fig);

%%   resize
img=imread(save_o_path);
% img=imread( 'save.jpg');
% BW = im2bw(img,0.5);
BW = imbinarize(rgb2gray(img),0.5);

M=100;     % resize边长为M
% fig=figure,
% imshow(BW)
im = imresize(BW,[M M]); % Resize the picture for alexnet
% imshow(im)

% save_r_path=strcat(resizedir,fname(1:end-5),'.png');
imwrite(im,save_r_path,'png')
% imwrite(im,'save_res.png','png')
% close(fig);


%% read from pic
img=imread(save_r_path);
img_target=1-img;       % 黑白问题

%% 低通滤波
% 调用lybo.m中的"好像比较靠谱 但是这个没有转灰度，"
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

