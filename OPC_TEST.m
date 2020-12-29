
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
% % plot����(0,0)���½�  �����
%    plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
   patch(boundary(:,1), boundary(:,2),'k');  % ���
end

set(gcf,'unit','normalized','position',[0.1,0.1,0.25,0.454]);  % 1920*1080=16��9 �ֶ��趨���Ϊ����
set(gca,'position',[0 0 1 1]); % axes��figure�е���߽磬�±߽磬��ȣ��߶ȣ���С0�����1

% % ��ȡ��ͼ��Ӧ��XY��Χ
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

M=100;     % resize�߳�ΪM
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
img_target=1-img;       % �ڰ�����

%% ��ͨ�˲�
% ����lybo.m�е�"����ȽϿ��� �������û��ת�Ҷȣ�"
if test_show_im==1
figure,
subplot(2,2,1),imshow(img_target);
title('ԭͼ��');
end

% ����Ƶ�ƶ���ͼ������ģ����Ҳ����Ҫ
s=fftshift(fft2(img_target));

if test_show_im==1
subplot(2,2,3),imshow(log(abs(s)),[]);% imagesc(abs(s));
title('ͼ����Ҷ�任ȡ��������Ƶ��');
end

% ���任���ͼ������ģ����Ǻ����Ĵ���������ͼ���ϵĵ�������ĵ�ľ��������д���
[a,b] = size(img_target);
a0 = round(a/2);
b0 = round(b/2);
d = min(a0,b0)/12;      %12  �˴������������Ķ�Զ��Ƶ�ʲ�Ҫ  
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
title('��ͨ�˲�Ƶ��');
end

img_process = uint8(real(ifft2(ifftshift(low_filter))));

if test_show_im==1
subplot(2,2,2),imshow(img_process,[]);
title('��ͨ�˲����ͼ��');
end

%%
% img_process only 0/1
bw_p=1-logical(img_process);      % �ڰ���ɫ
imwrite(bw_p,save_p_path,'png');

%% show result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('name','result','color','w'),
% show_edge(img_process,'b');
% % ԭͼ������
% show_edge(img_target,'r');
% set(gca,'YDir','reverse');        %��x�᷽������Ϊ����(���ϵ��µ���)��
% axis equal;
% axis off;

%---------------function---------------%
% ��gds2ascii�Ľ���ļ��ж�ȡ���ݵ�cell��
function data_cell=read_data_from(fpath)
    % fpath='./json/2.json';
    fid=fopen(fpath,'r');
    
    %  �����ȡ���ɹ�
    if fid==-1
        fprintf("ERROR: can't read from %s\n",fpath);
        data_cell={};
        return
    end
    
    %  �����ȡ�ɹ�
    xy_find=0;
    data_cell={};
    data_cnt=1;
    xy_array=[];
    xy_array_cnt=1;

    while feof(fid) ~= 1                %�����ж��ļ�ָ��p������ָ���ļ��е�λ�ã�������ļ�ĩ����������1�����򷵻�0  
        line = fgetl(fid);            %��ȡ�ĵ���һ��    

        %    "XY",      [   100,    200,    300,    400     ]
        % ��ʼ��ȡXY
        if contains(line,'"XY"')
            xy_find=1;
            line = fgetl(fid);   % ��������
        end

        if xy_find==1 
            % ��ʽ��ʼ��ȡ
            if contains(line,'[')
                line = fgetl(fid);   % ��������
            end

            % ������ε�XY��ȡ
            if contains(line,']')
                xy_find=0;
                data_cell{data_cnt}=xy_array;
                data_cnt=data_cnt+1;
                xy_array=[];
                xy_array_cnt=1;
            else
                % ���뵽cell
                x=str2double(line(isstrprop(line,'digit')));
                line = fgetl(fid);            %��ȡY
                y=str2double(line(isstrprop(line,'digit')));

                xy_array(xy_array_cnt,:)=[x,y];
                xy_array_cnt=xy_array_cnt+1;
            end     % end of real start getting XY
        end     % end of start getting XY
    end     % end of while

    fclose(fid);
end         % end of function :  read_data_from

%������ȡ
function show_edge(img,color)
    BW = imbinarize(img);
    % imshow(BW);
    [B,L] = bwboundaries(BW);
    hold on
    for k = 1:length(B)
       boundary = B{k};     %��������ľ��ο�
       plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1)    %plot����(0,0)���½�
       skip=10;
       scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1))    %������
    end
end

