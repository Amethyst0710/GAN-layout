[toc]

# 步骤





+ 



# 主要函数

```matlab
function OPC(img_source) % main 
function [img_process,bs,flag,EPE]=opc_process(bs,img_source,img_process_base,EPE_min,knum) % start opc, foreach part pof image do opc
function [img_process,bs,flag,EPE_min]=opc_process_k(bs,img_source,img_process_base,EPE_min,k) % for one part opc, do

function [img,bs]=cal_opc(img_source,bs,type,k) % modified img according to type
function EPE=cal_EPE(img_source,img_process) % calculate EPE

function img_process=filtering(img_target)
function p=show_edge(img,color,is_scatter) 
% following 2 func are used to do filtering in full img
function img=cut_center_img(img_source)
function img_source=restore_center_img(img_source,img_filled)
```



## OPC

<font size=2>**input:** img_source   <font color=green>% 原始未裁剪过的正方形图 </font></font>

主要做开始的前期准备。

先把原图切割，获取图形轮廓坐标（**？这里有个问题是不是很一致**）并存储在元组bs中，调用opc_process开始实行opc，得到[img_process,bs,flag,EPE]。opc_process里对分开的图形依次调用opc_process_k。

##opc_process_k

<font size=2>**input:** bs,img_source,img_process_base,EPE_min,k   </font>

<font size=2 color=green>  % 轮廓相关, origin_not_cut_img, best processed img_cut & EPE, k^th  part graphic</font>

<font size=2>**output:** img_process,bs,flag,EPE_min   </font>

<font size=2 color=green>  % best processed img_cut & EPE,  处理后的轮廓相关,  是否直接达到预期EPE_min</font>

对于第k个图形。首先看是否idx已到0，即所有的轮廓处理完毕，完毕则直接退出。

对于type=0，1，-1分别计算EPE并比较哪个使得EPE小，就选择哪个。

+ 如果是0，直接判断是否小于指定最小值（则记录w, EPE, flag然后退出该循环），或者小于当前最小值，否则继续循环。
+ 如果是1/-1，则对于不同的width判断哪个使得EPE最小并记录
  + 以img_process_base为基础，type,k,w,bs等为输入调用cal_opc_w进行处理。处理完还原到原图中再进行滤波后调用cal_EPE进行计算。

如果此时最小值已经使得超过指定EPE，则记录bs{k,3}中并退出，否则用达到最小值的type和width再进行计算<font size=2 color=gray>（这里可替换成每次最小值记录，?会不会麻烦了呢）</font>得到图片处理结果并以此作为下一次opc_process_k的base，记录bs，然后进行下一轮廓线段的移动<font size=2 color=gray>（轮廓线是由最开始的来确定的，而没有每次进行迭代，？后序可考虑）</font>。















## cal_opc

只根据current_idx和type对现有image进行修改。事实上current_idx还是基于最初图像的划分。

## cal_EPE

感觉有比较大的问题。现在是切割线段的一个端点对应的在滤波后的图形的点到目标的距离。

<img src=".\README.assets\5700B4BD2E566F509FBCA1C856D526A9.png" alt="EPE" style="zoom:50%;" />





<img src=".\README.assets\边缘放置误差.jpg" style="zoom: 80%;" />



# 注意事项

+ OPC只对完整图中的一部分做，但每次filter都是把处理好的的填回去然后进行滤波处理。
+ bs是一个轮廓相关数据的cell，k表示第k个轮廓、图形，
  + bs(k,1)=xs_list,    		bs(k,2)=ys_list, 
  + bs(k,3)=choose_type_list （每条边向内/外移动的量，边缘轮廓不做处理--0）, 
  + bs(k,4)=current_idx，初始值为length。
+ 做cal_opc时候，事实上current_idx还是基于最初图像的划分，而不是根据新图像重新进行线段划分，这点要注意。
+ 所有的xy都是第x行，第y列的意思，因为矩阵存储时img(x,y)表示试试这样的。

存在问题

√ 目前的错是因为cal_opc传入进来的图像不是初始图像，是处理后的，所以judge corner使用对象错误。

⭕ 错位，差一格没操作？,实际画的时候应该画一个点以此避免重叠的问题，现在的做法可能会导致EPE时候算错吗？

❓ 一直都是最后一个w的值的问题？代码错误，已更正。

❓仍然有5：-1：2与2：5结果不一样的问题？是因为出现两个EPE相等情况时，是选取先计算得到最小值的，所以两个顺序造成结果不同。

❓左下角没受边缘控制不画图，因为不是直线type，TODO？