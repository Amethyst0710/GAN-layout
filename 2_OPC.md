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

主要做开始的前期准备。

先把原图切割，获取图形边缘坐标（？）并存储在元组bs中，调用opc_process开始实行opc。

##opc_process_k

对于第k个图形，

（？）没有用迭代，

## cal_opc



## cal_EPE

感觉有比较大的问题



# 注意事项

+ OPC只对完整图中的一部分做，但每次filter都是把处理好的的填回去然后进行滤波处理。
+ 

