# this file for create coordinate,
# input x1,y1,x2,y2,interval then create coordinate /
x1 = 2696
y1 = 3504
x2 = 2747
y2 = 3555
### x1~x2 y1~y2  ～50*50
# need to cut bigger, for in matlab, it will cut again
## radius = 5
target_x = 2  # 50 # traget number only int
target_y = target_x
# interval = 20 #2
# num=int(target/interval)
with open('coor', 'w', encoding='utf-8') as f:
    f.write("X  Y\n")
    for x in [x1+i*(x2-x1)/target_x for i in range(target_x)]:
        for y in [y1 + j * (y2 - y1) / target_y for j in range(target_y)]:
            f.write(str(x) + " " + str(y) + "\n")
