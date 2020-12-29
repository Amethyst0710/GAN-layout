# guments:
# input_file coor_file radius output
import pya
import time

###########
store_format="gdsii"

print("Start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
ly = pya.Layout()
# print("Parse layout file start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
ly.read(input_file)

# print("Parse layout file end at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

topCellCount = len(ly.top_cells())
dbunit = ly.dbu

print("Input Layout file name:", input_file)
# print("OriginalTopCellCount:", topCellCount)
# print("OriginalTopCellList:")
# print("Reorganize layout and get info start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


### add_cell into "DMO_TOP" cell
dmo_top = ly.add_cell("DMO_TOP")
# process multiple top cells
# modified there [top_index? top_cell?]
for top_index in ly.each_top_cell():
    t = pya.Trans(pya.Trans.R0, 0, 0)
    top_cell = ly.cell(top_index)
    if top_cell.name == "DMO_TOP":
        # print("y")
        continue
    else:
        # print("Top_cell_index:", top_index, "Top_cell_name:", top_cell.name, "Top_cell_bbox:", top_cell.bbox())
        ly.cell(dmo_top).insert(pya.CellInstArray(top_cell.cell_index(), t))

# print("NewTopCellName: DMO_TOP")
# print("NewTopCellBoudingBox:", ly.cell(dmo_top).bbox())  # dmo_top.bbox()
print("DBUnit:", dbunit)
print("LayerInfo:")

layer_indexes = ly.layer_indexes()
layer_indexes.sort()
for idx in layer_indexes:
    layer_info = ly.get_info(idx)
    print(layer_info)

print("Reorganize layout and get info end at", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

# clipping and save result
boxes = []
# print("Read clipping boxes start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
with open(coor_file, 'r') as f:
    for line in f:
        if 'X' in line:
            continue
        # coors and radius 's unit is um
        coors = line.strip().split()
        x = int(float(coors[0]) / dbunit)
        y = int(float(coors[1]) / dbunit)
        r = int(float(radius) / dbunit)
        clip_box = pya.Box(x - r, y - r, x + r, y + r)
        boxes.append(clip_box)
# print("Read clipping boxes end at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
print("Clipping start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

# print("Store "+store_format+" files start at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# clip box in one area one layer and save to file
for box in boxes:
    # print(box)
    c_idx = ly.clip(dmo_top, box)
    # if the cut result is empty then next
    if ly.cell(c_idx).bbox().empty()==True:
        # print("there is nothing in this area this layer:",output[:-4] + "_" + str(format(x,'.2f')) + "-" + str(format(y,'.2f')))
        print("there is nothing in this area this layout's layer:",
              output[:-4] + "_" + str(int(x)) + "-" + str(int(y)) + "-" + str(int(float(radius) / dbunit)) + "."+store_format[:-2])
        continue
    # cliped box into new cell CLIP_TOP
    ############## change name to CLIP
    top = ly.add_cell("CLIP")
    t = pya.Trans(pya.Trans.R0, 0, 0)
    ly.cell(top).insert(pya.CellInstArray(c_idx, t))

    ############ if flatten:
    ly.cell(top).flatten(True)

    # save output
    options = pya.SaveLayoutOptions()
    options.add_cell(top)
    # options.set_format_from_filename(output)
    # options.format = "OASIS"

    # # um
    # x = (ly.cell(top).bbox().left + ly.cell(top).bbox().right) / 2.0 * dbunit
    # y = (ly.cell(top).bbox().bottom + ly.cell(top).bbox().top) / 2.0 * dbunit
    #
    # outputFile = output[:-4] + "_" + str(format(x,'.2f')) + "-" + str(format(y,'.2f')) + "."+store_format[:-2]

    # (x,y)
    x = (ly.cell(top).bbox().left + ly.cell(top).bbox().right) / 2.0
    y = (ly.cell(top).bbox().bottom + ly.cell(top).bbox().top) / 2.0

    outputFile = output[:-4] + "_" + str(int(x)) + "-" + str(int(y)) + "-" + str(int(float(radius) / dbunit)) + "."+store_format[:-2]

    ly.write(outputFile, options)
    print("Output Layout file name:", outputFile)
    # print("ResultTopCellName: CLIP_TOP")

    ly.delete_cell(top)    #

# print("Store "+store_format+" files end at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
print("End at:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
