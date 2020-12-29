import pya

ly = pya.Layout()
ly.read(input_file)
# to get layer info
print("LAYER")
layer_indexes = ly.layer_indexes()
layer_indexes.sort()
for idx in layer_indexes:
    layer_info = ly.get_info(idx)
    print(layer_info)

