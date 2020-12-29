# arguments:
# input_file layer datatype output_file
import pya

ly = pya.Layout()
layer_map = pya.LayerMap()

spec = pya.LayerInfo(int(layer), int(datatype))
layer_map.map(spec, ly.insert_layer(spec))

opt = pya.LoadLayoutOptions()
opt.set_layer_map(layer_map, False)
ly.read(input_file, opt)

ly.write(output_file)
