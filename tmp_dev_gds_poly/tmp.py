import klayout.db as kdb
import os
from shapely import Polygon

file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(file_dir)

target_layer = (210,0)

layout = kdb.Layout()
layout.read("./test.gds")
top_cell = layout.top_cell()

shapes = top_cell.shapes(layout.layer(target_layer[0], target_layer[1]))

for shape in shapes:
	kdb_poly = shape.polygon
	shapely_poly = Polygon([(point.x * layout.dbu, point.y * layout.dbu) for point in kdb_poly.each_point_hull()])