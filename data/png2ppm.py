from PIL import Image
import sys

image = Image.open(sys.argv[1])
w, h = image.size
f_out = open(sys.argv[2], "w")
f_out.write("P3\n%d %d\n255\n" % (w, h));
for i in range(h):
    for j in range(w):
        f_out.write("%d %d %d " % image.getpixel((j, i)))
