from PIL import Image
import numpy 


img = Image.open('horse.jpg').convert('1')  # convert image to 8-bit grayscale
WIDTH, HEIGHT = img.size
img.save('bwimage.png')

data = list(img.getdata()) # convert image data to a list of integers
# convert that to 2D list (list of lists of integers)
data = [data[offset:offset+WIDTH] for offset in range(0, WIDTH*HEIGHT, WIDTH)]

# At this point the image's pixels are all in memory and can be accessed
# individually using data[row][col].
arr = numpy.array(data)
arr = arr.clip(max=1)

f = open("binaryImage.txt", "a")

for row in arr:
    # f.write(' '.join('{:3}'.format(value) for value in row))
    # f.write(str(row))
    for element in row:
        f.write(str(element))
    f.write('\n')

f.close()

