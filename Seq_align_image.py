from PIL import Image, ImageDraw
from Bio import SeqIO
from operator import itemgetter
from itertools import groupby


def diff(a, b):
    index_list = []
    for index, (char1, char2) in enumerate(zip(a, b)):
        if char1 != char2:
            index_list.append(index)
    return index_list


fasta_sequences = SeqIO.parse(open("aligned_seqs.fasta"),'fasta')

seq1 = ''
seq2 = ''
count = 0
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    if count == 0:
        seq1 = sequence
    elif count == 1:
        seq2 = sequence
    count +=1

d = diff(seq1, seq2)

ranges = []
for k, g in groupby(enumerate(d), lambda x: x[0]-x[1]):
    group = (map(itemgetter(1), g))
    group = list(map(int, group))
    ranges.append((group[0], group[-1]))

if len(seq1) >= len(seq2):
    image_len = len(seq1)
else:
    image_len = len(seq2)
# size of image
canvas = (image_len, 600)

# scale ration
scale = 6
thumb = canvas[0]/scale, canvas[1]/5

# init canvas
im = Image.new('RGBA', canvas, (255, 255, 255, 255))
draw = ImageDraw.Draw(im)

draw.rectangle([0, 0, (len(seq1)), 150], fill="gray")
draw.rectangle([0, 250, (len(seq2)), 400], fill="gray")

# draw rectangles that don't match the reference sequence
for start, end in ranges:
    draw.rectangle([start, 250, end, 400], fill="red")

# make thumbnail
im.thumbnail(thumb)

# save image
im.save('./seq_align_PIL.png')
