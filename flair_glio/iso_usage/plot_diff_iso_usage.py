import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import sys
import seaborn as sns

import matplotlib.patches as mplpatches

try:
    psl = open(sys.argv[1])
    counts_matrix = open(sys.argv[2])
    genename = sys.argv[3]
except:
    sys.stderr.write('python script.py isoforms.psl counts_matrix.txt genename \n')
    sys.exit(1)

try:
    plt.style.use('./BME163_arial')
except:
    # sys.stderr.write('cannot find stylesheet\n')
    pass

def parse_psl(psl, names=False, plotany=False,keepiso=set()):
    info = []
    usednames = []
    lowbound, upbound = 1e9, 0
    chrom = ''
    color_order = []
    for line in psl:
        line = line.rstrip().split('\t')
        if line[9] not in keepiso:
            continue
        color_order += [keepiso[line[9]]]
        lowbound = min(lowbound, int(line[15]))
        upbound = max(upbound, int(line[16]))
        blocksizes = [int(n) for n in line[18].split(',')[:-1]]
        blockstarts = [int(n) for n in line[20].split(',')[:-1]]
        info += [[blocksizes, blockstarts, 0]]
        # info += [[blocksizes, blockstarts, int(line[22])]] # [list of blocksizes, list of blockstarts, productivity_bool]
        usednames += [line[9]]
    strand = line[8]
    upbound += 100
    lowbound -= 100
    for i in range(len(info)):
        info[i][1] = [n - lowbound for n in info[i][1]]
    if names:
        return info, usednames
    else:
        return info, lowbound, upbound, strand, color_order, usednames

def pack(data, rev=True, color=False, tosort = True):
    starts = [max(d[1]) for d in data] if rev else [min(d[1]) for d in data] # sort by right or left end
    if tosort:
        data = [d for (s,d) in sorted(zip(starts, data))]
    else:
        data = [d for s, d in zip(starts, data)]
    packed = [[ data[0] ]]
    ends = [max(data[0][1]) + data[0][0][data[0][1].index(max(data[0][1]))]]
    for i in range(1, len(data)):
        min_start = min(data[i][1])
        end = max(data[i][1]) + data[i][0][data[i][1].index(max(data[i][1]))]
        pos = -1
        for j in range(len(packed)):
            if ends[j] + 5 < min_start:  # added the plus 5 for spacing between reads
                pos = j  # pack this read with the read at position j
                break
        if pos >= 0:
            packed[pos] += [data[i]]
            ends[pos] = end
        else:
            packed += [[data[i]]]
            ends += [end]
    return packed

def plot_blocks(data, panel, names, utr=False, height=.5, l=0.8, color=False):
    panel.set_xlim(1, upper - lower)
    if strand == '-':
        panel.set_xlim(upper - lower, 1)
    panel.set_ylim(-1, len(data) * 2)
    # panel.set_ylim(-1, 10)

    panel.tick_params(axis='both', which='both',\
                       bottom=False, labelbottom=False,\
                       left=False, labelleft=False,\
                       right=False, labelright=False,\
                       top=False, labeltop=False)

    di = 0  # data index
    ni = 0  # nameindex
    for i in range(len(data)*2):  # each line
        if i % 2 == 0:
            if strand == '-':
                panel.text(data[di][0][1][0], i-height/2 + 1, names[ni], fontsize=6, ha='right', va='center')
            else:
                panel.text(data[di][0][1][0], i-height/2 + 1, names[ni], fontsize=6, ha='left', va='center')
            ni += 1
            continue
        read = data[di]
        # di += 1
        for j in range(len(read)):  # each read on each line
            line = read[j]
            sizes = line[0]
            starts = line[1]
            panel.plot([min(starts)+sizes[0], max(starts)], [i - 1]*2, 'k-', lw=l)
            for k in range(len(sizes)):  # each block of each read
                if utr and line[2][k]:
                    rectangle1 = mplpatches.Rectangle([starts[k], i-height/2], \
                        sizes[k], height, facecolor='white', linewidth=0, zorder=11)
                    rectangle2 = mplpatches.Rectangle([starts[k], i-.25/2],\
                        sizes[k], 0.25, facecolor='black', linewidth=0, zorder=12)
                    panel.add_patch(rectangle1)
                    panel.add_patch(rectangle2)
                elif color:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2 - 1], \
                        sizes[k], height, facecolor=color[di], linewidth=0, zorder=10)
                    panel.add_patch(rectangle)
         
                else:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2], \
                        sizes[k], height, facecolor='black', linewidth=0)
                    panel.add_patch(rectangle)
        di += 1


flatui = ["#3498db", "#34495e", "#2ecc71", "#e74c3c"]
colors = ["windows blue", "faded green", "dusty purple",'amber']
color_palette = sns.xkcd_palette(colors) + sns.color_palette(flatui)
gray = sns.xkcd_palette(["greyish"])
# color_palette = sns.color_palette("BrBG", 12)

plt.figure(figsize=(7,6))

panel = plt.axes([0.11, 0.11, .88, 0.88], frameon=False)  # plotting the proportion of expression
panel.tick_params(axis='both',which='both', \
                   bottom=True, labelbottom=True, \
                   left=True, labelleft=True, \
                   right=False, labelright=False, \
                   top=False, labeltop=False, labelsize=8)

colori = 0
minreads = 50
keepiso = {}  # isoforms that they have a sufficient proportion of reads mapping to them
sample_ids = counts_matrix.readline().rstrip().split('\t')[1:]
proportions = [['lowexpr']+[0]*len(sample_ids)+gray]  # the minor isoform bar
totals = [0]*len(sample_ids)
for line in counts_matrix:
    line = line.rstrip().split('\t')
    if genename not in line[0]:
        continue
    counts = [float(x) for x in line[1:]]
    for i in range(len(sample_ids)):
        totals[i] += counts[i]
    if all(x < minreads for x in counts):
        for i in range(len(sample_ids)):
            proportions[0][i+1] += counts[i]   # add to gray bar
        continue
    if colori == len(color_palette):
        sys.stderr.write('more isoforms than there are colors\n')
        sys.exit()
    keepiso[line[0]] = color_palette[colori]
    proportions += [[line[0]]+counts+[color_palette[colori]]]
    colori += 1
xlim = len(sample_ids) + 0.5


#if len(proportions) == 1:
#    sys.stderr.write('need more than 1 isoform with sufficient representation\n')
#    sys.exit()

proportions = sorted(proportions, key=lambda x:x[1])[::-1]
heights = [0]*len(sample_ids)
for iso in proportions:
    for i in range(len(sample_ids)):
        percentage = iso[1+i]/totals[i]*100
        rectangle = mplpatches.Rectangle([i+1-0.4, heights[i]], \
                    width=0.8, height=percentage, facecolor=iso[-1], linewidth=0)
        panel.add_patch(rectangle)
        if percentage >= 8:  # percentage
            panel.text(i+1, heights[i] + percentage/2, str(round(percentage,1))+'%',\
             fontsize=10, ha='center',va='center',color='white')
        if percentage >= 12.25:  # read num
            if iso[1+i] == 1:
                panel.text(i+1.36, heights[i]+1.75, str(int(iso[1+i]))+' read',\
                 fontsize=6, ha='right',va='center',color='white')
            else:
                panel.text(i+1.35, heights[i]+1.75, str(int(iso[1+i]))+' reads',\
                 fontsize=6, ha='right',va='center',color='white')

        heights[i] += percentage

panel.plot([-1, xlim], [0, 0], 'r--', lw=.75,color='black',alpha=.3,zorder=0)
panel.plot([-1, xlim], [20, 20], 'r--', lw=.5,color='black',alpha=.3,zorder=0)
panel.plot([-1, xlim], [40, 40], 'r--', lw=.5,color='black',alpha=.3,zorder=0)
panel.plot([-1, xlim], [60, 60], 'r--', lw=.5,color='black',alpha=.3,zorder=0)
panel.plot([-1, xlim], [80, 80], 'r--', lw=.5,color='black',alpha=.3,zorder=0)
panel.plot([-1, xlim], [100, 100], 'r--', lw=.75,color='black',alpha=.3,zorder=0)

panel.set_xticks(np.arange(1,xlim))
panel.set_xticklabels(sample_ids, rotation=20, ha='right')
panel.set_xlim(0.5, xlim)
panel.set_ylim(0, 100)

plt.savefig(genename+'_proportion.pdf', transparent=True)

# isoform structures
fig_0 = plt.figure(figsize=(10, 2))
panel = plt.axes([0.005, 0.015, .99, 0.97], frameon=True)  # annotation

isoforms, lower, upper, strand, color_order, names = parse_psl(psl,keepiso=keepiso)
packed = pack(isoforms, rev=False, tosort=False)

plot_blocks(packed, panel, names, color=color_order, l=1)  # plot as is

if strand == '-':
    panel.text(lower, 0.5, "5' end (-)", fontsize=6, ha='left',va='center')
    panel.invert_xaxis()
else:
    panel.text(lower, 0.5, "5' end (+)", fontsize=6, ha='left',va='center')

plt.savefig(genename+'_bars.pdf', transparent=True)
