import matplotlib
matplotlib.use('Agg') # hack to avoid DISPLAY error when in the console
import matplotlib.pyplot as plt
from matplotlib import cm
from drdcommon import *
import pylab

def boxplot(data, title="title here", ofn="boxplot.png", y_limit=None):
  fig = plt.figure()
  fig.set_size_inches(12,8)
  ax = fig.add_subplot(1,1,1)
  if y_limit:
    xlim(0, y_limit)
  ax.set_title(title)
  ax.boxplot(data)
  plt.savefig(ofn, dpi=150)

def barplot(y, labels, title="title_here", ofn="output.png"):
  fig = plt.figure()
  fig.set_size_inches(12,8)
  ax  = fig.add_subplot(1,1,1)

  N   = len(y)
  ind = range(N)

  #ax.set_yscale('log')
  ax.bar(ind, y,  align='center')
  ax.set_ylabel('Counts')
  ax.set_title(title, fontstyle='italic')
  ax.set_xticks(ind)
  ax.set_xticklabels(labels, rotation=45)
  #fig.autofmt_xdate()
  #p.show()
  #plt.savefig(ofn, dpi=150, bbox_inches='tight')
  plt.savefig(ofn, dpi=150)

def scatter_plot(ofn, x, y, title="title here", xlabel="xlabel", ylabel="ylabel", dot_size=10):
  fig = plt.figure()
  ax = fig.add_subplot(1, 1, 1)
  #ax.set_xticks(ticks=[i for i in range(1,61) if i % 5 == 0 or i == 1])
  #ax1.subplots(1)
  ax.scatter(x, y, dot_size)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.grid(True)
  plt.title(title)

  #fig.subplots_adjust(wspace=1, hspace=0.3)
  if in_ipython():
    print "Use show(%d) to display plot. " % 0
  else:
    plt.savefig(ofn, dpi=400, bbox_inches='tight')

class Heatmap:
  def __init__(self, values, cb_labels, xticks=[], yticks=[], fn_out="heat.png"):
    self.values = values
    self.xticks = xticks
    self.yticks = yticks
    self.fn_out = fn_out
    self.cb_labels = cb_labels

  def plot(self):
    #f   = pylab.Figure()
    #ax = f.add_subplot(111)
    ##plt.gcf().subplots_adjust(bottom=0.25, left=0.9, right=1)
    #m = ax.matshow(self.values)
    #pylab.xticks(range(len(self.xticks)), self.xticks, rotation=90)
    #pylab.yticks(range(len(self.yticks)), self.yticks)
    ##pylab.margins()
    #cb = f.colorbar(m, ticks=[0,1,2,3,4])
    ##cb.set_label(self.cb_label)
    #cb.ax.set_yticklabels(self.cb_labels)
    ##pylab.text(50, 0, 'yeah!', fontsize=12)
    #pylab.savefig(self.fn_out, dpi=100)

    font = {'family' : 'normal', 'weight' : 'normal', 'size' : 5}
    plt.rc('font', **font)

    fig = plt.figure()
    #plt.gcf().subplots_adjust(top=0.2)

    ax = fig.add_subplot(111)

    #cax = ax.matshow(self.values, interpolation='nearest', cmap=cm.afmhot)
    cax = ax.matshow(self.values, interpolation='nearest')
    pylab.xticks(range(len(self.xticks)), self.xticks, rotation=90)
    pylab.yticks(range(len(self.yticks)), self.yticks)
    #ax.set_title('xxxxxxxxx')

    cbar = fig.colorbar(cax, ticks=[0,1,2,3,4], orientation='horizontal')
    cbar.ax.set_xticklabels(self.cb_labels)# horizontal colorbar

    pylab.savefig(self.fn_out, dpi=200)
