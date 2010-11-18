from _conan import *
import os, subprocess

# Deep copy
def _undirected_graph_deepcopy(self, memo={}):
  new_g = undirected_graph_deep_copy(self)
  memo[id(self)] = new_g
  return new_g
undirected_graph.__deepcopy__ = _undirected_graph_deepcopy

def _directed_graph_deepcopy(self, memo={}):
  new_g = directed_graph_deep_copy(self)
  memo[id(self)] = new_g
  return new_g
directed_graph.__deepcopy__ = _directed_graph_deepcopy


# Helper functions for undirected graphs
def _undirected_graph_degree_distribution(self):
  return undirected_degree_distribution(self)
undirected_graph.degree_distribution = _undirected_graph_degree_distribution

def _undirected_graph_shortest_paths(self):
  return undirected_shortest_paths(self)
undirected_graph.shortest_paths = _undirected_graph_shortest_paths

def _undirected_graph_betweenness_centrality(self):
  return undirected_betweenness_centrality(self)
undirected_graph.betweenness_centrality = _undirected_graph_betweenness_centrality

def _undirected_graph_eigenvector_centrality(self):
  return undirected_eigenvector_centrality(self)
undirected_graph.eigenvector_centrality = _undirected_graph_eigenvector_centrality


# Helper functions for directed graphs
def _directed_graph_degree_distribution(self):
  return directed_degree_distribution(self)
directed_graph.degree_distribution = _directed_graph_degree_distribution

def _directed_graph_shortest_paths(self):
  return directed_shortest_paths(self)
directed_graph.shortest_paths = _directed_graph_shortest_paths

def _directed_graph_betweenness_centrality(self):
  return directed_betweenness_centrality(self)
directed_graph.betweenness_centrality = _directed_graph_betweenness_centrality

def _directed_graph_eigenvector_centrality(self):
  return directed_eigenvector_centrality(self)
directed_graph.eigenvector_centrality = _directed_graph_eigenvector_centrality


# Save graph as image
def graph_save_as(self, filename, program = 'neato'):
  image_type = filename[-3:]
  base_filename = filename[:-4]
  dot_filename = base_filename + '.dot'

  self.write_dot(dot_filename)
  subprocess.Popen(program + ' -T' + image_type + ' -o ' + filename + ' ' + dot_filename, \
                   shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
undirected_graph.save_as = graph_save_as
directed_graph.save_as = graph_save_as


# Default values
use_matplotlib = True # to display degree distribution
use_pyqt4 = False # to display the graph in a PyQt4 window
use_pygtk = True # to display the graph in a PyGtk window


# Parse config file
try:
  config_file = open(os.environ['HOME'] + '/.conanrc', 'r')

  for line in config_file:

    if line.startswith('USE_MATPLOTLIB'):
      value = line.split('=')[1].strip().lower()
      if value == 'false':
        use_matplotlib = False

    elif line.startswith('USE_PYQT4'):
      value = line.split('=')[1].strip().lower()
      if value == 'true':
        use_pyqt4 = True

    elif line.startswith('USE_PYGTK'):
      value = line.split('=')[1].strip().lower()
      if value == 'false':
        use_pygtk = False

  config_file.close()

except:
  pass


# Save distribution as image
if use_matplotlib:
  try:
    import math
    import numpy
    import matplotlib
    if 'DISPLAY' not in os.environ:
      matplotlib.use("Agg")
    import pylab
    
    def configure_graph(self):
      pylab.clf()

      nzd = self.non_zero_degrees()
      max_nzd = max(nzd)
      min_nzd = min(nzd)

      N = max_nzd - min_nzd # number of degrees between the minimum and maximum degree with a non-zero value
      xlocations = numpy.arange(N + 1) + 0.5

      width = .9 # the width of the bars

      P = [self.P(k) for k in range(min_nzd, max_nzd + 1)]

      pylab.bar(xlocations, P, width)

      if N < 15:
        pylab.xticks(xlocations + (width / 2), range(min_nzd, max_nzd + 1))
      else:
        pylab.xticks(numpy.arange(0, N + 1, int(N) / 15) + 0.5 + (width / 2), range(min_nzd, max_nzd + 1, int(N) / 15))

      pylab.ylim(ymin=0, ymax=max(P) * 1.2)
      xmin = -(width / 2)
      xmax = N + 2 + (width / 2)
      pylab.xlim(xmin=xmin, xmax=xmax)

      return
    
    undirected_degree_distribution.configure_graph = configure_graph
    directed_degree_distribution.configure_graph = configure_graph
    
    def save_distribution_as(self, filename):
      self.configure_graph()
      pylab.savefig(filename)
      return
    
    undirected_degree_distribution.save_as = save_distribution_as
    directed_degree_distribution.save_as = save_distribution_as
  
  except:
    pass


if 'DISPLAY' not in os.environ and (use_matplotlib or use_pyqt4 or use_pygtk):
  use_pyqt4 = False
  use_pygtk = False

  def display_not_found(self, program='', height=0, width=0):
    print 'Could not find DISPLAY environment variable... graphical output disabled.'

  undirected_graph.view = display_not_found
  directed_graph.view = display_not_found
  undirected_degree_distribution.view = display_not_found
  directed_degree_distribution.view = display_not_found

else: # DISPLAY is available

  if use_pyqt4 and not use_pygtk:
    try:
      import sys
      from PyQt4 import QtGui, QtCore, QtSvg

      class ViewImageWindow(QtGui.QMainWindow):
        def __init__(self, image_filename, width=0, heigth=0, parent=None):
          QtGui.QWidget.__init__(self, parent)

          self.setWindowTitle('Conan Graph Viewer')
          self.resize(width, height)

          if image_filename[-3:] == 'svg':
            svg_renderer = QtSvg.QSvgRenderer(image_filename)

            if svg_renderer.isValid():
              self.pixmap = QtGui.QPixmap(width, height)
              painter = QtGui.QPainter(self.pixmap)
              svg_renderer.render(painter)
            else:
              self.statusBar().showMessage('Svg invalid')

          else:
            self.pixmap = QtGui.QPixmap(image_filename)

          self.label = QtGui.QLabel(self)
          self.label.setPixmap(self.pixmap)
          self.setCentralWidget(self.label)

          self.exit = QtGui.QAction(
              QtGui.QIcon('/usr/share/icons/oxygen/22x22/actions/application-exit.png'), 'Exit', self)
          self.exit.setShortcut('Ctrl+Q')
          self.connect(self.exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))

          self.toolbar = self.addToolBar('Main')
          self.toolbar.addAction(self.exit)


      def view_graph(self, program='neato', width=None, height=None, keep_ratio=True):
        image_filename = '/tmp/conan_graph.svg'
        self.save_as(image_filename, program)

        output, err = subprocess.Popen('grep -e "<svg width" ' + image_filename, shell=True, stdout=subprocess.PIPE).communicate()
        real_height = output.split('"')[3]
        real_width = output.split('"')[1]
        ratio = real_height / real_width

        screen_height = gtk.gdk.screen_height()
        screen_width = gtk.gdk.screen_width()

        if height == None:
          if real_height < screen_height:
            height = real_height
          else:
            height = screen_height * .9

        if width == None:
          if real_width < screen_width:
            width = real_width
          else:
            width = screen_width * .9

        if keep_ratio:
          if height < ratio * width:
            height = ratio * width
          else:
            width = ratio * height

        app = QtGui.QApplication(sys.argv)
        window = ViewImageWindow(image_filename, width, height)
        window.show()
        return window, app

      undirected_graph.view = view_graph
      directed_graph.view = view_graph

    except:
      pass

  if use_pygtk:
    try:
      import pygtk
      pygtk.require('2.0')
      import gtk

      class ViewImageWindow:
        def __init__(self, image_filename, width=0, heigth=0):
          self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
          self.window.connect("destroy", gtk.main_quit)
          self.window.set_border_width(10)
          self.window.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse('black'))

          self.main_vbox = gtk.VBox()
          self.window.add(self.main_vbox)

          self.image = gtk.Image()
          if width == 0 or heigth == 0:
            self.image.set_from_file(image_filename)
          else:
            self.pixbuf = gtk.gdk.pixbuf_new_from_file(image_filename)
            self.pixbuf = self.pixbuf.scale_simple(width, heigth, gtk.gdk.INTERP_BILINEAR)
            self.image.set_from_pixbuf(self.pixbuf)
          self.main_vbox.pack_start(self.image)

          self.window.show_all()

        def main(self):
          gtk.main()

      def view_graph(self, program='neato', width=None, height=None, keep_ratio=True):
        image_filename = '/tmp/conan_graph.svg'
        self.save_as(image_filename, program)

        output, err = subprocess.Popen('grep -e "<svg width" ' + image_filename, shell=True, stdout=subprocess.PIPE).communicate()
        real_height = int(output.split('"')[3].strip('pt'))
        real_width = int(output.split('"')[1].strip('pt'))
        ratio = real_height / real_width

        screen_height = gtk.gdk.screen_height()
        screen_width = gtk.gdk.screen_width()

        if height == None:
          if real_height < screen_height:
            height = real_height
          else:
            height = screen_height * .9

        if width == None:
          if real_width < screen_width:
            width = real_width
          else:
            width = screen_width * .9

        if keep_ratio:
          if height < ratio * width:
            height = ratio * width
          else:
            width = ratio * height

        window = ViewImageWindow(image_filename, width, height)
        window.main()

        return

      undirected_graph.view = view_graph
      directed_graph.view = view_graph

    except:
      pass


  if use_matplotlib:
    try:
      def view_distribution(self):
        self.configure_graph()
        pylab.show()
        return

      undirected_degree_distribution.view = view_distribution
      directed_degree_distribution.view = view_distribution

    except:
      pass

