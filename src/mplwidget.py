# this is needed to display the plot in a Qt widget
# Python application for automation of tandem aerosol classifier experiments and data inversion
# Copyright (C) 2023 Cambustion, The University of Alberta, NRC Canada
# Original Authors: Jonathan Symonds, Morteza Kiasadegh, Tim Sipkens
# See LICENSE for details

import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
    """
    Object-oriented canvas for Matplotlib drawing.
    """
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        # Create a specific Figure instance, completely separate from pyplot
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111)
        
        super().__init__(self.fig)
        self.setParent(parent)
        
        # Ensure the canvas expands to fill the available space
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.updateGeometry()


class MplWidget(QtWidgets.QWidget):
    """
    A standalone Qt widget containing a Matplotlib canvas and interactive toolbar.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Instantiate canvas and toolbar
        self.canvas = MplCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Layout setup
        layout = QtWidgets.QVBoxLayout()
        # layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        # Remove default layout margins for a flush, professional look
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    def setup_axes(self, xlabel="", ylabel="", xscale='linear', yscale='linear', title=""):
        """Utility method to quickly configure the plot from the main app."""
        self.canvas.ax.clear()
        self.canvas.ax.set_title(title)
        self.canvas.ax.set_xlabel(xlabel)
        self.canvas.ax.set_ylabel(ylabel)
        self.canvas.ax.set_xscale(xscale)
        self.canvas.ax.set_yscale(yscale)
        
        # Let Matplotlib calculate the tightest bounds automatically
        self.canvas.fig.set_layout_engine('tight') 
        self.canvas.draw_idle()


