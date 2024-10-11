# this is needed to display the plot in a Qt widget
# Python application for automation of tandem aerosol classifier experiments and data inversion
# Copyright (C) 2023 Cambustion, The University of Alberta, NRC Canada
# Original Authors: Jonathan Symonds, Morteza Kiasadegh, Tim Sipkens
# See LICENSE for details

from PyQt5 import QtWidgets
from matplotlib.pyplot import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')

# Matplotlib canvas class to create figure


class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        self.fig.tight_layout()
        self.ax = self.fig.add_subplot(111)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(
            self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)

# Matplotlib widget


class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)


