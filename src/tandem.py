# Main UI file
# Python application for automation of tandem aerosol classifier experiments and data inversion
# Copyright (C) 2024 Cambustion, The University of Alberta, NRC Canada
# Original Authors: Jonathan Symonds, Morteza Kiasadegh, Tim Sipkens
# See LICENSE for details

from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, Qt, QSettings
from PyQt5.uic import loadUi
import sys
import instruments
import mplwidget
import numpy as np
import time
import datetime
import csv
import platform
import ctypes
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.colors as mcolors
import matplotlib as mpl


version = "0.2"

app = QtWidgets.QApplication(sys.argv)


class Ui(QtWidgets.QMainWindow):  # main GUI object
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('tandem.ui', self)
        if (platform.system() == 'Windows'):  # makes taskbar icon work in Windows 7 onwards
            appid = 'Cambustion.Tandem.'+version
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                appid)
        self.setFixedSize(self.width(), self.height())
        self.setWindowFlags(self.windowFlags() |
                            Qt.MSWindowsFixedSizeDialogHint)
        self.settings = QSettings("Cambustion", "Tandem")
        self.getSettings()
        self.graph = self.plotWidget.canvas.ax
        self.plot = self.MplWidget.canvas.ax
        self.plot_live = self.MplWidget_live.canvas.ax
        self.plot_live.set_xscale('log')
        self.plot_live.set_xlabel('Second Classifier setpoints')
        self.plot_live.set_ylabel('Particle Concentration [Counts]')
        self.plot.set_xscale('log')
        self.plot.set_xlabel('Second Classifier setpoints')
        self.plot.set_ylabel('Particle Concentration [Counts]')
        self.plot2 = self.plot.twinx()
        self.plot2.set_ylabel('Bypass Particle Concentration [Counts]')
        self.graph.set_yscale('log')
        self.graph.set_xscale('log')
        self.graph.set_xlabel('First Classifier')
        self.graph.set_ylabel('Second Classifier')
        self.fig = self.graph.get_figure()  # Get the Figure object
        self.isScanning = False
        self.show()

    def getSettings(self):
        self.actionSave_Settings_on_Exit.setChecked(
            self.settings.value("autoSaveSettings", False, type=bool))

        self.firstClassifierList.setCurrentRow(
            self.settings.value("first/Classifier", 0, type=int))
        self.firstIsEth.setChecked(
            self.settings.value("first/IsEth", True, type=bool))
        self.firstIsSerial.setChecked(
            not self.settings.value("first/IsEth", True, type=bool))
        self.firstIP.setPlainText(
            self.settings.value("first/IP", "192.168.1.2"))
        self.firstPort.setPlainText(self.settings.value("first/Port", "1"))
        self.firstQa.setPlainText(self.settings.value("first/Qa", "1.5"))
        self.firstRorQsh.setPlainText(
            self.settings.value("first/RorQsh", "15.0"))
        self.firstFrom.setPlainText(self.settings.value("first/From", "200.0"))
        self.firstTo.setPlainText(self.settings.value("first/To", "1000.0"))
        self.firstPerDecade.setPlainText(
            self.settings.value("first/PerDecade", "16"))
        self.firstDelay.setPlainText(self.settings.value("first/Delay", "2.0"))
        self.firstPolarity.setChecked(
            self.settings.value("first/Polarity", True, type=bool))

        self.secondClassifierList.setCurrentRow(
            self.settings.value("second/Classifier", 1, type=int))
        self.secondIsEth.setChecked(
            self.settings.value("second/IsEth", True, type=bool))
        self.secondIsSerial.setChecked(
            not self.settings.value("second/IsEth", True, type=bool))
        self.secondIP.setPlainText(
            self.settings.value("second/IP", "192.168.1.3"))
        self.secondPort.setPlainText(self.settings.value("second/Port", "2"))
        self.secondQa.setPlainText(self.settings.value("second/Qa", "1.5"))
        self.secondRorQsh.setPlainText(
            self.settings.value("second/RorQsh", "10.0"))
        self.secondFrom.setPlainText(
            self.settings.value("second/From", "200.0"))
        self.secondTo.setPlainText(self.settings.value("second/To", "1000.0"))
        self.secondPerDecade.setPlainText(
            self.settings.value("second/PerDecade", "16"))
        self.secondScanner.setChecked(
            self.settings.value("second/Scanner", False, type=bool))
        self.Bypass.setChecked(self.settings.value(
            "second/doBypass", False, type=bool))
        self.secondDelay.setPlainText(
            self.settings.value("second/Delay", "2.0"))
        self.secondHighFlow.setChecked(
            self.settings.value("second/HighFlow", True, type=bool))
        self.secondLowFlow.setChecked(
            not self.settings.value("second/HighFlow", True, type=bool))
        self.secondPolarity.setChecked(
            self.settings.value("second/Polarity", True, type=bool))
        self.secondScanUpTime.setPlainText(
            self.settings.value("second/ScanUpTime", "240"))
        self.secondLowerRange.setCurrentIndex(
            self.settings.value("second/LowerRange", 0, type=int))
        self.secondUpperRange.setCurrentIndex(
            self.settings.value("second/UpperRange", 191, type=int))

        self.scanDelay.setPlainText(
            self.settings.value("second/ScanDelay", "8.0"))
        self.scanAve.setPlainText(self.settings.value("second/ScanAve", "4"))

        self.cpcList.setCurrentRow(
            self.settings.value("cpc/Type", 1, type=int))
        self.cpcIP.setPlainText(self.settings.value("cpc/IP", "192.168.1.4"))
        self.cpcIsEth.setChecked(
            self.settings.value("cpc/IsEth", False, type=bool))
        self.cpcIsSerial.setChecked(
            not self.settings.value("cpc/IsEth", False, type=bool))
        self.cpcPort.setPlainText(self.settings.value("cpc/Port", "3"))
        self.average.setPlainText(self.settings.value("average", "1"))

        self.fileName = self.settings.value("RawFile", "../data/RawData.txt")
        self.rawFileDisplay.setText("Raw Data File: "+self.fileName)

    def saveSettings(self):
        self.settings.setValue(
            "autoSaveSettings", self.actionSave_Settings_on_Exit.isChecked())

        self.settings.setValue(
            "first/Classifier", self.firstClassifierList.currentRow())
        self.settings.setValue("first/IsEth", self.firstIsEth.isChecked())
        self.settings.setValue("first/IP", self.firstIP.toPlainText())
        self.settings.setValue("first/Port", self.firstPort.toPlainText())
        self.settings.setValue("first/Qa", self.firstQa.toPlainText())
        self.settings.setValue("first/RorQsh", self.firstRorQsh.toPlainText())
        self.settings.setValue("first/From", self.firstFrom.toPlainText())
        self.settings.setValue("first/To", self.firstTo.toPlainText())
        self.settings.setValue(
            "first/PerDecade", self.firstPerDecade.toPlainText())
        self.settings.setValue("first/Delay", self.firstDelay.toPlainText())

        self.settings.setValue(
            "first/Polarity", self.firstPolarity.isChecked())

        self.settings.setValue("second/Classifier",
                               self.secondClassifierList.currentRow())
        self.settings.setValue("second/IsEth", self.secondIsEth.isChecked())
        self.settings.setValue("second/IP", self.secondIP.toPlainText())
        self.settings.setValue("second/Port", self.secondPort.toPlainText())
        self.settings.setValue("second/Qa", self.secondQa.toPlainText())
        self.settings.setValue(
            "second/RorQsh", self.secondRorQsh.toPlainText())
        self.settings.setValue("second/From", self.secondFrom.toPlainText())
        self.settings.setValue("second/To", self.secondTo.toPlainText())
        self.settings.setValue(
            "second/PerDecade", self.secondPerDecade.toPlainText())
        self.settings.setValue(
            "second/Scanner", self.secondScanner.isChecked())
        self.settings.setValue("second/doBypass", self.Bypass.isChecked())
        self.settings.setValue("second/Delay", self.secondDelay.toPlainText())
        self.settings.setValue(
            "second/HighFlow", self.secondHighFlow.isChecked())
        self.settings.setValue(
            "second/LowFlow", self.secondLowFlow.isChecked())
        self.settings.setValue(
            "second/Polarity", self.secondPolarity.isChecked())
        self.settings.setValue("second/ScanUpTime",
                               self.secondScanUpTime.toPlainText())
        self.settings.setValue(
            "second/ScanDelay", self.scanDelay.toPlainText())
        self.settings.setValue("second/ScanAve", self.scanAve.toPlainText())
        self.settings.setValue("cpc/Type", self.cpcList.currentRow())
        self.settings.setValue("cpc/IsEth", self.cpcIsEth.isChecked())
        self.settings.setValue("cpc/IP", self.cpcIP.toPlainText())
        self.settings.setValue("cpc/Port", self.cpcPort.toPlainText())
        self.settings.setValue("average", self.average.toPlainText())

        self.settings.setValue("RawFile", self.fileName)

    def reset(self):
        if (QtWidgets.QMessageBox.question(self, "Tandem", "Do you want to return all settings to their defaults?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes):
            self.settings.clear()
            self.getSettings()

    def closeEvent(self, event):
        if (self.actionSave_Settings_on_Exit.isChecked()):
            self.saveSettings()
        try:
            del self.firstClass
        except:
            pass
        try:
            del self.secondClass
        except:
            pass
        try:
            del self.cpc
        except:
            pass
        event.accept()  # let the window close

    def runScan(self):
        if (os.path.isfile(self.fileName)):
            if (QtWidgets.QMessageBox.question(self, "Tandem", "Raw data file exists. Overwrite?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No):
                return
        self.header_3082 = True
        self.makesecondClassX = True
        self.Scatterplot = False
        # Clear the existing plots
        self.plot.clear()
        self.plot2.clear()
        if self.plot is not None and self.plot.get_legend() is not None:
            self.plot.get_legend().remove()
        # Adjust the layout to accommodate the legend outside the box
        self.plotWidget.canvas.figure.tight_layout(
            rect=[0, 0, 1, 1])  # Adjust the right margin
        self.MplWidget.canvas.figure.tight_layout(
            rect=[0.03, 0.03, 0.87, 1])   # Adjust the right margin
        self.graph.clear()
        self.graph = self.plotWidget.canvas.ax
        self.plot = self.MplWidget.canvas.ax
        self.plot_live = self.MplWidget_live.canvas.ax
        self.plot_live.set_xscale('log')
        self.plot_live.set_xlabel('Second Classifier setpoints')
        self.plot_live.set_ylabel('Particle Concentration [Counts]')
        self.plot.set_xscale('log')
        self.plot.set_xlabel('Second Classifier setpoints')
        self.plot.set_ylabel('Particle Concentration [Counts]')
        self.plot2.set_ylabel('Bypass Particle Concentration [Counts]')
        self.graph.set_yscale('log')
        self.graph.set_xscale('log')
        self.graph.set_xlabel('First Classifier')
        self.graph.set_ylabel('Second Classifier')
        self.fig = self.graph.get_figure()  # Get the Figure object
        self.stopButton.setEnabled(True)
        self.startButton.setEnabled(False)
        self.statusbar.showMessage("Initialising first classifier...")
        app.processEvents()
        c = self.firstClassifierList.currentItem().text()
        if (c == "AAC"):
            self.firstClass = instruments.AAC
        elif (c == "CPMA"):
            self.firstClass = instruments.CPMA
        elif (c == "3081 DMA"):
            self.firstClass = instruments.TSI3080
        elif (c == "3082 DMA"):
            self.firstClass = instruments.TSI3082
        if (c == "3082 DMA"):
            self.firstClass = self.firstClass(self.firstIsSerial.isChecked(), self.firstIP.toPlainText(), self.firstPort.toPlainText(), float(self.firstFrom.toPlainText()), float(self.firstTo.toPlainText()),  int(
                self.firstPerDecade.toPlainText()), float(self.firstQa.toPlainText()), float(self.firstRorQsh.toPlainText()),  False, self.firstPolarity.isChecked(), 0, 0, 0, 0, 0, 0, 0, False, False, False, False, False)
        elif (c == "AAC"):
            self.firstClass = self.firstClass(self.firstIsSerial.isChecked(), self.firstIP.toPlainText(), self.firstPort.toPlainText(), float(self.firstFrom.toPlainText()), float(
                self.firstTo.toPlainText()), int(self.firstPerDecade.toPlainText()), float(self.firstQa.toPlainText()), float(self.firstRorQsh.toPlainText()), 0, 0, 0, 0, 0, 0, 0, False, False, False, False, False)
        else:
            self.firstClass = self.firstClass(self.firstIsSerial.isChecked(), self.firstIP.toPlainText(), self.firstPort.toPlainText(), float(self.firstFrom.toPlainText(
            )), float(self.firstTo.toPlainText()),  int(self.firstPerDecade.toPlainText()), float(self.firstQa.toPlainText()), float(self.firstRorQsh.toPlainText()))
        if (not self.firstClass.connected):
            self.barf("Connection Error to First Classifier")
            del self.firstClass
            return
        self.statusbar.showMessage("Initialising second classifier...")
        app.processEvents()
        self.isScanner = False
        self.varBins = False
        c = self.secondClassifierList.currentItem().text()
        if (c == "AAC"):
            self.isScanner = self.secondScanner.isChecked()
            self.varBins = self.variableBins.isChecked()
            self.secondClass = instruments.AAC
        elif (c == "CPMA"):
            self.secondClass = instruments.CPMA
        elif (c == "3081 DMA"):
            self.secondClass = instruments.TSI3080
        elif (c == "3082 DMA"):
            self.isScanner = self.secondScanner.isChecked()
            self.varBins = self.variableBins.isChecked()
            self.secondClass = instruments.TSI3082
        self.doBypass = self.Bypass.isChecked()

        if (c == "3082 DMA"):
            self.secondClass = self.secondClass(self.secondIsSerial.isChecked(), self.secondIP.toPlainText(), self.secondPort.toPlainText(), float(self.secondFrom.toPlainText()), float(self.secondTo.toPlainText()), int(self.secondPerDecade.toPlainText()), float(self.secondQa.toPlainText()), float(self.secondRorQsh.toPlainText()), self.secondHighFlow.isChecked(
            ), self.secondPolarity.isChecked(), int(self.secondScanUpTime.toPlainText()), int(self.secondLowerRange.currentIndex()), int(self.secondUpperRange.currentIndex()), self.PreFactor.toPlainText(), self.MM_Exponent.toPlainText(), self.f_lower.toPlainText(), self.f_upper.toPlainText(), self.particleIsWater.isChecked(), self.particleIsSoot.isChecked(), self.particleIsNone.isChecked(), self.varBins, self.isScanner)
        elif (c == "AAC"):
            self.secondClass = self.secondClass(self.secondIsSerial.isChecked(), self.secondIP.toPlainText(), self.secondPort.toPlainText(), float(self.secondFrom.toPlainText()), float(self.secondTo.toPlainText()), int(
                self.secondPerDecade.toPlainText()), float(self.secondQa.toPlainText()), float(self.secondRorQsh.toPlainText()), int(self.secondScanUpTime.toPlainText()), int(self.scanAve.toPlainText()), float(self.scanDelay.toPlainText()), self.PreFactor.toPlainText(), self.MM_Exponent.toPlainText(), self.f_lower.toPlainText(), self.f_upper.toPlainText(), self.particleIsWater.isChecked(), self.particleIsSoot.isChecked(), self.particleIsNone.isChecked(), self.varBins, self.isScanner)
        else:
            self.secondClass = self.secondClass(self.secondIsSerial.isChecked(), self.secondIP.toPlainText(), self.secondPort.toPlainText(), float(self.secondFrom.toPlainText(
            )), float(self.secondTo.toPlainText()), int(self.secondPerDecade.toPlainText()), float(self.secondQa.toPlainText()), float(self.secondRorQsh.toPlainText()))

        if (not self.secondClass.connected):
            self.barf("Connection Error to Second Classifier")
            del self.secondClass
            del self.firstClass
            return
        self.statusbar.showMessage("")

        self.arduino = instruments.Arduino
        if (self.Bypass.isChecked()):
            self.arduino = self.arduino(
                True, "", self.arduinoPort.toPlainText())

        app.processEvents()
        if (self.isScanner):
            self.cpc = instruments.DummyCPC
        else:
            self.statusbar.showMessage("Initialising CPC...")
            c = self.cpcList.currentItem().text()
            if (c == "Cambustion 5210"):
                self.cpc = instruments.CambustionCPC
            elif (c == "3022/25"):
                self.cpc = instruments.TSI30xx
            elif (c == "3775/76"):
                self.cpc = instruments.TSI377x
            elif (c == "375x"):
                self.cpc = instruments.TSI375x
            elif (c == "Dummy"):
                self.cpc = instruments.DummyCPC
            elif (c == "Magic"):
                self.cpc = instruments.MagicCpc

        self.cpc = self.cpc(self.cpcIsSerial.isChecked(
        ), self.cpcIP.toPlainText(), self.cpcPort.toPlainText())

        if (not self.cpc.connected):
            self.barf("Connection Error to CPC")
            del self.secondClass
            del self.firstClass
            del self.cpc
            return

        if (not self.firstClass.run()):
            self.barf("Can't start first classifier: " +
                      self.firstClass.lastResponse)
            del self.secondClass
            del self.firstClass
            del self.cpc
            return

        if (not self.secondClass.run()):
            self.barf("Can't start second classifier: " +
                      self.secondClass.lastResponse)
            del self.secondClass
            del self.firstClass
            del self.cpc
            return

        self.graph.set_xlabel(
            f"{self.firstClass.label}[{self.firstFromUnit.text()}]")
        self.graph.set_ylabel(
            f"{self.secondClass.label}[{self.secondFromUnit.text()}]")
        self.plot.set_xlabel(
            f"{self.secondClass.label}[{self.secondFromUnit.text()}]")

        if not (self.isScanner):
            self.plot.set_ylabel("Particle Concentration")
            self.X, self.Y = np.meshgrid(
                self.firstClass.X, self.secondClass.X, indexing='ij')
            self.Z = np.zeros(
                (self.firstClass.points, self.secondClass.points))
            self.G = np.zeros(self.secondClass.points)
        elif (self.secondClassifierList.currentItem().text() == "AAC"):
            self.plot.set_ylabel("N [Counts]")
            self.X, self.Y = np.meshgrid(
                self.firstClass.X, self.secondClass.X, indexing='ij')
            self.Z = np.zeros((1, self.secondClass.points))
            self.G = np.zeros(1)
        else:
            self.plot.set_ylabel("N [Counts]")

        # Create a colormap with a gradient from light to dark color
        light_color = np.array([0.0, 1.0, 0.0])
        dark_color = np.array([0.0, 0.0, 0.1])
        # Define the number of intermediate colors
        num_intermediate_colors = self.firstClass.points
        # Create a list of colors by linear interpolation
        colors = [light_color + (dark_color - light_color) * i /
                  num_intermediate_colors for i in range(num_intermediate_colors + 1)]
        # Create a LinearSegmentedColormap
        self.cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", colors)
        self.isScanning = True
        self.statusbar.showMessage("Scanning")
        self.timer = QTimer()
        self.timer.timeout.connect(self.check)
        if (not self.doBypass):
            if (not self.firstClass.next()):
                self.barf("Failed to set first classifier. Reason: " +
                          self.firstClass.lastResponse)
                return
        # else:
        #    if (not self.firstClass.enableBypass(3) and not self.firstClass.enableBypass(2)):
        #        self.barf("Neither classifier supports sending a bypass signal")
        #        return
        if (not self.secondClass.next()):
            self.barf("Failed to set second classifier. Reason: " +
                      self.secondClass.lastResponse)
            return
        if (self.doBypass):
            # self.firstClass.enableBypass(1)  # When using the analog output of the instrument
            self.arduino.enableBypass()  # When using the Arduino

        self.startRawLog()
        self.startRow = True
        self.finalBypass = False
        self.firstBypass = True
        self.disBypass = True
        self.timer.start(1)

    def startRawLog(self):
        self.file = open(self.fileName, 'w', newline='')
        self.file.write("Cambustion / University of Alberta\tTandem\tv"+version+"\t"+datetime.datetime.now(
        ).strftime("%Y-%m-%d%t%H:%M:%S")+"\tBypass scans:\t" + str(self.doBypass)+"\n")
        firstHeaders = self.firstClass.getHeader()
        firstHeaders['Classifier 1'] = self.firstClassifierList.currentItem().text()
        firstHeaders['Data points'] = self.firstClass.points
        firstHeaders['Data length'] = len(self.firstClass.fileFields)
        self.writer = csv.DictWriter(self.file, fieldnames=['Classifier 1']+['Data points']+[
                                     'Data length']+self.firstClass.headerFields, dialect='excel-tab')
        self.writer.writeheader()
        self.writer.writerow(firstHeaders)

        secondHeaders = self.secondClass.getHeader()
        if (self.secondClassifierList.currentItem().text() == "3082 DMA" and self.secondScanner.isChecked()):
            secondHeaders['Classifier 2'] = self.secondClassifierList.currentItem(
            ).text()
            secondHeaders['Start (nm)'] = self.secondLowerRange.itemText(
                int(self.secondClass.LowerSizeRange()))
            secondHeaders['End (nm)'] = self.secondLowerRange.itemText(
                int(self.secondClass.UpperSizeRange()))
            self.writer = csv.DictWriter(self.file, fieldnames=[
                                         'Classifier 2']+self.secondClass.headerFields+['Start (nm)']+['End (nm)'], dialect='excel-tab')
            self.writer.writeheader()
            self.writer.writerow(secondHeaders)
        else:
            secondHeaders['Classifier 2'] = self.secondClassifierList.currentItem(
            ).text()
            secondHeaders['Data points'] = self.secondClass.points
            secondHeaders['Data length'] = len(self.secondClass.fileFields)+1
            self.writer = csv.DictWriter(self.file, fieldnames=['Classifier 2']+['Data points']+[
                                         'Data length']+self.secondClass.headerFields, dialect='excel-tab')
            self.writer.writeheader()
            self.writer.writerow(secondHeaders)

        cpcHeaders = {'CPC': self.cpcList.currentItem().text()}
        self.writer = csv.DictWriter(self.file, fieldnames=[
                                     'CPC'], dialect='excel-tab')
        self.writer.writeheader()
        self.writer.writerow(cpcHeaders)

        self.fileFields = self.getFileFields()
        self.writer = csv.DictWriter(
            self.file, fieldnames=self.fileFields, dialect='excel-tab')
        self.writer.writeheader()
        self.fileData = {}

    def check(self):  # called every second to check readyness of classifiers
        self.timer.stop()

        if (self.startRow):
            if (self.doBypass or self.firstClass.isReady()):
                self.startRow = False
                self.updateStatus("Delay 1")
                for i in range(int(float(self.firstDelay.toPlainText())*10.0)):
                    time.sleep(0.1)
                    app.processEvents()
                    if (not self.isScanning):
                        return
                # run the scanning 2nd classifier after first classifier is ready
                if (self.isScanner):
                    if (not self.doBypass):
                        if (self.varBins):
                            self.secondClass.VarDiameter(
                                self.firstClass.X[self.firstClass.point])

                        self.updateStatus("Scanning 2")
                        self.secondClass.StartScan()
                        if (self.secondClassifierList.currentItem().text() == "AAC"):
                            if (not self.secondClass.isScanning):
                                QtWidgets.QMessageBox.warning(self, "Tandem Error", "AAC scan error (range?)")
                                self.endScan()
                            self.secondClass.next()
            else:
                self.updateStatus("Stabilising 1")
                if (self.isScanner and self.secondClassifierList.currentItem().text() == "3082 DMA"):
                    self.secondClass.StopScan()
                self.timer.start(1)
                return
        if (not self.startRow):
            # write the measured data to the file and start the next setpoint
            if (self.secondClass.isReady()):
                self.timer.stop()
                if (not (self.isScanner)):
                    self.updateStatus("Delay 2")
                    for i in range(int(float(self.secondDelay.toPlainText())*10.0)):
                        time.sleep(0.1)
                        app.processEvents()
                        if (not self.isScanning):
                            return
                self.next()
                if (not self.isScanning):
                    return
            else:
                self.updateStatus("Stabilising 2")
        self.timer.start(1)

    def next(self):
        self.plotIt()
        if (not self.isScanning):
            return
        if (self.secondClassifierList.currentItem().text() == "3082 DMA" and self.secondScanner.isChecked()):
            if (self.secondClass.DataReady()):  # 2nd class processing have completed
                if (self.doBypass):
                    self.doBypass = False
                if (not self.firstClass.moreToCome()):
                    if (self.finalBypass or not self.Bypass.isChecked()):
                        self.endScan()
                        return
                    else:  # Final bypass starts
                        self.doBypass = True
                        self.finalBypass = True
                        self.disBypass = True
                        # self.firstClass.enableBypass(3)  # When using the analog output of the instrument
                        self.arduino.enableBypass()  # When using the Arduino

                        for i in range(int(float(self.firstDelay.toPlainText())*10.0)):
                            time.sleep(0.1)
                        self.secondClass.StartScan()
                        if (self.secondClassifierList.currentItem().text() == "AAC"):
                            if (not self.secondClass.isScanning):
                                QtWidgets.QMessageBox.warning(self, "Tandem Error", "AAC scan error (range?)")
                                self.endScan()

                self.secondClass.reset()
                self.fileData = {}
                self.startRow = True
                if (not self.doBypass):
                    if (not self.firstClass.next()):
                        self.barf(
                            "Failed to set first classifier. Reason: "+self.firstClass.lastResponse)
                        return

        else:
            if (self.isScanner and self.secondClassifierList.currentItem().text() == "AAC"):
                self.secondClass.next()
            if (not self.secondClass.moreToCome()):  # end of 2nd class sweep, end of file row
                self.Scatterplot = True
                if (self.doBypass):
                    self.doBypass = False
                if (not self.firstClass.moreToCome()):
                    if (self.finalBypass or not self.Bypass.isChecked()):
                        # plot final step
                        self.plot.axes.plot(
                            self.secondClass.X, self.G, label=f"{self.firstClass.label} = {round(self.firstClass.X[self.firstClass.point - 1], 2)} {self.firstFromUnit.text()}", color=self.color)
                        self.plot.axes.legend(
                            loc='upper left', fontsize='small')
                        # Combine legends from both plots and place it outside the plot area
                        legend_lines, legend_labels = self.plot.axes.get_legend_handles_labels()
                        legend_lines2, legend_labels2 = self.plot2.axes.get_legend_handles_labels()
                        combined_legend_lines = legend_lines + legend_lines2
                        combined_legend_labels = legend_labels + legend_labels2
                        # remove the bypass legend
                        if self.plot.get_legend():
                            self.plot.get_legend().remove()
                        self.plot2.axes.legend(
                            combined_legend_lines, combined_legend_labels,
                            loc='upper left', bbox_to_anchor=(1.1, 1),
                            fontsize='small', frameon=False)
                        # self.plot.axes.set_title('Second Classifier Data')
                        # Set non-scientific notation for x-axis tick labels
                        self.plot.axes.xaxis.set_major_formatter(
                            mpl.ticker.ScalarFormatter())
                        self.plot.axes.yaxis.set_major_formatter(
                            mpl.ticker.ScalarFormatter())
                        self.plot.axes.yaxis.set_minor_formatter(
                            mpl.ticker.ScalarFormatter())
                        self.plot2.axes.yaxis.set_major_formatter(
                            mpl.ticker.ScalarFormatter())
                        self.MplWidget.canvas.draw()
                        self.endScan()
                        return
                    else:  # Final bypass starts
                        self.doBypass = True
                        self.finalBypass = True
                        self.disBypass = True
                        # self.firstClass.enableBypass(3)  # When using the Arduino
                        self.arduino.enableBypass()  # When using the Arduino
                        for i in range(int(float(self.firstDelay.toPlainText())*10.0)):
                            time.sleep(0.1)
                        self.secondClass.run()

                self.secondClass.reset()
                self.fileData = {}
                self.startRow = True
                if (not self.doBypass):
                    if (not self.firstClass.next()):
                        self.barf(
                            "Failed to set first classifier. Reason: "+self.firstClass.lastResponse)
                        return
            if (not self.isScanner):
                if (not self.secondClass.next()):
                    self.barf(
                        "Failed to set second classifier. Reason: " + self.secondClass.lastResponse)
                    return

    def plotIt(self):  # do plot, and add data log row
        if (self.isScanner and self.secondClassifierList.currentItem().text() == "AAC"):
            self.G.resize(self.secondClass.points)
            self.X, self.Y = np.meshgrid(
                self.firstClass.X, self.secondClass.X, indexing='ij')
            self.Z.resize((self.firstClass.points, self.secondClass.points))
        ave = int(self.average.toPlainText())
        self.updateStatus("Averaging")
        conc = 0.0
        tzero = time.time()
        # Calculate colour based on first classifer points
        color_value = self.firstClass.point / \
            (self.firstClass.points)  # Range from 0 to 1
        self.color = self.cmap(color_value)

        # Stop the 2nd class # set the pinch valve for the bypass # Start the bypass
        if (self.doBypass and self.firstBypass):
            self.secondClass.StopScan()
            # self.firstClass.enableBypass(3)  # When using the analog output of the instrument
            self.arduino.enableBypass()  # When using Arduino
            for i in range(int(float(self.firstDelay.toPlainText())*10.0)):
                time.sleep(0.1)
            self.secondClass.StartScan()
            if (self.secondClassifierList.currentItem().text() == "AAC"):
                if (not self.secondClass.isScanning):
                    QtWidgets.QMessageBox.warning(self, "Tandem Error", "AAC scan error (range?)")
                    self.endScan()
            self.firstBypass = False
        if (self.isScanner and self.secondClassifierList.currentItem().text() == "AAC"):
            if (not self.doBypass):
                self.firstClass.monitor()
        else:
            self.cpc.startpoll()
            for i in range(ave):
                if (not self.isScanning):
                    return
                if (not self.doBypass):
                    self.firstClass.monitor()  # get feedback from classifier, add to average
                app.processEvents()
                if (not self.isScanning):
                    return
                self.secondClass.monitor()
                app.processEvents()
                conc += self.cpc.conc()  # CPC averaging handled here
                while (time.time()-tzero < 1.0):
                    app.processEvents()
                    if (not self.isScanning):
                        return
                    tzero = time.time()
            self.cpc.endpoll()
            conc /= float(ave)
        # get fed back averages from classifer
        secData = self.secondClass.getFileData()
        for item, data in secData.items():
            self.fileData[item + str(2)] = data

        if (not (self.secondClassifierList.currentItem().text() == "3082 DMA" and self.secondScanner.isChecked())):
            if (self.isScanner):
                conc = self.secondClass.conc
            self.fileData["Conc "] = conc
            self.cpcNumber.display(conc)
            if (self.doBypass):
                bypassDummy = {self.firstClass.fileFields[i]: "Bypassed" for i in range(
                    len(self.firstClass.fileFields))}
                self.fileData = {**self.fileData, **bypassDummy}
                # mapping the CPC data to the list G, and then plot the scatter

                self.G[self.secondClass.point] = conc
                if (not self.finalBypass):
                    # self.plot2.clear()
                    self.plot2.axes.plot(self.secondClass.X, self.G,
                                         label='First Bypass', marker='x', linestyle='-')
                else:
                    self.plot2.axes.plot(self.secondClass.X, self.G,
                                         label='Final Bypass', marker='x', linestyle='-')
                self.plot2.axes.legend(loc='upper right',  fontsize='small')
            else:
                self.fileData = {**self.fileData, **
                                 self.firstClass.getFileData()}
                self.Z[self.firstClass.point, self.secondClass.point] = conc
                # plot a contour
                self.graph.axes.pcolor(
                    self.X, self.Y, self.Z, shading='auto')
                # # Set non-scientific notation for x-axis and y-axis tick labels
                self.graph.axes.xaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.graph.axes.yaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.graph.axes.xaxis.set_minor_formatter(
                    mpl.ticker.ScalarFormatter())
                self.graph.axes.yaxis.set_minor_formatter(
                    mpl.ticker.ScalarFormatter())

                if (self.Scatterplot):
                    # plot a scatter
                    self.plot.axes.plot(self.secondClass.X, self.G,
                                        label=f"{self.firstClass.label} = {round(self.firstClass.X[self.firstClass.point-1], 3)} {self.firstFromUnit.text()}", color=self.color)
                    self.plot.axes.legend(loc='upper left', fontsize='small')
                    # Set non-scientific notation for x-axis tick labels
                    self.plot.axes.xaxis.set_major_formatter(
                        mpl.ticker.ScalarFormatter())
                    self.plot.axes.yaxis.set_major_formatter(
                        mpl.ticker.ScalarFormatter())
                    self.plot.axes.xaxis.set_minor_formatter(
                        mpl.ticker.ScalarFormatter())
                    self.plot.axes.yaxis.set_minor_formatter(
                        mpl.ticker.ScalarFormatter())
                    # Combine legends from both plots and place it outside the plot area
                    legend_lines, legend_labels = self.plot.axes.get_legend_handles_labels()
                    legend_lines2, legend_labels2 = self.plot2.axes.get_legend_handles_labels()
                    combined_legend_lines = legend_lines + legend_lines2
                    combined_legend_labels = legend_labels + legend_labels2
                    # remove the bypass legend
                    if self.plot.get_legend():
                        self.plot.get_legend().remove()
                    self.plot2.axes.legend(
                        combined_legend_lines, combined_legend_labels,
                        loc='upper left', bbox_to_anchor=(1.04, 1),
                        fontsize='small', frameon=False)
                    # self.plot.axes.set_title('Second Classifier Data')
                    self.MplWidget.canvas.draw()

                    # Set non-scientific notation for x-axis tick labels
                    self.plot_live.axes.xaxis.set_major_formatter(
                        ScalarFormatter())
                    # Adjust the layout to accommodate the legend outside the box
                    self.MplWidget.canvas.figure.tight_layout(
                        rect=[0, 0, 1, 1])   # Adjust the right margin
                    self.Scatterplot = False
                    self.G = np.zeros(self.secondClass.points)

                # mapping the CPC data to the list G, and then plot the  scatter
                self.G[self.secondClass.point] = conc
                # Clear the existing live-plot
                self.plot_live.clear()
                # plot a live scatter
                self.plot_live.set_xscale('log')
                self.plot_live.set_xlabel(
                    f"{self.secondClass.label}[{self.secondFromUnit.text()}]")
                self.plot_live.set_ylabel("N [Counts]")
                self.plot_live.axes.plot(self.secondClass.X, self.G,
                                         label=f"{self.firstClass.label} = {round(self.firstClass.x, 2)} {self.firstFromUnit.text()}")
                self.plot_live.axes.legend(loc='upper left', fontsize='small')
                # Set non-scientific notation for x-axis tick labels
                self.plot_live.axes.xaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot_live.axes.xaxis.set_minor_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot_live.axes.yaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot_live.axes.yaxis.set_minor_formatter(
                    mpl.ticker.ScalarFormatter())
                # Draw the plot
                self.plotWidget.canvas.draw_idle()
                # Adjust the layout to prevent cutoff of axis labels
                self.plotWidget.canvas.figure.tight_layout()

            # Adjust the layout to accommodate the legend outside the box
            self.plotWidget.canvas.figure.tight_layout(
                rect=[0, 0, 1, 1])  # Adjust the right margin
            # self.plot_live.axes.set_title('Second Classifier Live Data')
            self.MplWidget_live.canvas.draw()
            # Adjust the layout to prevent cutoff of axis labels
            self.MplWidget_live.canvas.figure.tight_layout()

            self.writer.writerow(self.fileData)  # commit log row to file

        # If the 2nd class is DMA-3082 and the scanning is done, commit log row to file
        elif (self.secondClassifierList.currentItem().text() == "3082 DMA" and self.secondScanner.isChecked()):
            # Get rid of Scientific notation
            self.plot.axes.xaxis.set_major_formatter(
                mpl.ticker.ScalarFormatter())
            self.plot.axes.yaxis.set_major_formatter(
                mpl.ticker.ScalarFormatter())
            self.plot.axes.yaxis.set_minor_formatter(
                mpl.ticker.ScalarFormatter())
            self.plot2.axes.yaxis.set_major_formatter(
                mpl.ticker.ScalarFormatter())
            if (self.secondClass.DataReady()):
                LowerSizeRange = self.secondClass.LowerSizeRange()
                UpperSizeRange = self.secondClass.UpperSizeRange()
                Conc_3082 = self.secondClass.output3082()
                conc3082 = {}
                # set the concentration of 3082,
                conc3082 = {i: Conc_3082[i] for i in range(193)}
                for i in range(int(LowerSizeRange), int(UpperSizeRange)):
                    self.fileData["Dm (nm)2"] = self.secondLowerRange.itemText(
                        i)   # read the textvalue of particle diameter
                    self.fileData["Conc "] = conc3082[i+1]
                    if (self.doBypass):
                        bypassDummy = {self.firstClass.fileFields[i]: "Bypassed" for i in range(
                            len(self.firstClass.fileFields))}
                        self.fileData = {**self.fileData, **bypassDummy}
                        if (self.disBypass):  # return the pinch valve to  normal scanning
                            # self.firstClass.disableBypass(3)  # When using the analog output of the instrument

                            self.arduino.disableBypass()  # When using Arduino
                            for i in range(int(float(self.firstDelay.toPlainText())*10.0)):
                                time.sleep(0.1)
                            self.disBypass = False
                    else:
                        self.fileData = {**self.fileData, **
                                         self.firstClass.getFileData()}
                    # commit log row to file
                    self.writer.writerow(self.fileData)

                if self.makesecondClassX:  # making 3082 DMA points as an array
                    self.secondClassX = np.array([float(self.secondLowerRange.itemText(
                        i)) for i in range(int(LowerSizeRange), int(UpperSizeRange))])
                    self.Z = np.zeros(
                        (self.firstClass.points, self.secondClassX.shape[0]))
                    self.z = np.zeros(self.secondClassX.shape[0])
                    # self.makesecondClassX = False
                # self.z = np.zeros(self.secondClassX.shape[0])
                self.X, self.Y = np.meshgrid(
                    self.firstClass.X, self.secondClassX, indexing='ij')

                # mapping the CPC data to the matrix Z and list z, and then plot the contour and scatter, respectively
                j = 0

                for i in range(int(LowerSizeRange), int(UpperSizeRange)):
                    self.Z[self.firstClass.point, j] = conc3082[i+1]
                    self.z[j] = conc3082[i+1]
                    j += 1
                if (not self.doBypass):
                    # plot a contour
                    self.graph.axes.pcolor(
                        self.X, self.Y, self.Z, shading='auto')

                    # plot a scatter
                    self.plot.axes.plot(self.secondClassX, self.z,
                                        label=f"{self.firstClass.label} = {round(self.firstClass.x, 3)} {self.firstFromUnit.text()}", color=self.color)
                    self.plot.axes.legend(loc='upper left', fontsize='small')

                    # Clear the existing live-plot
                    self.plot_live.clear()
                    self.plot_live.set_xscale('log')
                    # plot a live scatter
                    self.plot_live.set_xlabel(
                        f"{self.secondClass.label}[{self.secondFromUnit.text()}]")
                    self.plot_live.set_ylabel("N [Counts]")
                    self.plot_live.axes.plot(self.secondClassX, self.z,
                                             label=f"{self.firstClass.label} = {round(self.firstClass.x, 3)} {self.firstFromUnit.text()}")
                    self.plot_live.axes.legend(
                        loc='upper left', fontsize='small')
                    # Draw the plot
                    self.plotWidget.canvas.draw_idle()
                    # Adjust the layout to prevent cutoff of axis labels
                    self.plotWidget.canvas.figure.tight_layout()

                else:  # plot ByPass results
                    if (not self.finalBypass):
                        # self.plot2.clear()
                        self.Z = np.zeros(
                            (self.firstClass.points, self.secondClassX.shape[0]))
                        self.plot2.axes.plot(self.secondClassX, self.z,
                                             label='First Bypass', marker='x', linestyle='-')
                    else:
                        self.plot2.axes.plot(self.secondClassX, self.z,
                                             label='Final Bypass', marker='x', linestyle='-')
                    self.plot2.axes.legend(
                        loc='upper right',  fontsize='small')

                # Set non-scientific notation for x-axis tick labels
                self.plot.axes.xaxis.set_major_formatter(ScalarFormatter())
                self.plot_live.axes.xaxis.set_major_formatter(
                    ScalarFormatter())

                # Get rid of Scientific notation

                # Combine legends from both plots and place it outside the plot area
                legend_lines, legend_labels = self.plot.axes.get_legend_handles_labels()
                legend_lines2, legend_labels2 = self.plot2.axes.get_legend_handles_labels()

                combined_legend_lines = legend_lines + legend_lines2
                combined_legend_labels = legend_labels + legend_labels2
                # remove the bypass legend
                if self.plot.get_legend():
                    self.plot.get_legend().remove()
                self.plot2.axes.legend(
                    combined_legend_lines, combined_legend_labels,
                    loc='upper left', bbox_to_anchor=(1.1, 1),
                    fontsize='small', frameon=False)

                # Get rid of Scientific notation
                self.plot.axes.xaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot.axes.yaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot.axes.yaxis.set_minor_formatter(
                    mpl.ticker.ScalarFormatter())
                self.plot2.axes.yaxis.set_major_formatter(
                    mpl.ticker.ScalarFormatter())
                self.MplWidget.canvas.draw()

                # Adjust the layout to accommodate the legend outside the box
                self.plotWidget.canvas.figure.tight_layout()  # Adjust the right margin
                self.MplWidget.canvas.figure.tight_layout()   # Adjust the right margin

                # self.plot.axes.set_title('Second Classifier Data')
                # self.plot_live.axes.set_title('Second Classifier Live Data')

                self.MplWidget_live.canvas.draw()
                # Adjust the layout to prevent cutoff of axis labels
                self.MplWidget_live.canvas.figure.tight_layout()

        self.updateStatus()

    def updateStatus(self, message=None):
        if (message != None):
            message = " : "+message
        else:
            message = ""
        if (self.isScanner):
            if (self.doBypass):
                self.statusbar.showMessage("Bypass scan: ")
            else:
                self.statusbar.showMessage(
                    "Scanning: "+str(self.firstClass.point+1)+"/"+str(self.firstClass.points)+message)
        else:
            if (self.doBypass):
                self.statusbar.showMessage(
                    "Bypass scan: "+str(self.secondClass.point+1)+"/"+str(self.secondClass.points)+message)
            else:
                self.statusbar.showMessage("Scanning: "+str(self.firstClass.point+1)+"/"+str(
                    self.firstClass.points)+" & "+str(self.secondClass.point+1)+"/"+str(self.secondClass.points)+message)
        app.processEvents()

    def barf(self, message):
        self.endScan()
        QtWidgets.QMessageBox.warning(self, "Tandem Error", message)

    def endScan(self):
        self.statusbar.clearMessage()
        self.stopButton.setEnabled(False)
        self.startButton.setEnabled(True)
        if (not self.isScanning):
            return
        try:
            self.file.close()
        except:
            pass
        self.isScanning = False
        self.timer.stop()
        app.processEvents()
        time.sleep(1)
        del self.firstClass
        del self.secondClass
        del self.cpc

    def firstClassifierChanged(self, c):
        if (c == "CPMA"):
            self.firstFromUnit.setText("fg")
            self.firstToUnit.setText("fg")
            self.firstRorQshLabel.setText("Rm")
            self.firstRorQshUnit.setText("")

        else:
            self.firstFromUnit.setText("nm")
            self.firstToUnit.setText("nm")
            self.firstRorQshLabel.setText("Sheath")
            self.firstRorQshUnit.setText("lpm")

        if (c == "AAC" or c == "CPMA"):
            self.firstIsSerial.setEnabled(True)
            self.firstIsEth.setEnabled(True)

        if (c == "3081 DMA"):
            self.firstIsSerial.setEnabled(True)
            self.firstIsSerial.setChecked(True)
            self.firstIsEth.setDisabled(True)

        if (c == "3082 DMA"):
            self.firstIsEth.setEnabled(True)
            self.firstIsEth.setChecked(True)
            self.firstIsSerial.setDisabled(True)
            self.firstPolarity.setEnabled(True)

        else:
            self.firstPolarity.setDisabled(True)

    def secondClassifierChanged(self, c):
        if (c == "CPMA"):
            self.secondFromUnit.setText("fg")
            self.secondToUnit.setText("fg")
            self.secondRorQshLabel.setText("Rm")
            self.secondRorQshUnit.setText("")
            self.secondPerDecade.setEnabled(True)
            self.secondPerDecade.setDisabled(True)
            self.secondScanner.setChecked(False)
            self.secondScanner.setEnabled(False)
        else:
            self.secondFromUnit.setText("nm")
            self.secondToUnit.setText("nm")
            self.secondRorQshLabel.setText("Sheath")
            self.secondRorQshUnit.setText("lpm")

        if (c == "AAC" or c == "CPMA"):
            self.secondIsSerial.setEnabled(True)
            self.secondIsEth.setEnabled(True)

        if (c == "3081 DMA"):
            self.secondIsSerial.setEnabled(True)
            self.secondIsSerial.setChecked(True)
            self.secondIsEth.setDisabled(True)
            self.secondScanner.setChecked(False)
            self.secondScanner.setEnabled(False)

        if (c == "3082 DMA"):
            self.secondIsEth.setEnabled(True)
            self.secondIsEth.setChecked(True)
            self.secondIsSerial.setDisabled(True)
            self.secondScanner.setEnabled(True)
            self.secondPolarity.setEnabled(True)

        else:
            self.secondPolarity.setDisabled(True)

        if (c == "AAC"):
            self.secondScanner.setEnabled(True)
            self.secondFrom.setEnabled(True)
            self.secondTo.setEnabled(True)

        if (self.secondScanner.isChecked()):
            if (c == "AAC"):
                self.secondFrom.setEnabled(False)
                self.secondTo.setEnabled(False)
            else:
                self.secondFrom.setEnabled(True)
                self.secondTo.setEnabled(True)

            if (c == "AAC" or c == "3082 DMA"):
                self.variableBins.setEnabled(True)

            else:
                self.variableBins.setChecked(False)
                self.variableBins.setDisabled(True)

            if (c == "3082 DMA"):
                self.secondHighFlow.setDisabled(False)
                self.secondLowFlow.setDisabled(False)
                self.secondScanUpTime.setEnabled(True)
                self.SutUnit.setEnabled(True)
                self.SutText.setEnabled(True)
                self.secondPerDecade.setEnabled(False)
                self.secondClassText.setEnabled(False)
                self.secondClassUnit.setEnabled(False)
                self.secondFrom.setVisible(False)
                self.secondTo.setVisible(False)
                self.secondLowerRange.setVisible(True)
                self.secondUpperRange.setVisible(True)
                self.scanDelay.setEnabled(False)
                self.scanAve.setEnabled(False)
                self.secondDelay.setEnabled(False)
                self.WaterPF_lable.setVisible(False)
                self.WaterMME_lable.setVisible(False)
                self.SootPF_lable.setVisible(False)
                self.SootMME_lable.setVisible(False)

            else:
                self.secondHighFlow.setDisabled(True)
                self.secondLowFlow.setDisabled(True)
                self.secondScanUpTime.setEnabled(True)
                self.SutUnit.setEnabled(True)
                self.SutText.setEnabled(True)
                self.secondPerDecade.setEnabled(True)
                self.secondClassText.setEnabled(True)
                self.secondClassUnit.setEnabled(True)
                self.secondFrom.setVisible(True)
                self.secondTo.setVisible(True)
                self.secondLowerRange.setVisible(False)
                self.secondUpperRange.setVisible(False)
                self.scanDelay.setEnabled(True)
                self.scanAve.setEnabled(True)
                self.secondDelay.setEnabled(False)
                self.mobilityDiameterBox.setEnabled(False)
                self.WaterPF_lable.setVisible(False)
                self.WaterMME_lable.setVisible(False)
                self.SootPF_lable.setVisible(False)
                self.SootMME_lable.setVisible(False)

        else:
            self.variableBins.setChecked(False)
            self.SecondBinsChanged(False)
            self.secondHighFlow.setDisabled(True)
            self.secondLowFlow.setDisabled(True)
            self.secondScanUpTime.setDisabled(True)
            self.SutUnit.setDisabled(True)
            self.SutText.setDisabled(True)
            self.secondPerDecade.setEnabled(True)
            self.secondClassText.setEnabled(True)
            self.secondClassUnit.setEnabled(True)
            self.secondLowerRange.setVisible(False)
            self.secondUpperRange.setVisible(False)
            self.scanDelay.setEnabled(False)
            self.scanAve.setEnabled(False)
            self.secondDelay.setEnabled(True)
            self.aveText.setEnabled(False)
            self.secondDelayText.setEnabled(True)
            self.scanDelayText.setEnabled(False)
            self.scanDelayUnit.setEnabled(False)
            self.variableBins.setChecked(False)
            self.variableBins.setDisabled(True)
            self.mobilityDiameterBox.setEnabled(False)
            self.WaterPF_lable.setVisible(False)
            self.WaterMME_lable.setVisible(False)
            self.SootPF_lable.setVisible(False)
            self.SootMME_lable.setVisible(False)

    def cpcChanged(self, c):
        if (c == "3022/25" or c == "3775/76" or c == "Magic" or c=="375x"):
            self.cpcIP.setEnabled(False)
            self.cpcPort.setEnabled(True)
            self.cpcIsSerial.setEnabled(True)
            self.cpcIsSerial.setChecked(True)
            self.cpcIsEth.setDisabled(True)

        elif (c == "Cambustion 5210"):
            self.cpcIsSerial.setEnabled(True)
            self.cpcIsEth.setEnabled(True)
            self.cpcIP.setEnabled(self.cpcIsEth.isChecked())
            self.cpcPort.setEnabled(self.cpcIsSerial.isChecked())

    def secondScanChanged(self, c):

        if (c):
            self.secondHighFlow.setDisabled(True)
            self.secondLowFlow.setDisabled(True)
            self.secondPerDecade.setEnabled(False)
            self.secondClassText.setEnabled(False)
            self.secondClassUnit.setEnabled(False)
            self.secondDelay.setEnabled(False)
            self.secondDelayText.setEnabled(False)
            self.secondDelayUnit.setEnabled(False)
            self.secondScanUpTime.setEnabled(True)
            self.SutUnit.setEnabled(True)
            self.SutText.setEnabled(True)
            self.variableBins.setChecked(False)
            self.variableBins.setDisabled(True)
            self.mobilityDiameterBox.setEnabled(False)
            self.secondFrom.setVisible(True)
            self.secondTo.setVisible(True)

            if (self.secondClassifierList.currentItem().text() == "3082 DMA" or self.secondClassifierList.currentItem().text() == "AAC"):
                self.variableBins.setEnabled(True)
                self.secondFrom.setDisabled(False)
                self.secondTo.setDisabled(False)
            else:
                self.variableBins.setChecked(False)
                self.variableBins.setDisabled(True)
                self.secondFrom.setDisabled(True)
                self.secondTo.setDisabled(True)

            if (self.secondClassifierList.currentItem().text() == "3082 DMA"):
                self.secondHighFlow.setDisabled(False)
                self.secondLowFlow.setDisabled(False)
                self.secondFrom.setVisible(False)
                self.secondTo.setVisible(False)
                self.secondLowerRange.setVisible(True)
                self.secondUpperRange.setVisible(True)
                self.scanDelay.setEnabled(False)
                self.scanDelayText.setEnabled(False)
                self.scanDelayUnit.setEnabled(False)
                self.scanAve.setEnabled(False)
                self.aveText.setEnabled(False)

            else:
                self.secondHighFlow.setDisabled(True)
                self.secondLowFlow.setDisabled(True)
                self.secondLowerRange.setVisible(False)
                self.secondUpperRange.setVisible(False)
                self.scanDelay.setEnabled(True)
                self.scanAve.setEnabled(True)
                self.scanDelayText.setEnabled(True)
                self.scanDelayUnit.setEnabled(True)
                self.aveText.setEnabled(True)
                self.mobilityDiameterBox.setEnabled(False)

        else:
            self.secondHighFlow.setDisabled(True)
            self.secondLowFlow.setDisabled(True)
            self.secondPerDecade.setEnabled(True)
            self.secondClassText.setEnabled(True)
            self.secondClassUnit.setEnabled(True)
            self.secondFrom.setVisible(True)
            self.secondTo.setVisible(True)
            self.secondLowerRange.setVisible(False)
            self.secondUpperRange.setVisible(False)
            self.scanDelay.setEnabled(False)
            self.scanDelayText.setEnabled(False)
            self.scanDelayUnit.setEnabled(False)
            self.scanAve.setEnabled(False)
            self.aveText.setEnabled(False)
            self.secondDelayText.setEnabled(True)
            self.secondDelay.setEnabled(True)
            self.secondDelayUnit.setEnabled(True)
            self.SutUnit.setEnabled(False)
            self.SutText.setEnabled(False)
            self.variableBins.setChecked(False)
            self.variableBins.setDisabled(True)
            self.mobilityDiameterBox.setEnabled(False)

    def SecondBinsChanged(self, c):
        if (c):
            self.mobilityDiameterBox.setEnabled(True)
            self.secondFrom.setEnabled(True)
            self.secondTo.setEnabled(True)

            if (self.secondClassifierList.currentItem().text() == "3082 DMA" or self.secondClassifierList.currentItem().text() == "AAC"):
                self.secondFrom.setDisabled(True)
                self.secondTo.setDisabled(True)
                self.secondLowerRange.setDisabled(True)
                self.secondUpperRange.setDisabled(True)

            if self.particleIsWater.isChecked():
                self.WaterPF_lable.setVisible(True)
                self.WaterMME_lable.setVisible(True)
                self.PreFactor.setVisible(False)
                self.MM_Exponent.setVisible(False)

            else:
                self.WaterPF_lable.setVisible(False)
                self.WaterMME_lable.setVisible(False)

            if self.particleIsSoot.isChecked():
                self.SootPF_lable.setVisible(True)
                self.SootMME_lable.setVisible(True)
                self.PreFactor.setVisible(False)
                self.MM_Exponent.setVisible(False)

            else:
                self.SootPF_lable.setVisible(False)
                self.SootMME_lable.setVisible(False)

            if self.particleIsNone.isChecked():
                self.PreFactor.setVisible(True)
                self.MM_Exponent.setVisible(True)

        else:
            self.mobilityDiameterBox.setEnabled(False)
            self.secondFrom.setEnabled(True)
            self.secondTo.setEnabled(True)
            if (self.secondClassifierList.currentItem().text() == "3082 DMA"):
                self.secondLowerRange.setDisabled(False)
                self.secondUpperRange.setDisabled(False)

    def selectFile(self):
        self.fileName, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Log to", "", "Text Files (*.txt)")
        self.rawFileDisplay.setText("Raw Data File: "+self.fileName)

    def getFileFields(self):
        secFields = []
        secFields = [field + str(2)
                     for field in self.secondClass.fileFields]+["Conc "]
        # note we _prepend_ the firstClass headers. This will set the order in the file.
        return self.firstClass.fileFields+secFields

    def testCPC(self):
        self.statusbar.showMessage("Initialising CPC...")
        c = self.cpcList.currentItem().text()
        if (c == "Cambustion 5210"):
            self.cpc = instruments.CambustionCPC
        elif (c == "3022/25"):
            self.cpc = instruments.TSI30xx
        elif (c == "3775/76"):
            self.cpc = instruments.TSI377x
        elif (c == "375x"):
            self.cpc = instruments.TSI375x
        elif (c == "Dummy"):
            self.cpc = instruments.DummyCPC
        elif (c == "Magic"):
            self.cpc = instruments.MagicCpc

        self.cpc = self.cpc(self.cpcIsSerial.isChecked(), self.cpcIP.toPlainText(), self.cpcPort.toPlainText())

        if (not self.cpc.connected):
            self.barf("Connection Error to CPC")
            del self.cpc
            return

        self.cpcNumber.display(self.cpc.conc())

        self.cpc.disconnect()

    def invert(self):
        self.barf("Inversion feature not implemented yet.")
        

app.setApplicationName("Tandem")
app.setOrganizationName("Cambustion and The University of Alberta")
window = Ui()
app.exec_()
