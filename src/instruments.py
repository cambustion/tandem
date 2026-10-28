# Defines all the instruments (classifiers & CPCs)
# Python application for automation of tandem aerosol classifier experiments and data inversion
# Copyright (C) 2024 Cambustion, The University of Alberta, NRC Canada
# Original Authors: Jonathan Symonds, Morteza Kiasadegh, Tim Sipkens
# See LICENSE for details

import math
import socket
import sys
import time
import serial
import platform
import json
import random
import numpy as np


class Instrument:  # Everything is an instrument

    def __init__(self, isSerial, ip, serPort):
        self.isSerial = isSerial
        self.ip = ip
        self.serPort = serPort
        self.ser = serial.Serial()

        self.ser.timeout = 1
        self.connected = self.connect()
        self.lastResponse = ""

    def __del__(self):
        self.disconnect()
        self.connected = False

    def connect(self):
        if (self.isSerial):
            return self.connectSerial()
        else:
            return self.connectEth()

    def disconnect(self):
        if (self.isSerial):
            return self.disconnectSerial()
        else:
            return self.disconnectEth()

    def connectSerial(self):
        self.ser.baudrate = self.baudrate
        self.ser.bytesize = self.databits
        self.ser.stopbits = self.stopbits
        self.ser.parity = self.parity
        if (self.serPort.isnumeric()):
            if (platform.system() == 'Windows'):
                self.ser.port = "COM"+str(self.serPort)  # for PC
            else:
                # for WSL and some UNIX - note WSL2 does not yet support serial ports!
                self.ser.port = "/dev/ttyS"+str(self.serPort)
        else:
            self.ser.port = self.serPort  # for everything else
        try:
            self.ser.open()
        except:
            return False
        return self.ser.is_open

    def connectEth(self):
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            server_address = (self.ip, self.ipPort)
            self.sock.connect(server_address)
        except:
            return False
        return True

    def disconnectSerial(self):
        try:
            self.ser.close()
        except:
            pass

    def disconnectEth(self):
        try:
            self.sock.close()
        except:
            pass

    def sendCmd(self, cmd):
        if (self.isSerial):
            return self.sendCmdSerial(cmd + "\r\n")
        else:
            return self.sendCmdEth(cmd + "\r\n")

    def getRsp(self):
        if (self.isSerial):
            return self.getRspSerial()
        else:
            return self.getRspEth()

    def getRspEth(self):
        try:
            return self.sock.recv(4096).decode()
        except:
            self.connected = False
            return ""

    def getRspSerial(self):
        try:
            return self.ser.readline().decode()
        except:
            self.connected = False
            return ""

    def sendCmdSerial(self, cmd):
        try:
            self.ser.write(cmd.encode())
        except:
            self.connected = False

    def sendCmdEth(self, cmd):
        try:
            self.sock.sendall(cmd.encode())
        except:
            self.connected = False

    def sendQuery(self, cmd):
        self.sendCmd(cmd)
        time.sleep(self.querydelay)
        return self.getRsp()

    def getFloat(self, cmd):
        return float(self.sendQuery(cmd))

    def match(self, s, v):
        return (s[0:len(v)] == v)


class Classifier(Instrument):

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh):
        self.fileFields = [self.quantity]+self.valueFields
        self.start = start
        self.end = end
        self.perdec = perdec
        self.points = math.floor(
            ((perdec)*(math.log10(end)-math.log10(start)))+0.00001)+1
        self.point = -1
        self.flow = flow
        self.RorQsh = RorQsh

        self.x = 0.0
        self.X = np.zeros(self.points)
        for i in range(self.points):
            self.X[i] = self.getPoint(i)
        super().__init__(isSerial, ip, serPort)
        self.setRorQsh(RorQsh)

    def next(self):
        self.values = {self.valueFields[i]: 0.0 for i in range(
            len(self.valueFields))}
        self.tally = 0
        self.point += 1
        self.x = self.X[self.point]
        return self.set(self.x)

    def moreToCome(self):
        return (self.point+1 < self.points)

    def reset(self):
        self.point = -1

    def getPoint(self, p):
        return self.start*((10.0**(1.0/self.perdec))**p)

    def monitor(self):  # gets the feedback values which we store in the file
        r = self.doMonitorCmd()
        self.values = {self.valueFields[i]: (self.values[self.valueFields[i]]+float(
            r[self.valueFields[i]])) for i in range(len(self.valueFields))}
        self.tally += 1

    def getFileData(self):
        return {**{self.quantity: self.x}, **{self.valueFields[i]: (self.values[self.valueFields[i]]/float(self.tally)) for i in range(len(self.valueFields))}}

    def bypassDummy(self):
        return {self.fileFields[i]: "Bypassed" for i in range(len(self.fileFields))}

    def doMonitorCmd(self):  # sends the command which gets the feedback values
        r = self.sendQuery(self.monitorCmd).split(self.monitorSep)
        return {self.monitorFields[i]: r[i] for i in range(len(r))}

    def run(self):  # default for classifiers which don't "run or stop"
        return True

    def stop(self):
        return True

# AACs and CPMAs are Cambustion classifiers - quite a lot of commonality


class CambustionClass(Classifier):
    ipPort = 23
    timeout = 3.0
    querydelay = 1.8
    baudrate = 19200
    databits = 8
    stopbits = 1
    parity = 'N'
    monitorCmd = "monitor"  # command to get feedback values
    monitorSep = " "  # separator for data returned
    bypassVoltage = 5.0
    noBypassVoltage = 0.0
    bypassChannel = 3

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh):
        super().__init__(isSerial, ip, serPort, start, end, perdec, flow, RorQsh)
        self.setFlow(flow)
        self.set(self.start)

    def __del__(self):
        self.stop()
        super().__del__()

    def connect(self):
        if (not super().connect()):
            return False
        if (self.isSerial):
            self.ser.write(bytes([4]))
        return True

    def getRsp(self):
        chunks = ""
        start = time.time()
        while (chunks.find("\r\n>") < 0 and time.time()-start < self.timeout):
            try:
                if (self.isSerial):
                    chunk = self.ser.read().decode()
                else:
                    chunk = self.sock.recv(1024).decode()
            except:
                self.connected = False
                return ""
            chunks += chunk
        return chunks[0:chunks.find("\r\n>")].replace(">", "")

    def sendFloat(self, cmd, val):
        return (self.sendAndCheck(cmd + " " + "%0.4E" % val))

    def sendAndCheck(self, cmd):
        self.lastResponse = self.sendQuery(cmd)
        return (self.match(self.lastResponse, "OK"))

    def setFlow(self, flow):
        self.sampleFlow = flow
        return self.sendFloat("SetSampleFlow", flow)

    def run(self):
        #          return True # uncomment to use "hardwareless" CPMA/AAC screen
        return self.sendAndCheck("start")

    def stop(self):
        return self.sendAndCheck("stop")

    def isReady(self):
        #          return True # uncomment to use "hardwareless" CPMA/AAC screen
        return (self.match(self.sendQuery("Status"), "Running"))

    def disconnectEth(self):
        try:
            self.sock.send(bytes([13, 10, 4]))  # ctrl-D
            time.sleep(0.05)
            super().disconnectEth()
        except:
            pass

    def enableBypass(self, bypassChannel):
        # set analogue out to "fixed" mode
        self.sendAndCheck("SetAOFunc "+str(bypassChannel)+" 1")
        return self.sendAndCheck("SetAOV "+str(bypassChannel)+" "+"%0.4E" % self.bypassVoltage)

    def disableBypass(self, bypassChannel):
        self.sendAndCheck("SetAOFunc "+str(bypassChannel)+" 1")
        return self.sendAndCheck("SetAOV "+str(bypassChannel)+" "+"%0.4E" % self.noBypassVoltage)


class CPMA(CambustionClass):
    monitorFields = ["Speed (rad/s)", "Voltage (V)", "HT current", "brake temp", "motor current", "power stage temp", "ref pressure", "Temperature (C)", "un-scaled voltage", "voltage PID demand", "voltage gain range", "door locked", "door closed", "detector 1", "detector 2", "serial detector", "AI1 V",
                     "AI2 V", "AI3 V", "AO1 V", "AO2 V", "AO3 V", "run status", "motor temperature", "output function", "PSU voltage", "Pressure (Pa)", "outer speed feedback", "inner speed feedback", "vibration level", "sensor board ref voltage", "lock actuator position"]  # fields returned by monitor command
    # which of those fields we actually want to record in our file
    valueFields = ["Speed (rad/s)", "Voltage (V)",
                   "Pressure (Pa)", "Temperature (C)"]
    quantity = "m* (fg)"
    label = "m*"

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh):
        super().__init__(isSerial, ip, serPort, start, end, perdec, flow, RorQsh)
        self.calK = 0.0  # Default fallback if version is < 2.410 or query fails


    def connect(self):
        if (not super().connect()):
            return False
        
        rsp = self.getRsp()
        if not self.match(rsp, "Cambustion CPMA"):
            return False

        # Check Version and Retrieve K-Factor
        self.update_k_factor()
        return True

    def parse_version(self, ver_str):
        """Extracts numerical version (e.g., '2.410') from CPMA response string."""
        import re
        match = re.search(r'(\d+\.\d+)', ver_str)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                return 0.0
        return 0.0

    def update_k_factor(self):
        """Checks if CPMA UI >= 2.410 and fetches active K-factor."""
        ver_rsp = self.sendQuery("Version")
        ver_num = self.parse_version(ver_rsp)

        if ver_num >= 2.410:
            which_cal = self.sendQuery("GetWhichCal").strip().lower()
            if "user" in which_cal:
                k_rsp = self.sendQuery("GetUserCal")
            else:
                k_rsp = self.sendQuery("GetFactoryCal")
            
            try:
                self.calK = float(k_rsp.strip())
            except ValueError:
                self.calK = 0.0
        else:
            self.calK = 0.0
    
    def GetRealMass(self, V, omega):
        """
        Calculates actual mass in fg using the K-factor adjusted CPMA geometry.
        Formula matches CPMA UI v2.410+ specification.
        """
        if omega == 0.0:
            return 0.0

        e = 1.6e-19
        r_inner_nominal = 60.0e-3
        r_outer_nominal = 61.0e-3

        gap_nominal = r_outer_nominal - r_inner_nominal
        rInner = r_inner_nominal + (self.calK * gap_nominal) / 2.0
        rOuter = r_outer_nominal - (self.calK * gap_nominal) / 2.0
        rAve = (rOuter + rInner) / 2.0

        # Output in fg (1.0E18 conversion multiplier)
        mass_fg = (1.0e18 * e * V) / (omega * omega * rAve * rAve * np.log(rOuter / rInner))
        return mass_fg

    def getFileData(self):
        # Get the standard dictionary of averaged values from the parent class
        data = super().getFileData()
        
        # Extract the dynamically averaged speed and voltage for this specific setpoint
        omega = data["Speed (rad/s)"]
        V = data["Voltage (V)"]
        
        # Calculate and append the actual mass to the data dictionary
        data[self.quantity] = self.GetRealMass(V, omega)
        
        return data
    
    def set(self, x):
        return self.sendFloat("SetMass", x)

    def setRorQsh(self, Rm):
        self.Rm = Rm
        return self.sendFloat("SetRm", Rm)

    def getHeader(self):  # information to put in header row for the classifier
        header = {}
        header['Serial number'] = self.sendQuery("serial")
        header['Sample flow (lpm)'] = self.sampleFlow
        header['Resolution'] = self.Rm
        header['Start (fg)'] = self.start
        header['End (fg)'] = self.end
        header['Per decade'] = self.perdec
        header['Calibration K Factor'] = self.calK
        self.headerFields = list(header.keys())
        return header
    


class AAC(CambustionClass):
    monitorFields = ["Speed (rad/s)", "Sheath flow (lpm)", "brake temp", "motor current", "power stage temp", "ref pressure", "Temperature (C)", "door locked", "door closed", "analogue detector",
                     "serial detector", "AI1 V", "AI2 V", "AI3 V", "AO1 V", "AO2 V", "AO3 V", "run status", "motor temperature", "main CPC", "PSU voltage", "Pressure (Pa)", "classifier speed feedback"]
    valueFields = ["Speed (rad/s)", "Sheath flow (lpm)",
                   "Pressure (Pa)", "Temperature (C)"]
    quantity = "da* (nm)"
    label = "da*"
    isScanning = False
    mtc = False

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh, ScanUpTime, resAve, delayTime, PreFactor, Exponent, f_lower, f_upper, IsWater, IsSoot, IsNone, variableBins, isScanner):
        # Create an instance-specific copy of the fields so we don't alter the global class
        self.valueFields = list(self.valueFields)
        
        super().__init__(isSerial, ip, serPort, start, end, perdec, flow, RorQsh)
        self.isScanner = isScanner
        
        if (self.isScanner):
            # Dynamically remove the Speed column from the file output arrays
            if "Speed (rad/s)" in self.valueFields:
                self.valueFields.remove("Speed (rad/s)")
            if "Speed (rad/s)" in self.fileFields:
                self.fileFields.remove("Speed (rad/s)")
                
            self.resAve = resAve
            self.delayTime = delayTime
            self.ScanUpTime = ScanUpTime
            self.X = np.zeros(0)
            self.points = 0
            self.conc = 0.0
        self.variableBins = variableBins
        self.IsWater = IsWater
        self.IsSoot = IsSoot
        self.IsNone = IsNone
        if f_lower and f_lower.strip():  # Check if not None and not an empty string
            self.f_lower = float(f_lower)
        else:
            self.f_lower = 0.85
        if f_upper and f_upper.strip():  # Check if not None and not an empty string
            self.f_upper = float(f_upper)
        else:
            self.f_upper = 1.1

        if self.IsWater:
            self.PreFactor = 523.6
            self.Exponent = 3
        elif self.IsSoot:
            self.PreFactor = 0.0612
            self.Exponent = 2.48
        elif self.IsNone:
            if PreFactor and PreFactor.strip():  # Check if not None and not an empty string
                self.PreFactor = float(PreFactor)
            else:
                self.PreFactor = 523.6

            if Exponent and Exponent.strip():  # Check if not None and not an empty string
                self.Exponent = float(Exponent)
            else:
                self.Exponent = 3

    def VarDiameter(self, mStar):  # m = K * dm ^ Dm based on Olfert and Rogak, 2019
        K = self.PreFactor
        Dm = self.Exponent
        dmStar = pow((mStar * (10 ** -18)) / K, 1 / Dm) * (10 ** 9)
        self.start = pow(dmStar, self.f_lower)
        self.end = pow(dmStar, self.f_upper)
        return self.start and self.end

    def connect(self):
        if (not super().connect()):
            return False
        return (self.match(self.getRsp(), "Cambustion AAC"))

    def set(self, x):
        return self.sendFloat("SetSize", x)

    def setRorQsh(self, Qsh):
        self.Qsh = Qsh
        return self.sendFloat("SetSheath", Qsh)

    def getLine(self):
        start = time.time()
        while (self.scanchunks.find("\r\n") < 0 and time.time()-start < 100):
            try:
                if (self.isSerial):
                    chunk = self.ser.read().decode()
                else:
                    chunk = self.sock.recv(1024).decode()
            except:
                self.connected = False
                return ""
            self.scanchunks += chunk
        pos = self.scanchunks.find("\r\n")
        line = self.scanchunks[0:pos]
        self.scanchunks = self.scanchunks[pos+2:]
        return line

    def StartScan(self):

        self.sendCmd("SassScan s %d" % self.ScanUpTime + " %0.4E" % self.start+" %0.4E" %
                     self.end + " %0.4E" % self.RorQsh+" %0.4E" % self.delayTime+" %0.4E" % self.resAve+" u 1")
        self.scanchunks = ""
        while True:
            dummy = self.getLine()
            if (not self.match(dummy, "Cambustion")):
                return;
            if (self.match(dummy, "SCAN")):
                dummy = self.getLine()
                dummy = self.getLine()
                break
        self.mtc = True
        self.isScanning = True

    def StopScan(self):
        self.isScanning = False
        if (self.isSerial):
            self.ser.write(bytes([3, 13, 10]))  # ctrl-C
        else:
            self.sock.send(bytes([3, 13, 10]))  # ctrl-C

    def run(self):
        if (self.isScanner):
            return True
        else:
            return super().run()

    def stop(self):
        return super().stop()

    def getHeader(self):
        header = {}
        header['Serial number'] = self.sendQuery("serial")
        header['Sample flow (lpm)'] = self.sampleFlow
        header['Sheath flow SP (lpm)'] = self.Qsh
        header['Start (nm)'] = self.start
        header['End (nm)'] = self.end
        if (self.isScanner):
            header['Scan Time (s)'] = self.ScanUpTime
            header['Averaging'] = self.resAve
            header['Scan Delay Time (s)'] = self.delayTime
        else:
            header['Per decade'] = self.perdec
        self.headerFields = list(header.keys())
        return header

    def next(self):
        if (not self.isScanner):
            return super().next()
        self.values = {self.valueFields[i]: 0.0 for i in range(
            len(self.valueFields))}
        if (not self.isScanning):
            self.tally = 1
            return True
        line = self.getLine()
        aacfileline = line.split('\t')
        if (self.match(aacfileline[0], "END OF SCAN")):
            while True:
                dummy = self.getLine()
                if (self.match(dummy, "OK")):
                    break
            self.isScanning = False
            self.mtc = False
            return False
        self.x = float(aacfileline[2])
        self.point += 1
        if (self.point >= self.points):
            self.points += 1
            self.X.resize(self.points)
        self.X[self.point] = self.x
        self.values[self.quantity] = self.x
        self.conc = float(aacfileline[14])
        self.values["Speed (rad/s)"] = float(aacfileline[16])
        self.values["Sheath flow (lpm)"] = float(aacfileline[17])
        self.values["Pressure (Pa)"] = float(aacfileline[19])
        self.values["Temperature (C)"] = float(aacfileline[18])
        return True

    def moreToCome(self):
        if (self.isScanner):
            return self.mtc
        else:
            return super().moreToCome()

    def isReady(self):
        if (not self.isScanner):
            return super().isReady()
        return True

    def monitor(self):
        if (self.isScanner and self.isScanning):
            return True
        else:
            return super().monitor()


class TSIDma(Classifier):
    querydelay = 0.1
    quantity = "dm* (nm)"
    label = "dm*"

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh):
        super().__init__(isSerial, ip, serPort, start, end, perdec, flow, RorQsh)

    def sendAndCheck(self, cmd):
        self.lastResponse = self.sendQuery(cmd)
        return (not self.match(self.lastResponse, "ERROR"))

    def sendFloat(self, cmd, val):
        return (self.sendAndCheck(cmd + "%0.1F" % val))

    def setRorQsh(self, Qsh):
        self.Qsh = Qsh
        return self.sendFloat(self.sheathCommand, Qsh)

    def set(self, x):
        return self.sendFloat(self.sizeCommand, x)


class TSI3080(TSIDma):
    baudrate = 9600
    databits = 7
    stopbits = 1
    parity = "E"
    sizeCommand = "SPD"
    sheathCommand = "SQS"
    monitorCmd = "RMV"
    monitorFields = ["Dp", "Voltage (V)", "Sheath flow (lpm)", "Bypass flow", "Pressure (mbar)", "Temperature (C)", "Case temp", "Impactor flow", "E-mobility",
                     "Control mode", "Flow mode", "Sheath status", "Bypass status", "HV status", "Impactor Dp", "DMA model", "Gas", "Impactor", "Dp act", "Min Dp", "Max Dp"]
    valueFields = ["Sheath flow (lpm)", "Voltage (V)",
                   "Temperature (C)", "Pressure (mbar)"]
    monitorSep = ','

    def isReady(self):
        r = self.sendQuery("RFL")
        return (self.match(r, "1,1,1") or self.match(r, "1,0,1"))

    def getHeader(self):
        header = {}
        header['Sample flow (lpm)'] = self.sampleFlow
        header['Sheath flow SP (lpm)'] = self.Qsh
        header['Start (nm)'] = self.start
        header['End (nm)'] = self.end
        header['Per decade'] = self.perdec
        self.headerFields = list(header.keys())
        return header


class TSI3082(TSIDma):
    ipPort = 3602
    querydelay = 0.1

    def __init__(self, isSerial, ip, serPort, start, end, perdec, flow, RorQsh, HighFlow, Polarity, ScanUpTime, LowerRange, UpperRange, PreFactor, Exponent, f_lower, f_upper, IsWater, IsSoot, IsNone, variableBins, isScanner):
        super().__init__(isSerial, ip, serPort, start, end, perdec, flow, RorQsh)
        self.isScanner = isScanner
        self.sampleFlow = flow
        self.HighFlow = HighFlow
        self.Polarity = Polarity
        self.setPolarity(self.Polarity)
        self.ScanUpTime = ScanUpTime
        self.LowerRange = LowerRange
        self.UpperRange = UpperRange
        self.variableBins = variableBins
        self.IsWater = IsWater
        self.IsSoot = IsSoot
        self.IsNone = IsNone
        if f_lower and f_lower.strip():  # Check if not None and not an empty string
            self.f_lower = float(f_lower)
        else:
            self.f_lower = 0.75
        if f_upper and f_upper.strip():  # Check if not None and not an empty string
            self.f_upper = float(f_upper)
        else:
            self.f_upper = 1.25

        if self.IsWater:
            self.PreFactor = 523.6
            self.Exponent = 3
        elif self.IsSoot:
            self.PreFactor = 0.0612
            self.Exponent = 2.48
        elif self.IsNone:
            if PreFactor and PreFactor.strip():  # Check if not None and not an empty string
                self.PreFactor = float(PreFactor)
            else:
                self.PreFactor = 523.6

            if Exponent and Exponent.strip():  # Check if not None and not an empty string
                self.Exponent = float(Exponent)
            else:
                self.Exponent = 3

        if (self.isScanner):
            self.setHighFlow(self.HighFlow)
            self.setScanUpTime(self.ScanUpTime)
            self.setSampleFlow(self.sampleFlow)
            self.PurgeTime()
            self.SmpsUnits()
            self.SmpsWeighting()
            self.ReadBinBoundary()
            if (not self.variableBins):
                self.setLowerRange(self.LowerRange)
                self.setUpperRange(self.UpperRange)

    sizeCommand = "WSPARTICLEDIAM "  # note space at end
    scanUpTimeCommand = "WSSCANUPTIME "
    sampleflowCommand = "WSAEROSOLFLOW "
    sheathCommand = "WSSHFLOW "
    LowerRangeCommand = "WSLOWERSIZE "
    UpperRangeCommand = "WSUPPERSIZE "

    def __del__(self):
        if (self.isScanner):
            self.StopScan()
        super().__del__()

    def PurgeTime(self):
        return self.sendAndCheck("WSPURGETIME 0")

    def SmpsUnits(self):
        return self.sendAndCheck("WSSMPSUNITS 5")

    def SmpsWeighting(self):
        return self.sendAndCheck("WSSMPSWEIGHTS 0")

    def ReadBinBoundary(self):
        return self.sendAndCheck("RDSMPSDATA3")

    def setSampleFlow(self, flow):
        self.flow = flow
        return self.sendFloat("WSAEROSOLFLOW ", flow)

    def setHighFlow(self, HighFlow):
        self.HighFlow = HighFlow
        if (self.HighFlow):
            return self.sendAndCheck("WSDETINFLOW 1")
        else:
            return self.sendAndCheck("WSDETINFLOW 0")

    def setPolarity(self, Polarity):
        self.Polarity = Polarity
        if (self.Polarity):
            return self.sendAndCheck("WSHVPOL 0")
        else:
            return self.sendAndCheck("WSHVPOL 1")

    def setScanUpTime(self, ScanUpTime):
        return self.sendFloat("WSSCANUPTIME ", ScanUpTime)

    def setLowerRange(self, LowerRange):
        return self.sendAndCheck(f"WSLOWERSIZE {int(LowerRange)}")

    def setUpperRange(self, UpperRange):
        return self.sendAndCheck(f"WSUPPERSIZE {int(UpperRange)}")

    def VarDiameter(self, mStar):  # m = K * dm ^ Dm based on Olfert and Rogak, 2019
        K = self.PreFactor
        Dm = self.Exponent
        dmStar = pow((mStar * 1e-18) / K, 1 / Dm) * 1e9

        Lower = dmStar * self.f_lower
        Upper = dmStar * self.f_upper
        
        bin_edges = np.logspace(np.log10(1), np.log10(1000), num=193)
        bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
        TSIsetPoints = np.round(bin_centers, 2).tolist()

        LowerIndices = [i for i, value in enumerate(TSIsetPoints) if value <= Lower]
        LowerIndex = max(LowerIndices) if LowerIndices else 0

        UpperIndices = [i for i, value in enumerate(TSIsetPoints) if value < Upper]
        
        # Clamp maximum index to 191 (192 channels total)
        UpperIndex = min(191, max(UpperIndices) + 1) if UpperIndices else 191

        return self.setLowerRange(LowerIndex) and self.setUpperRange(UpperIndex)

    def LowerSizeRange(self):
        return self.sendQuery("RSLOWERSIZE")

    def UpperSizeRange(self):
        return self.sendQuery("RSUPPERSIZE")

    def ErrorStatus(self):
        return self.sendQuery("RMERRORS")

    def output3082(self):
        self.sendCmd("RDSMPSDATA4")
        time.sleep(0.3)  # Give the 3082 time to format the payload
        
        full_rsp = ""
        start_time = time.time()
        
        # Keep reading until 193+ lines arrive or timeout occurs (3 seconds)
        while time.time() - start_time < 3.0:
            if self.isSerial:
                chunk = self.ser.readline().decode(errors='ignore')
            else:
                try:
                    chunk = self.sock.recv(8192).decode(errors='ignore')
                except:
                    break
            
            full_rsp += chunk
            lines = full_rsp.strip().split('\r\n')
            if len(lines) >= 193:
                break
                
        return full_rsp.strip().split('\r\n')

    def StartScan(self):
        self.StopScan()  # Issues DOABORTSCAN to unlock settings buffer
        time.sleep(0.2)  # Grace period for firmware state change
        self.sendQuery("RMERRORS")
        self.Dataready3082 = 0
        return self.sendAndCheck("DOSCAN")

    def StopScan(self):
        return self.sendAndCheck("DOABORTSCAN")

    def DataReady(self):  # Gets the DMA status and check if it is at the end of scanning
        try:
            rsp = self.sendQuery("RDSMPSDATA1")
            self.ScanStatus = rsp.split('\r\n')
            
            if len(self.ScanStatus) > 8 and len(self.ScanStatus[8]) > 11:
                self.Dataready3082 = int(self.ScanStatus[8][11])
                return self.Dataready3082 == 1
        except (IndexError, ValueError, AttributeError):
            self.Dataready3082 = 0
            
        return False

    tolerance = 0.05
    # 3082 doesn't have a "return everything" command, so need to send a series of commands, which are in this dictionary of field:command
    monitorCmds = {"Sheath flow SP (lpm)": "RSSHFLOW", "Sheath flow (lpm)": "RMSHFLOW", "Voltage SP (V)": "RSHV",
                   "Voltage (V)": "RMHV", "Temperature (C)": "RMSHTEMP", "Pressure (kPa)": "RMSHAP", "Error Status": "RMERRORS"}
    # the fields we actually want in our file
    valueFields = ["Temperature (C)", "Pressure (kPa)"]

    def isReady(self):  # 3082 does't have nice status check command like 3080 does so we check if sheath and voltage are in tolerance
        if (self.isScanner):
            return True

        r = self.doMonitorCmd()

        ss = float(r['Sheath flow SP (lpm)'])
        sm = float(r['Sheath flow (lpm)'])
        vs = float(r['Voltage SP (V)'])
        vm = float(r['Voltage (V)'])

        return (abs((ss-sm)/ss) < self.tolerance and abs((vs-vm)/vs) < self.tolerance)

    def doMonitorCmd(self):
        output = self.monitorCmds.copy()
        for key, cmd in self.monitorCmds.items():
            output[key] = self.sendQuery(cmd)
        return output

    def getHeader(self):  # 3
        header = {}
        if (self.isScanner):
            header['Data points'] = int(
                self.UpperSizeRange()) - int(self.LowerSizeRange())
            header['Data length'] = len(self.fileFields) + 1
            header['Sample flow (lpm)'] = self.sampleFlow
            header['Sheath flow SP (lpm)'] = self.Qsh
            header['Error status'] = self.ErrorStatus()
        else:
            header['Sample flow (lpm)'] = self.sampleFlow
            header['Sheath flow SP (lpm)'] = self.Qsh
            header['Start (nm)'] = self.start
            header['End (nm)'] = self.end
            header['Per decade'] = self.perdec
        self.headerFields = list(header.keys())
        return header

    def reset(self):
        super().reset()
        if (self.isScanner):
            del self.ScanStatus


class CPC(Instrument):
    def conc(self):
        pass

    def startpoll(self):
        pass

    def endpoll(self):
        pass


class DummyCPC(CPC):
    def conc(self):
        return random.random()*1000.0

    def connect(self):
        return True

    def disconnect(self):
        return True

class CambustionCPC(CPC):
    ipPort = 23
    querydelay = 0.1
    timeout = 3.0
    baudrate = 115200
    databits = 8
    stopbits = 1
    parity = "N"
    
    def connect(self):
        if (not super().connect()):
            return False
        if (self.isSerial):
            self.ser.write(bytes([4]))
        return (self.match(self.getRsp(), "Cambustion CPC"))

    def getRsp(self):
        chunks = ""
        start = time.time()
        while (chunks.find("\r\n>") < 0 and time.time()-start < self.timeout):
            try:
                if (self.isSerial):
                    chunk = self.ser.read().decode()
                else:
                    chunk = self.sock.recv(1024).decode()
            except:
                self.connected = False
                return ""
            chunks += chunk
        return chunks[0:chunks.find("\r\n>")].replace(">", "")

    def sendFloat(self, cmd, val):
        return (self.sendAndCheck(cmd + " " + "%0.4E" % val))

    def sendAndCheck(self, cmd):
        self.lastResponse = self.sendQuery(cmd)
        return (self.match(self.lastResponse, "OK"))

    def disconnectEth(self):
        try:
            self.sock.send(bytes([13, 10, 4]))  # ctrl-D
            time.sleep(0.05)
            super().disconnectEth()
        except:
            pass
        
    def conc(self):
        return float(self.sendQuery("GCS"))
    
class MagicCpc(CPC):
    querydelay = 0.1
    baudrate = 115200
    databits = 8
    stopbits = 1
    parity = "N"

    def conc(self):
        return float(self.getRsp().split(',')[1])

    def startpoll(self):
        self.SendCmd("log,1")

    def endpoll(self):
        self.SendCmd("log,0")


class TSICpc(CPC):
    querydelay = 0.1

    def conc(self):
        if (self.isSerial):
            # base CPC class can only cope with serial.
            return float(self.sendQuery("RD"))


class TSI30xx(TSICpc):
    baudrate = 9600
    databits = 7
    stopbits = 1
    parity = "E"


class TSI377x(TSICpc):
    baudrate = 115200
    databits = 8
    stopbits = 1
    parity = "N"


class TSI375x(TSICpc):
    ipPort = 3603
    baudrate = 115200
    databits = 8
    stopbits = 1
    parity = "N"


class Arduino(Instrument):
    baudrate = 115200
    databits = 8
    stopbits = 1
    parity = "N"

    def enableBypass(self):
        if (self.connected):
            self.sendCmdSerial("on")
        return

    def disableBypass(self):
        if (self.connected):
            self.sendCmdSerial("off")
        return
