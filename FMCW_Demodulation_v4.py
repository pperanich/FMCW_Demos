from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget, QApplication, QSizePolicy

import pyqtgraph as pg
import numpy as np
import sys

FMCW_Demodulation_UI, FMCW_Demodulation_Base = uic.loadUiType("FMCW_Demodulation_v4.ui")

class FMCW_Demodulation_v4(FMCW_Demodulation_Base, FMCW_Demodulation_UI):
    __c = 3e8
    __fc = 79.0e9
    __lambdac = __c / __fc

    def __init__(self, parent: QWidget = None):
        # uic boilerplate
        FMCW_Demodulation_Base.__init__(self, parent=parent)
        self.setupUi(self)

        # UI items event connection
        self.update_button.clicked.connect(self.__update_graph)
        self.obj_1_range.editingFinished.connect(self.__obj_1_range_edited)
        self.obj_2_range.editingFinished.connect(self.__obj_2_range_edited)
        self.obj_3_range.editingFinished.connect(self.__obj_3_range_edited)
        self.chirp_bw.editingFinished.connect(self.__chirp_bw_edited)
        self.chirp_st.editingFinished.connect(self.__chirp_st_edited)
 
        # Setup plotting
        self.plot_window.setBackground(background="w")  # White background
        self.plot_window.setMenuEnabled(False)
        self.plot_window.setMouseEnabled(x=False, y=False)
        self.plot_window.showGrid(x=True, y=True)
        self.plot_window_legend = None

        # Radar control (not in GUI)
        self.gui_radar_max_range = 250.0               # Max range of radar in meters
        
        # Ranges of variables
        self.range_min = 0
        self.range_max = 200

        self.radar_chirp_bw_min = 0.2
        self.radar_chirp_bw_max = 4.0

        self.radar_chirp_sweep_time_min = 1.0
        self.radar_chirp_sweep_time_max = 100.0

    # Set range control functions
    @pyqtSlot(name="__obj_1_range_edited")
    def __obj_1_range_edited(self):
            try:
                temp = float(self.obj_1_range.text())
                if (temp >= self.range_min) & (temp <= self.range_max):
                    print("Set object 1 range to", self.obj_1_range.text(), " meters.")
                else:
                    print("Object 1 range must be in [",
                          self.range_min,
                          ",",
                          self.range_max,
                          "] meters.")
                    self.obj_1_range.setText(str((self.range_min+self.range_max)/2))
            except ValueError:
                print("Object 1 range must be a number.")
                self.obj_1_range.setText(str((self.range_min+self.range_max)/2))

    @pyqtSlot(name="__obj_2_range_edited")
    def __obj_2_range_edited(self):
        try:
            temp = float(self.obj_2_range.text())
            if (temp >= self.range_min) & (temp <= self.range_max):
                print("Set object 2 range to", self.obj_2_range.text(), " meters.")
            else:
                print("Object 2 range must be in [",
                      self.range_min,
                      ",",
                      self.range_max,
                      "] meters.")
                self.obj_2_range.setText(str((self.range_min+self.range_max)/2))
        except ValueError:
            print("Object 2 range must be a number.")
            self.obj_2_range.setText(str((self.range_min+self.range_max)/2))
            
    @pyqtSlot(name="__obj_3_range_edited")
    def __obj_3_range_edited(self):
            try:
                temp = float(self.obj_3_range.text())
                if (temp >= self.range_min) & (temp <= self.range_max):
                    print("Set object 3 range to", self.obj_3_range.text(), " meters.")
                else:
                    print("Object 3 range must be in [",
                          self.range_min,
                          ",",
                          self.range_max,
                          "] meters.")
                    self.obj_3_range.setText(str((self.range_min+self.range_max)/2))
            except ValueError:
                print("Object 3 range must be a number.")
                self.obj_3_range.setText(str((self.range_min+self.range_max)/2))

    # Radar control
    @pyqtSlot(name="__chirp_bw_edited")
    def __chirp_bw_edited(self):
        try:
            temp = float(self.chirp_bw.text())
            if (temp >= self.radar_chirp_bw_min) & (temp <= self.radar_chirp_bw_max):
                print("Set radar chirp bandwidth to", self.chirp_bw.text(), "GHz.")
            else:
                print("Radar chirp BW must be in [",
                      self.radar_chirp_bw_min,
                      ",",
                      self.radar_chirp_bw_max,
                      "] GHz.")
                self.chirp_bw.setText(str((self.radar_chirp_bw_min+self.radar_chirp_bw_max)/2))
        except ValueError:
            print("Radar chirp BW must be a number.")
            self.chirp_bw.setText(str((self.radar_chirp_bw_min+self.radar_chirp_bw_max)/2))

    @pyqtSlot(name="__chirp_st_edited")
    def __chirp_st_edited(self):
        try:
            temp = float(self.chirp_st.text())
            if (temp >= self.radar_chirp_sweep_time_min) & (temp <= self.radar_chirp_sweep_time_max):
                print("Set radar chirp sweep time to", self.chirp_st.text(), " microseconds.")
            else:
                print("Radar chirp sweep time must be in [",
                      self.radar_chirp_sweep_time_min,
                      ",",
                      self.radar_chirp_sweep_time_max,
                      "] microsec.")
                self.chirp_st.setText(str((self.radar_chirp_sweep_time_min+self.radar_chirp_sweep_time_max)/2))
        except ValueError:
            print("Radar chirp sweep time must be a number.")
            self.chirp_st.setText(str((self.radar_chirp_sweep_time_min+self.radar_chirp_sweep_time_max)/2))

    # The uber-function to update the graph
    @pyqtSlot(name="__update_graph")
    def __update_graph(self):
    
        # Enable
        obj_1_en = self.obj_1_en.isChecked()
        obj_2_en = self.obj_2_en.isChecked()
        obj_3_en = self.obj_3_en.isChecked()

        # Range
        obj_1_range = float(self.obj_1_range.text())
        obj_2_range = float(self.obj_2_range.text())
        obj_3_range = float(self.obj_3_range.text())

        # Radar control
        chirp_bw = float(self.chirp_bw.text())
        chirp_st = int(float(self.chirp_st.text()))

        # Plot control
        graph_zoom = self.graph_zoom.isChecked()

        # Constants for radar

        flo = FMCW_Demodulation_v4.__fc - chirp_bw * 1.0e9 / 2.0
        fhi = FMCW_Demodulation_v4.__fc + chirp_bw * 1.0e9 / 2.0
        
        c = FMCW_Demodulation_v4.__c

        target_range = np.array([obj_1_range, obj_2_range, obj_3_range])
        target_enable = np.array([obj_1_en, obj_2_en, obj_3_en])

        Tsw = chirp_st * 1.0e-6      # Wavelength in meters
        B = chirp_bw * 1.0e9                 # Chirp BW in Hz

        # Plot limits and number of points
        tmax = 100.0e-6
        npoints = 1000
        delta_t = tmax / npoints

        # TX Chirp: time limits are (t >= 0) and (t <= Tsw)
        tx_t = np.arange(0.0, Tsw+delta_t, delta_t)                   # Time in seconds
        tx_f = (flo + B / Tsw * tx_t) / 1.0e9                         # Frequency in GHz
        tx_t = tx_t * 1e6                                             # Time now in microseconds

        # RX signal, object #1: time limits are (t >= t1) & (t <= (t1 + Tsw)
        # Only if object 1 is enabled.
        if target_enable[0]:
            t1 = 2 * target_range[0] / c                              # Start time for rx signal 1 in microseconds
            rx_t1 = np.arange(t1, t1+Tsw+delta_t, delta_t)            # Time in seconds
            rx_f1 = (flo + B / Tsw * (rx_t1 - t1)) / 1.0e9            # Frequency in GHz
            rx_t1 = rx_t1 * 1e6                                       # Time now in microseconds
        else:
            rx_t1 = []
            rx_f1 = []

        # RX signal, object # 2: time limits are (t >= t2) & (t <= (t2 + Tsw)
        # Only if object 1 is enabled.
        if target_enable[1]:
            t2 = 2 * target_range[1] / c                              # Start time for rx signal 1 in microseconds
            rx_t2 = np.arange(t2, t2+Tsw+delta_t, delta_t)            # Time in seconds
            rx_f2 = (flo + B / Tsw * (rx_t2 - t2)) / 1.0e9            # Frequency in GHz
            rx_t2 = rx_t2 * 1e6                                       # Time now in microseconds
        else:
            rx_t2 = []
            rx_f2 = []

        # RX signal, object # 3: time limits are (t >= t2) & (t <= (t2 + Tsw)
        # Only if object 1 is enabled.
        if target_enable[2]:
            t3 = 2 * target_range[2] / c                              # Start time for rx signal 1 in microseconds
            rx_t3 = np.arange(t3, t3+Tsw+delta_t, delta_t)            # Time in seconds
            rx_f3 = (flo + B / Tsw * (rx_t3 - t3)) / 1.0e9            # Frequency in GHz
            rx_t3 = rx_t3 * 1e6                                       # Time now in microseconds
        else:
            rx_t3 = []
            rx_f3 = []

        # Set zoom
        if graph_zoom:
            zoom = [0.0, Tsw/10.0*1.0e6, flo/1.0e9, flo/1.0e9 + 0.5]
        else:
            zoom = [0.0, tmax*1.0e6, flo/1.0e9, flo/1.0e9 + 4.0]

        # Update the plot
        self.plot_window.clear()
        if self.plot_window_legend is not None:
            self.plot_window_legend.scene().removeItem(self.plot_window_legend)

        self.plot_window_legend = self.plot_window.addLegend(offset=(-5, 5))
        self.plot_window.setLabel('bottom', "t (microseconds)")
        self.plot_window.setLabel('left', "f (GHz)")
        self.plot_window.setTitle("Frequency vs. time Plot")
        self.plot_window.setXRange(zoom[0], zoom[1])
        self.plot_window.setYRange(zoom[2], zoom[3])
        
        tx_pen = pg.mkPen(color='b', width=3)
        self.plot_window.plot(tx_t, tx_f, name='TX', pen=tx_pen)
        if len(rx_t1) > 0:
            rx_t1_pen = pg.mkPen(color='r', width=3)
            self.plot_window.plot(rx_t1, rx_f1, name='RX1', pen=rx_t1_pen)
        if len(rx_t2) > 0:
            rx_t2_pen = pg.mkPen(color='g', width=3)
            self.plot_window.plot(rx_t2, rx_f2, name='RX2', pen=rx_t2_pen)
        if len(rx_t3) > 0:
            rx_t3_pen = pg.mkPen(color='y', width=3)
            self.plot_window.plot(rx_t3, rx_f3, name='RX3', pen=rx_t3_pen)
        print("Plot updated.\n")


if __name__ == "__main__":
    # Initialize application
    app = QApplication(sys.argv)
    mainWindow = FMCW_Demodulation_v4()
    mainWindow.show()
    sys.exit(app.exec_())
