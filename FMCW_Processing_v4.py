from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget, QApplication, QSizePolicy

import matplotlib.pyplot as plt
import pyqtgraph as pg
import numpy as np
import sys

FMCW_Processing_UI, FMCW_Processing_Base = uic.loadUiType("FMCW_Processing_v4.ui")

class FMCW_Processing_v4(FMCW_Processing_Base, FMCW_Processing_UI):

    def __init__(self, parent: QWidget = None):
        # uic boilerplate
        FMCW_Processing_Base.__init__(self, parent=parent)
        self.setupUi(self)

        # UI items event connection
        self.UpdatePlot.clicked.connect(self.graph_control_update_plot)
        
        self.Object1_Enable.toggled.connect(self.object1_enable_toggled)
        self.Object1_Range.editingFinished.connect(self.object1_range_edited)
        self.Object1_Radial_Velocity.editingFinished.connect(self.object1_radial_velocity_edited)
        self.Object1_RCS.editingFinished.connect(self.object1_rcs_edited)
        
        self.Object2_Enable.toggled.connect(self.object2_enable_toggled)
        self.Object2_Range.editingFinished.connect(self.object2_range_edited)
        self.Object2_Radial_Velocity.editingFinished.connect(self.object2_radial_velocity_edited)
        self.Object2_RCS.editingFinished.connect(self.object2_rcs_edited)
        
        self.Object3_Enable.toggled.connect(self.object3_enable_toggled)
        self.Object3_Range.editingFinished.connect(self.object3_range_edited)
        self.Object3_Radial_Velocity.editingFinished.connect(self.object3_radial_velocity_edited)
        self.Object3_RCS.editingFinished.connect(self.object3_rcs_edited)
        
        self.ThermalNoise_Enable.toggled.connect(self.thermal_noise_enable_toggled)
        
        self.Radar_Max_Velocity.editingFinished.connect(self.radar_max_velocity_edited)
        self.Radar_Chirp_BW.editingFinished.connect(self.radar_chirp_BW_edited)
        self.Radar_Chirps_per_Frame.editingFinished.connect(self.radar_chirps_per_frame_edited)
        self.Radar_Range_FFT_Size.editingFinished.connect(self.radar_range_fft_size_edited)
        self.Radar_Doppler_FFT_Size.editingFinished.connect(self.radar_doppler_fft_size_edited)
  
        # Setup plotting
        self.Range_IF_Plot.setBackground(background="w")  # White background
        self.Range_IF_Plot.setMenuEnabled(False)
        self.Range_IF_Plot.setMouseEnabled(x=False, y=False)
        self.Range_IF_Plot.showGrid(x=True, y=True)
        
        self.Range_Doppler_Plot.setBackground(background="w")  # White background
        self.Range_Doppler_Plot.setMenuEnabled(False)
        self.Range_Doppler_Plot.setMouseEnabled(x=False, y=False)
        
        """
            Declare variables to control demo
        """
        self.object1_enabled = self.Object1_Enable.isChecked()
        self.object2_enabled = self.Object2_Enable.isChecked()
        self.object3_enabled = self.Object3_Enable.isChecked()
        self.object_thermal_noise_enabled = self.ThermalNoise_Enable.isChecked()

        self.object1_range = float(self.Object1_Range.text())
        self.object2_range = float(self.Object2_Range.text())
        self.object3_range = float(self.Object3_Range.text())

        self.object1_radial_velocity = float(self.Object1_Radial_Velocity.text())
        self.object2_radial_velocity = float(self.Object2_Radial_Velocity.text())
        self.object3_radial_velocity = float(self.Object3_Radial_Velocity.text())

        self.object1_rcs = float(self.Object1_RCS.text())
        self.object2_rcs = float(self.Object2_RCS.text())
        self.object3_rcs = float(self.Object3_RCS.text())

        # Radar control
        self.radar_max_velocity = float(self.Radar_Max_Velocity.text())
        self.radar_chirp_bw = float(self.Radar_Chirp_BW.text())
        self.radar_chirps_per_frame = int(self.Radar_Chirps_per_Frame.text())
        self.radar_range_fft_size = int(self.Radar_Range_FFT_Size.text())
        self.radar_Doppler_fft_size = int(self.Radar_Doppler_FFT_Size.text())

        # Radar control (not in GUI)
        self.radar_max_range = 250.0               # Max range of radar in meters
        self.radar_num_receive_channels = int(21)  # Number of RX channels

        # Constants for radar
        self.c = 3e8
        self.fc = 79.0e9
        self.flo = self.fc - self.radar_chirp_bw * 1.0e9 / 2.0
        self.fhi = self.fc + self.radar_chirp_bw * 1.0e9 / 2.0
        self.lambdac = self.c / self.fc

        # Ranges of variables
        self.range_min = 0
        self.range_max = 200

        self.radial_velocity_min = -100
        self.radial_velocity_max = 100

        self.rcs_min = -20
        self.rcs_max = 40

        self.radar_max_velocity_min = 20
        self.radar_max_velocity_max = 200

        self.radar_chirp_bw_min = 20
        self.radar_chirp_bw_max = 400

        self.radar_chirps_per_frame_min = 1
        self.radar_chirps_per_frame_max = 100

        self.radar_range_fft_size_min = 8
        self.radar_range_fft_size_max = 1024

        self.radar_Doppler_fft_size_min = 8
        self.radar_Doppler_fft_size_max = 1024

    #
    # Object control
    #
    @pyqtSlot(name="object1_enable_toggled")
    def object1_enable_toggled(self):
        self.object1_enabled = self.Object1_Enable.isChecked()
        if (self.object1_enabled):
            print("Object 1 is enabled.")
        else:
            print("Object 1 is disabled.")

    @pyqtSlot(name="object2_enable_toggled")
    def object2_enable_toggled(self):
        self.object2_enabled = self.Object2_Enable.isChecked()
        if (self.object2_enabled):
            print("Object 2 is enabled.")
        else:
            print("Object 2 is disabled.")

    @pyqtSlot(name="object3_enable_toggled")
    def object3_enable_toggled(self):
        self.object3_enabled = self.Object3_Enable.isChecked()
        if (self.object3_enabled):
            print("Object 3 is enabled.")
        else:
            print("Object 3 is disabled.")

    @pyqtSlot(name="thermal_noise_enable_toggled")
    def thermal_noise_enable_toggled(self):
        self.object_thermal_noise_enabled = self.ThermalNoise_Enable.isChecked()
        if (self.object_thermal_noise_enabled):
            print("Thermal noise is enabled.")
        else:
            print("Thermal noise is disabled.")

    #
    # Range control functions
    #
    @pyqtSlot(name="object1_range_edited")
    def object1_range_edited(self):
        try:
            temp = float(self.Object1_Range.text())
            if (temp >= self.range_min) & (temp <= self.range_max):
                self.object1_range = temp
                print("Set object 1 range to", self.Object1_Range.text(), "meters.")
            else:
                self.Object1_Range.setText(str(self.object1_range))
                print("Object 1 range must be in [",
                      self.range_min,
                      ",",
                      self.range_max,
                      "] meters.")
                self.Object1_Range.setText(str(self.object1_range))
        except ValueError:
            self.Object1_Range.setText(str(self.object1_range))

    @pyqtSlot(name="object2_range_edited")
    def object2_range_edited(self):
        try:
            temp = float(self.Object2_Range.text())
            if (temp >= self.range_min) & (temp <= self.range_max):
                self.object2_range = temp
                print("Set object 2 range to", self.Object2_Range.text(), "meters.")
            else:
                print("Object 2 range must be in [",
                      self.range_min,
                      ",",
                      self.range_max,
                      "] meters.")
                self.Object2_Range.setText(str(self.object2_range))
        except ValueError:
            print("Object 2 range must be a number.")
            self.Object2_Range.setText(str(self.object2_range))

    @pyqtSlot(name="object3_range_edited")
    def object3_range_edited(self):
        try:
            temp = float(self.Object3_Range.text())
            if (temp >= self.range_min) & (temp <= self.range_max):
                self.object3_range = temp
                print("Set object 3 range to", self.Object3_Range.text(), "meters.")
            else:
                print("Object 3 range must be in [",
                      self.range_min,
                      ",",
                      self.range_max,
                      "] meters.")
                self.Object3_Range.setText(str(self.object3_range))
        except ValueError:
            print("Object 3 range must be a number.")
            self.Object3_Range.setText(str(self.object3_range))


    #
    # Radial velocity control functions
    #
    @pyqtSlot(name="object1_radial_velocity_edited")
    def object1_radial_velocity_edited(self):
        try:
            temp = float(self.Object1_Radial_Velocity.text())
            if (temp >= self.radial_velocity_min) & (temp <= self.radial_velocity_max):
                self.object1_radial_velocity = temp
                print("Set object 1 radial velocity to", self.Object1_Radial_Velocity.text(), "meters/second.")
            else:
                print("Object 1 radial velocity must be in [",
                      self.radial_velocity_min,
                      ",",
                      self.radial_velocity_max,
                      "] meters/second.")
                self.Object1_Radial_Velocity.setText(str(self.object1_radial_velocity))
        except ValueError:
            print("Object 1 radial velocity must be a number.")
            self.Object1_Radial_Velocity.setText(str(self.object1_radial_velocity))

    @pyqtSlot(name="object2_radial_velocity_edited")
    def object2_radial_velocity_edited(self):
        try:
            temp = float(self.Object2_Radial_Velocity.text())
            if (temp >= self.radial_velocity_min) & (temp <= self.radial_velocity_max):
                self.object2_radial_velocity = temp
                print("Set object 2 radial velocity to", self.Object2_Radial_Velocity.text(), "meters/second.")
            else:
                print("Object 2 radial velocity must be in [",
                      self.radial_velocity_min,
                      ",",
                      self.radial_velocity_max,
                      "] meters/second.")
                self.Object2_Radial_Velocity.setText(str(self.object2_radial_velocity))
        except ValueError:
            print("Object 2 radial velocity must be a number.")
            self.Object2_Radial_Velocity.setText(str(self.object2_radial_velocity))

    @pyqtSlot(name="object3_radial_velocity_edited")
    def object3_radial_velocity_edited(self):
        try:
            temp = float(self.Object3_Radial_Velocity.text())
            if (temp >= self.radial_velocity_min) & (temp <= self.radial_velocity_max):
                self.object1_radial_velocity = temp
                print("Set object 3 radial velocity to", self.Object3_Radial_Velocity.text(), "meters/second.")
            else:
                print("Object 3 radial velocity must be in [",
                      self.radial_velocity_min,
                      ",",
                      self.radial_velocity_max,
                      "] meters/second.")
                self.Object3_Radial_Velocity.setText(str(self.object3_radial_velocity))
        except ValueError:
            print("Object 3 radial velocity must be a number.")
            self.Object3_Radial_Velocity.setText(str(self.object3_radial_velocity))


    #
    # RCS control functions
    #
    @pyqtSlot(name="object1_rcs_edited")
    def object1_rcs_edited(self):
        try:
            temp = float(self.Object1_RCS.text())
            if (temp >= self.rcs_min) & (temp <= self.rcs_max):
                self.object1_rcs = temp
                print("Set object 1 RCS to ", temp, "dBsm.")
            else:
                print("Object 1 RCS must be in [",
                      self.rcs_min,
                      ",",
                      self.rcs_max,
                      "] dBsm.")
                self.Object1_RCS.setText(str(self.object1_rcs))
        except ValueError:
            print("Object 1 RCS must be a number.")
            self.Object1_RCS.setText(str(self.object1_rcs))

    @pyqtSlot(name="object2_rcs_edited")
    def object2_rcs_edited(self):
        try:
            temp = float(self.Object2_RCS.text())
            if (temp >= self.rcs_min) & (temp <= self.rcs_max):
                self.object2_rcs = temp
                print("Set object 2 RCS to ", temp, "dBsm.")
            else:
                print("Object 2 RCS must be in [",
                      self.rcs_min,
                      ",",
                      self.rcs_max,
                      "] dBsm.")
                self.Object2_RCS.setText(str(self.object2_rcs))
        except ValueError:
            print("Object 2 RCS must be a number.")
            self.Object2_RCS.setText(str(self.object2_rcs))

    @pyqtSlot(name="object3_rcs_edited")
    def object3_rcs_edited(self):
        try:
            temp = float(self.Object3_RCS.text())
            if (temp >= self.rcs_min) & (temp <= self.rcs_max):
                self.object3_rcs = temp
                print("Set object 3 RCS to ", temp, "dBsm.")
            else:
                print("Object 3 RCS must be in [",
                      self.rcs_min,
                      ",",
                      self.rcs_max,
                      "] dBsm.")
                self.Object3_RCS.setText(str(self.object3_rcs))
        except ValueError:
            print("Object 3 RCS must be a number.")
            self.Object3_RCS.setText(str(self.object3_rcs))


    #
    # Radar control
    #
    @pyqtSlot(name="radar_max_velocity_edited")
    def radar_max_velocity_edited(self):
        try:
            temp = float(self.Radar_Max_Velocity.text())
            if (temp >= self.radar_max_velocity_min) & (temp <= self.radar_max_velocity_max):
                self.radar_max_velocity = temp
                print("Set radar max velocity to", self.Radar_Max_Velocity.text(), "meters/second.")
            else:
                print("Radar max velocity must be in [",
                      self.radar_max_velocity_min,
                      ",",
                      self.radar_max_velocity_max,
                      "] meters/second.")
                self.Radar_Max_Velocity.setText(str(self.radar_max_velocity))
        except ValueError:
            print("Radar max velocity must be a number.")
            self.Radar_Max_Velocity.setText(str(self.radar_max_velocity))

    @pyqtSlot(name="radar_chirp_BW_edited")
    def radar_chirp_BW_edited(self):
        try:
            temp = float(self.Radar_Chirp_BW.text())
            if (temp >= self.radar_chirp_bw_min) & (temp <= self.radar_chirp_bw_max):
                self.radar_chirp_bw = temp
                print("Set radar chirp bandwidth to", self.Radar_Chirp_BW.text(), "MHz.")
            else:
                print("Radar chirp BW must be in [",
                      self.radar_chirp_bw_min,
                      ",",
                      self.radar_chirp_bw_max,
                      "] MHz.")
                self.Radar_Chirp_BW.setText(str(self.radar_chirp_bw))
        except ValueError:
            print("Radar chirp BW must be a number.")
            self.Radar_Chirp_BW.setText(str(self.radar_chirp_bw))

    @pyqtSlot(name="radar_chirps_per_frame_edited")
    def radar_chirps_per_frame_edited(self):
        try:
            temp = int(self.Radar_Chirps_per_Frame.text())
            if (temp >= self.radar_chirps_per_frame_min) & (temp <= self.radar_chirps_per_frame_max):
                self.radar_chirps_per_frame = temp
                print("Set radar chirps per frame to", self.Radar_Chirps_per_Frame.text(), ".")
            else:
                print("Radar chirp BW must be in [",
                      self.radar_chirps_per_frame_min,
                      ",",
                      self.radar_chirps_per_frame_max,
                      "] MHz.")
                self.Radar_Chirps_per_Frame.setText(str(self.radar_chirps_per_frame))
        except ValueError:
            print("Radar chirps per frame must be a number.")
            self.Radar_Chirps_per_Frame.setText(str(self.radar_chirps_per_frame))

    @pyqtSlot(name="radar_range_fft_size_edited")
    def radar_range_fft_size_edited(self):
        try:
            temp = int(self.Radar_Range_FFT_Size.text())
            if (temp >= self.radar_range_fft_size_min) & (temp <= self.radar_range_fft_size_max):
                self.radar_range_fft_size = temp
                print("Set radar range FFT size to", self.Radar_Range_FFT_Size.text(), ".")
            else:
                print("Radar range FFT size must be in [",
                      self.radar_range_fft_size_min,
                      ",",
                      self.radar_range_fft_size_max,
                      "].")
                self.Radar_Range_FFT_Size.setText(str(self.radar_range_fft_size))
        except ValueError:
            print("Radar range FFT size must be a number.")
            self.Radar_Range_FFT_Size.setText(str(self.radar_range_fft_size))

    @pyqtSlot(name="radar_doppler_fft_size_edited")
    def radar_doppler_fft_size_edited(self):
        try:
            temp = int(self.Radar_Doppler_FFT_Size.text())
            if (temp >= self.radar_Doppler_fft_size_min) & (temp <= self.radar_Doppler_fft_size_max):
                self.radar_Doppler_fft_size = temp
                print("Set radar Doppler FFT size to", self.Radar_Doppler_FFT_Size.text(), ".")
            else:
                print("Radar Doppler FFT size must be in [",
                      self.radar_Doppler_fft_size_min,
                      ",",
                      self.radar_Doppler_fft_size_max,
                      "].")
                self.Radar_Doppler_FFT_Size.setText(str(self.radar_range_fft_size))
        except ValueError:
            print("Radar Doppler FFT size must be a number.")
            self.Radar_Doppler_FFT_Size.setText(str(self.radar_range_fft_size))

    #
    # The uber-function to update the graph
    #
    @pyqtSlot(name="graph_control_update_plot")
    def graph_control_update_plot(self):

        target_range = np.array([self.object1_range,
                                 self.object2_range,
                                 self.object3_range])

        target_enable = np.array([self.object1_enabled,
                                  self.object2_enabled,
                                  self.object3_enabled])


        lambdac = self.lambdac                             # Wavelength in meters
        B = self.radar_chirp_bw * 1e6                  # Chirp BW in Hz
        c = self.c                                         # Speed of light in m/s

        fc = self.fc                                       # Radar center frequency in Hz
        flo = self.flo                                     # Radar low frequency in Hz
        fhi = self.fhi                                     # Radar high frequency in Hz
        lambdalo = c / flo                                 # Wavelength with flo
        lambdahi = c / fhi                                 # Wavelength with fhi
        lambdamin = min(lambdalo, lambdahi)                # Pick lambda for radar

        L = self.radar_num_receive_channels            # Number of RX channels
        dx = 0.5 * lambdamin                               # Receiver element spacing
        rmax = self.radar_max_range                    # Max radar range in meters
        vmax = self.radar_max_velocity                 # Max velocity in meters/second
        Tsw = 1.0 / (2.0 * vmax / lambdac)                 # Sweep time in seconds
        P = self.radar_chirps_per_frame                # Chirps per frame
        PRI = 1.0 / Tsw                                    # PRI
        fif_max = 2.0 * B * rmax / (c * Tsw)               # max IF frequency
        fs = 2.0 * fif_max                                 # sampling rate in Hz
        Nrange = self.radar_range_fft_size             # Range FFT size
        NDoppler = self.radar_Doppler_fft_size         # Doppler FFT size

        # Fast time indices
        dr = c / (2.0 * B)  # Range resolution
        N = int(np.ceil(rmax / dr))  # Number of rast-time range samples

        # FMCW constant
        K = B / Tsw  # FMCW constant

        # Create signal (complex)
        ys = np.zeros(shape=(L, N, P)) + (1j * np.zeros(shape=(L, N, P)))

        # Set up loop variables
        kloop = [1, 2, 3]
        lloop = range(1, L + 1)  # 1:L
        nloop = range(1, N + 1)  # 1:N
        ploop = range(1, P + 1)  # 1:P

        for k in kloop:

            # k is index for object number
            if k == 1:
                r = self.object1_range
                v = self.object1_radial_velocity
                if self.object1_enabled:
                    rcs = 10 ** (self.object1_rcs / 10.0)
                else:
                    rcs = 0.0
            elif k == 2:
                r = self.object2_range
                v = self.object2_radial_velocity
                if self.object2_enabled:
                    rcs = 10 ** (self.object2_rcs / 10.0)
                else:
                    rcs = 0.0
            else:
                r = self.object3_range
                v = self.object3_radial_velocity
                if self.object3_enabled:
                    rcs = 10 ** (self.object3_rcs / 10.0)
                else:
                    rcs = 0.0

            print("\n Simulating object", k, "return.")
            print("   r   =", r, 'm')
            print("   rcs =", rcs, 'sm')
            print("   v   =", v, 'm/s')

            # Doppler shift (fd) of q-th target
            fdq = 2.0 * v / lambdac

            # Beamforming - assume all objects just in front of radar
            u = np.sin(0.0)

            for l in lloop:
                for n in nloop:
                    for p in ploop:
                        # Signal - phase due to object range/Doppler
                        theta = 2 * np.pi * (
                                    (2 * K * r / c + fdq) * ((n - 1) / fs) + fc * (l - 1) * dx * u / c + fdq * (
                                        p - 1) * Tsw + 2 * fc * r / c)
                        ys[l - 1, n - 1, p - 1] = ys[l - 1, n - 1, p - 1] + rcs * np.exp(1j * theta)

            # Noise
            if self.object_thermal_noise_enabled:
                for l in lloop:
                    for n in nloop:
                        for p in ploop:
                            rand_real = np.random.randn(1, 1)
                            rand_imag = np.random.randn(1, 1)
                            ys[l - 1, n - 1, p - 1] = ys[l - 1, n - 1, p - 1] + rand_real + 1j * rand_imag

            # Beamform: Sum across first dimension
            yb = np.zeros(shape=(N, P)) + (1j * np.zeros(shape=(N, P)))
            for l in lloop:
                for n in nloop:
                    for p in ploop:
                        yb[n - 1, p - 1] = yb[n - 1, p - 1] + ys[l - 1, n - 1, p - 1]
            yb = yb / L

            # Range FFT - take 1D FFT across range bins
            # All zero complex matrix size(Nrange X P)
            yr = np.zeros(shape=(Nrange, P)) + (1j * np.zeros(shape=(Nrange, P)))
            yr_shift = np.zeros(shape=(Nrange, P)) + (1j * np.zeros(shape=(Nrange, P)))

            for p in ploop:
                temp = yb[:, p - 1]
                temp_fft = np.fft.fft(temp, Nrange)
                yr_shift[:, p - 1] = np.fft.fftshift(temp_fft)
                yr[:, p - 1] = temp_fft

            # Perform FFT on range processed data to get Doppler
            yr_rows = int(yr.shape[0])
            yrd = np.zeros(shape=(yr_rows, NDoppler)) + (1j * np.zeros(shape=(yr_rows, NDoppler)))
            yrd_shift = np.zeros(shape=(yr_rows, NDoppler)) + (1j * np.zeros(shape=(yr_rows, NDoppler)))
            yrloop = range(1, yr_rows + 1)  # 1:rows in yr

            for n in yrloop:
                temp = yr[n - 1, :]
                temp_fft = np.fft.fft(temp, NDoppler)
                yrd_shift[n - 1, :] = np.fft.fftshift(temp_fft)
                yrd[n - 1, :] = temp_fft

            # Perform FFT on range processed data to get Doppler
            yrd = np.fft.fftshift(yrd, 1)
            yrd_shift = np.fft.fftshift(yrd_shift, 1)


        #
        # Update the Range IF Plot
        #
        yr_dB = 10.0 * np.log10(abs(yr) + 1e-16)
        rangeplot_vals = np.linspace(0, ((Nrange/2) - 1) / Nrange, len(yr_dB))  # 0 to (Nrange/2)-1
        fft_range_vals = rmax * 2.0 * rangeplot_vals  
        
        self.Range_IF_Plot.clear()
        self.Range_IF_Plot.setLabel('bottom', "Range (m)")
        self.Range_IF_Plot.setLabel('left', "Level (dB)")
        self.Range_IF_Plot.setTitle("Range IF FFT Plot")
        self.Range_IF_Plot.setYRange(-10, 60)
        
        if_pen = pg.mkPen(color='b', width=3)
        self.Range_IF_Plot.plot(fft_range_vals, yr_dB[0:(len(fft_range_vals)), 1], pen=if_pen)                          

        #
        # Update Range Doppler Plot
        #
        RR = yrd.shape[0]  # Rows in yrd_shift_dB
        yrd_indx = range(1, int(np.floor(RR/2)))
        yrd_plot = yrd[yrd_indx, :]
        yrd_plot_dB = 10.0 * np.log10(abs(yrd_plot) + 1e-6)

        # Calculate the velocity labels and tick marks for plotting
        # In multiples of 10 m/s
        vmax_half_10 = 10.0 * np.floor(vmax / 2.0 / 10.0)
        fft_velocity_labels = np.arange(-vmax_half_10, vmax_half_10+1, 10.0)
        fft_velocity_ticks = fft_velocity_labels * NDoppler / vmax + (NDoppler / 2.0)

        fft_velocity_ticks_fixed = ["%d" % int(float(l)) for l in fft_velocity_labels]
        fft_velocity_tuple = []
        for idx in range(len(fft_velocity_ticks)):
            fft_velocity_tuple.append((fft_velocity_ticks[idx],fft_velocity_ticks_fixed[idx]))
        
        self.Range_Doppler_Plot.clear()
        self.Range_Doppler_Plot.setLabel('bottom', "Velocity(m/s)")
        self.Range_Doppler_Plot.setLabel('left', "Range (m)")
        self.Range_Doppler_Plot.setTitle("Range Doppler 2D Plot")
        
        # Get the colormap
        colors = plt.get_cmap('viridis')
        colors._init()
        lut = (colors._lut * 255).view(np.ndarray)  # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
        
        # Apply the colormap
        img = pg.ImageItem(abs(yrd_plot).T)
        img.setLookupTable(lut)
        
        # Add image to plot window
        self.Range_Doppler_Plot.addItem(img)
        
        # Correct axis ticks
        pltItem = self.Range_Doppler_Plot.getPlotItem()
        axis = pltItem.getAxis('bottom')
        axis.setTicks([fft_velocity_tuple])

        print("\nGUI Plot Updated.\n")

if __name__ == "__main__":
    # Initialize application
    app = QApplication(sys.argv)
    mainWindow = FMCW_Processing_v4()
    mainWindow.show()
    sys.exit(app.exec_())
