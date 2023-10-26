# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 10:05:12 2023

Data visualization for QNC 2023

@author: Kristen Park 

"""
import os
import sleepy
import matplotlib.pyplot as plt
import numpy as np
import sleepy_mini #for freq analysis 

ppath = r'C:\Users\ChungWeberPC_23\Desktop\Kristen\KP01_warmingexp'
recFile = 'KP01_filelistingCE.txt'

(C, E) = sleepy.load_recordings(ppath, recFile)
whichfile = E[0]

#plot recordings (delta/theta/sigma/emg)
sleepy.calculate_spectrum(ppath, whichfile)
sleepy.sleep_state(ppath, whichfile, pplot=True)

#plot slow wave activity 
sleepy.plot_swa(ppath, whichfile, 60, 0.5, swa_yrange=[0, 0.012])

#sleep state dependent power spectral density analysis
sleepy_mini.basic_dreadd_analysis(ppath, recordingFile = recFile) #will compare two groups 

