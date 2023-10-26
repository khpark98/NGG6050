# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:37:49 2019

@author: snowi
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 08 09:00:49 2018

@author: Joe S.
"""

import scipy.signal
import numpy as np
import numpy.random as rand
import scipy.io as so
import scipy.stats as stats
import os.path
import re #regular expressions module
import matplotlib.pylab as plt
import h5py
import matplotlib.patches as patches
import pdb
import pandas as pd
import seaborn as sns
from functools import reduce
import sleepy
#import pyphi_new as pyphi
import pyphi
import pingouin as pg
import gc

'''
Run all the basic plotting functions on a sleep file that has been run through sleep_annotation_qt.py
'''
def basic_opto_analysis(ppath, recordings):
    sleepy.laser_brainstate(ppath, recordings, 300, 600, start_time=0, end_time=-1)  #200, 250 , #first half of 6 hour recording 0, 10799
    # sleepy.laser_brainstate_bootstrap(ppath, recordings, 300, 600, nboots=1000,bootstrap_mode=0, ma_thr=10, alpha=0.1)
#    if type(recordings) == list:
#        sleepy.laser_triggered_eeg_avg(ppath, recordings, pre=200, post=250, f_max=30, laser_dur=120) #laser dur is number of seconds laser train is active
#        sleepy.laser_brainstate_bootstrap(ppath, recordings, 200, 320, nboots=1000,bootstrap_mode=0, ma_thr=10, alpha=0.1)
#        sleepy.laser_brainstate_bootstrap(ppath, recordings, 600, 720, nboots=1000,bootstrap_mode=0, ma_thr=10, alpha=0.1)
#        #command to get good, simple version of spectogram
#        sleepy.laser_triggered_eeg_avg(ppath, recordings, pre=200, post=320, f_max=30, laser_dur=120, pplot = 2)
#        sleepy.transition_analysis(ppath, recordings, 120, 120, 10, 30, after_laser=200, bootstrap_mode=0, paired_stats=True)
#        sleepy.transition_markov_strength(ppath, recordings, pre = 120, laser_tend = 120, tdown = 10, nboot=1000, stats_lastpoint=False, pstationary=False, paired_stats=True, ma_thr=10)
#    else:
#        sleepy.laser_triggered_eeg(ppath, recordings, 200, 250, f_max = 30)   
    # sleepy.sleep_spectrum_simple(ppath, recordings, istate = 1)
    # sleepy.sleep_spectrum_simple(ppath, recordings, istate = 2)
    # sleepy.sleep_spectrum_simple(ppath, recordings, istate = 3)
    sleepy.sleep_spectrum(ppath, recordings, istate = 1)
    sleepy.sleep_spectrum(ppath, recordings, istate = 2)
    sleepy.sleep_spectrum(ppath, recordings, istate = 3)
#    sleepy.sleep_spectrum_simple(ppath, recordings, istate=2, tstart=0, tend=-1, fmax=-1, mu=[10,100], ci='sd', pmode=1, pnorm = False, pplot=True, harmcs=10, pemg2=False, csv_files=['eeg_file.csv', 'emg_file.csv'])
#    sleepy.laser_triggered_eeg_avg(ppath, recordings, pre=200, post=250, f_max=30, laser_dur=120, pplot = 2) #good version of laser triggered eeg spectogram
    
    

'''
Run all the stats for dreadd recordings all at once so it's consolidated and I don't have to keep typing different functions 
''' 
def basic_dreadd_analysis(ppath, recordingFile = '', tstart = 0, tend = -1):
    #timecourse info
    # sleepy.sleep_timecourse(ppath, trace_file = recordingFile, tbin = 1800, n = 10, tstart = tstart, tend = tend, pplot=True, stats='perc', csv_file='')
    # sleepy.sleep_timecourse(ppath, trace_file = recordingFile, tbin = 1800, n = 10, tstart = tstart, tend = tend, pplot=True, stats='freq', csv_file='')
    # sleepy.sleep_timecourse(ppath, trace_file = recordingFile, tbin = 1800, n = 10, tstart = tstart, tend = tend, pplot=True, stats='dur', csv_file='')
    
    #spectrum comparison
    compare_spectrum_conditions(ppath, recordingFile = recordingFile, tstart = tstart, tend = tend)
    
    
'''
Function for running comparing dreadd spectrums using custom sleep_spectrum function present in this file
Will only plot 2 conditions. For more detailed script that can plot multiple conditions, use spectrum_conditions.py / spectrum_conditions_transitions.py
'''
############## C O M P A R E   D R E A D D S   W I T H   S L E E P S P E C T R U M ##########################
def compare_spectrum_conditions(ppath, controlRecordings = [], experimentalRecordings = [], recordingFile = '', tstart = 0, tend = -1, f_max = 16):    
    if recordingFile != '':
        controlRecordings, experimentalRecordingsDict = load_dose_recordings(ppath, recordingFile) 
        #load dose recordings returns the control recordings as a list and the experimental recordings as a dicitonary with doses
        #so we have to pull out the recordings associated with their 'dose'
        for key,value in experimentalRecordingsDict.items():
            experimentalRecordings += experimentalRecordingsDict[key] #lumps all the doses together into one spectrum
 
    sleep_spectrum_dreadds(ppath, controlRecordings, experimentalRecordings, istate=1, f_max = f_max, tstart = tstart, tend = tend)
    sleep_spectrum_dreadds(ppath, controlRecordings, experimentalRecordings, istate=2, f_max = f_max, tstart = tstart, tend = tend)
    sleep_spectrum_dreadds(ppath, controlRecordings, experimentalRecordings, istate=3, f_max = f_max, tstart = tstart, tend = tend)
#    sleep_spectrum_dreadds(ppath, controlRecordings, experimentalRecordings, istate=4, f_max = 30, tstart = tstart, tend = tend) #for transitions
    
'''
Run the basic photometry analysis scripts, each file must be run through sleep_annotation_qt beforehand, as per usual.
''' 
def basic_photometry_analysis(ppath, recordings, rawTracesForList = False):
    if type(recordings) == list:
        for eachRecording in recordings:
            pyphi.calculate_dff(ppath, eachRecording, perc=30) 
            if rawTracesForList:
                pyphi.plot_rawtraces(ppath, eachRecording)
        pyphi.avg_activity_recs(ppath, recordings, pzscore = True)
        pyphi.avg_activity_recs(ppath, recordings, pzscore = False)
        
    else:
        pyphi.calculate_dff(ppath, recordings, perc=30)
        pyphi.plot_rawtraces(ppath, recordings, pzscore = False)
        pyphi.plot_rawtraces(ppath, recordings, pzscore = True)
        pyphi.avg_activity(ppath, recordings)
        _=pyphi.spectralfield_highres(ppath, recordings, 4, 4, fmax=20, theta=[50, 100, 1000, 1000, 10000], states=[1])
    pyphi.activity_transitions(ppath, recordings, [(3,1), (1, 2), (3,2), (2,3)], 30, 30, [20, 60, 60], [20, 0, 0], pzscore = False, ylim=(0, 30)) #z-score false
    pyphi.activity_transitions(ppath, recordings, [(3,1), (1, 2), (3,2), (2,3)], 30, 30, [20, 60, 60], [20, 0, 0], pzscore = True, ylim=(-1,2)) #z-score true
    pyphi.dff_sleepcycle(ppath, recordings)
    pyphi.dff_vs_statedur(ppath, recordings, istate = 1, pzscore=True, sf=0, fs=1.5, ma_thr=0, backup='')


def online_homeostasis_stats(ppath, files):
    df = sleepy.online_homeostasis(ppath, files)
#    laser = df.loc[df['laser'] == 'y']
#    iREMLaser = laser.iREM
#    noLaser = df.loc[df['laser'] == 'n']
#    noLaseriREM = noLaser.iREM.tolist()
#    pdb.set_trace()
#    pvalues = pg.ttest(iREMLaser, noLaseriREM)
    
    x = scipy.stats.ttest_rel(df[df.laser=='y']['iREM'], df[df.laser=='n']['iREM'])
    return x


#Script for determining at what timepoints NREM --> REM and NREM --> WAKE transitions become significantly different
def photometry_transitions_stats(ppath, files):
    #run zscored version
    trans_act, trans_act_trials, t, df = pyphi.activity_transitions(ppath, files, [(3,1), (1, 2), (3,2), (2,3)], 60, 60, [20, 60, 60], [20, 10, 0], ylim=(-1,2), pzscore=True)
#    pyphi.activity_transitions_old(ppath, files, [(3,1), (1, 2), (3,2), (2,3)], 60, 60, [20, 60, 60], [20, 0, 0], ylim=(-1,2), pzscore=True)
#    pdb.set_trace()
    
    pvaluesList = []
    
    for i in range(len(t)):
        ttest_indResult = stats.ttest_ind(trans_act_trials['NR'][:,i], trans_act_trials['NW'][:,i])
        pvalue = ttest_indResult[1]
#        pdb.set_trace()
        
        pvalueAdjust = pvalue / trans_act_trials['NR'].shape[1]
        pvaluesList.append(pvalueAdjust)
        
    labels = ['Time','pvalue','Significant']
    df = pd.DataFrame(columns = labels )
    significance = 'NO'
#    print('\n')
#    print('Significance values for NREM --> REM and NREM --> WAKE timepoint comparison')
    i=0
    for pvalue in pvaluesList:
        if pvalue < .05 / trans_act_trials['NR'].shape[1]:
#            print('Time: ', round(t[i], 8), ' pvalue: ', round(pvalue, 8), ' Significant: YES')
            significance = 'YES'
            dfTemp = pd.DataFrame(data = [[round(t[i], 1), pvalue, significance]], columns = labels)
            df = df.append(dfTemp)
        else:
#            print('Time: ', round(t[i], 8), ' pvalue: ', round(pvalue, 8), ' Significant: NO')
            significance = 'NO'
            dfTemp = pd.DataFrame(data = [[round(t[i], 1), pvalue, significance]], columns = labels)
            df = df.append(dfTemp)            
        i+=1

    return df

#Timecourse ANOVAS for dreadds, just use Franz's script instead
def timecourse_rm_anova(ppath, recordingFile):
#    condition = 'dreaddCNO'
    TimeMxCtr, TimeMxExp, df = sleepy.sleep_timecourse(ppath, trace_file = recordingFile, tbin = 1800, n = 10, stats = 'perc')
    df_rem = df[df.state == 'NREM']
    anova = pg.rm_anova(df_rem, dv='perc', within=['time', 'dose'], subject='mouse', correction=True)
    posthoc = pg.pairwise_ttests(data=df_rem, dv='perc', within=['time', 'dose'], subject='mouse', padjust='Bonferroni')
    pd.set_option('display.max_columns', None)
    print(posthoc)
#    pdb.set_trace()
    pd.set_option('display.max_columns', None)
    return anova
    
#Find instances of REM --> NREM in open loop dataset.
def rem_to_nrem_finder(ppath, recordings):
    for rec in recordings:
        M,K = load_stateidx(ppath, rec)
#        pdb.set_trace()
        for i in range(len(M)):
            if M[i] == 1:
                if M[i+1] == 3:
                    print(rec, ' has REM --> NREM transition')                    
                    break   
                
def wake_to_rem_finder(ppath, recordings):
    for rec in recordings:
        M,K = load_stateidx(ppath, rec)
        Mlen = range(len(M))
        for i in Mlen:
            if M[i] == 2 and i != Mlen[-1]:
                if M[i+1] == 1:
                    print(rec, ' has Wake --> REM transition')                    
                    # break 
                
# In the main opto dataset, count how many REM episodes are triggered by laser vs those that don't overlap with laser at all or start before the laser starts
def rem_with_laser_finder(ppath, recordings):
    totalREMBouts = 0
    REMBoutsWithLaser = 0
    REMBoutsNoLaser = 0
    for rec in recordings:
        M,K = load_stateidx(ppath, rec)
        
        # get indices for REM episodes
        REMBouts = sleepy.get_sequences(np.where(M==1)[0])
        totalREMBouts += len(REMBouts) #get total number of REM bouts
        
        laser = load_laser(ppath, rec)
        SR = 1000
        (istart, iend) = sleepy.laser_start_end(laser, SR)
        istart = [int(i/SR) for i in istart] #replace fbin with sr_eeg which is sampling rate of eeg (1000)
        iend = [int(i/SR) for i in iend]
        
        numBouts = len(REMBouts)
        numLaserBouts = 0
                
        for bout in REMBouts:
            bout = bout*2.5
            for (start,end) in zip(istart,iend):
                if bout[-1] < start: #if the laser period we're looking at is past the REM period we're looking at, just break out of the for loop
                    break
                laser = range(start,end)
                if any(np.isin(bout, laser)):
                    if bout[0] >= start: #Double check that the REM bout starts AFTER the laser onset, so we get true REM triggered episodes
                        numLaserBouts += 1 #update count for current rec for bugchecking purposes
                        REMBoutsWithLaser += 1 #if we get an overlap of the REM bout indices with the laser indices, update the rem bout with laser count and break out of the laser for loop
                        break    
                    # use the following if you just want any REM bouts overlapping with laser episodes and not specifically 'triggered' by the laser
                    # numLaserBouts += 1 #update count for current rec for bugchecking purposes
                    # REMBoutsWithLaser += 1 #if we get an overlap of the REM bout indices with the laser indices, update the rem bout with laser count and break out of the laser for loop
                    # break 
        numNoLaserBouts = numBouts - numLaserBouts
        REMBoutsNoLaser += numNoLaserBouts
        
    percentBoutsDuringLaser = (REMBoutsWithLaser / totalREMBouts)*100
        
    print('Number of REM bouts with no laser is: ' + str(REMBoutsNoLaser))
    print('Number of REM bouts starting after laser is: ' + str(REMBoutsWithLaser))    
    print('Percent of all REM bouts starting after laser onset is: ' + str(round(percentBoutsDuringLaser)))
        
    return REMBoutsWithLaser, REMBoutsNoLaser, percentBoutsDuringLaser

        


### FOR NRT STATS, SEE SEPARATE SCRIPT, NRT_STATS.PY ###



#basic results for closed loop modified to compare cellbody and retro, outputting a csv file for mixed anovas
def rem_online_analysis_dmm(ppath, recordings, backup='', single_mode=False, fig_file=''):
    """
    analyze results from closed-loop experiments
    :param ppath: base folder
    :param recordings: list of strings, recordinds
    :param backup: string, potential second backup folder with recordings
    :param single_mode: boolean, if True, average across all REM periods (irrespective of mouse)
           and plot each single REM period as dot
    :return: df, pd.DataFrame, with control and experimental REM durations as data columns
    """
    import pandas as pd

    if type(recordings) != list:
        recordings = [recordings]

    paths = dict()
    for rec in recordings:
        if os.path.isdir(os.path.join(ppath, rec)):
            paths[rec] = ppath
        else:
            paths[rec] = backup


    cellbody = ['JS68','JS69','JS156','JS158','JS159']
    retro = ['JS171','JS172','JS174','JS246','JS247','JS277']
    mice = dict()
    for rec in recordings:
        idf = re.split('_', rec)[0]
        if not idf in mice:
            mice[idf] = 1
    mice = list(mice.keys())
    if len(mice) == 1:
        single_mode=True

    dur_exp = {m:[] for m in mice}
    dur_ctr = {m:[] for m in mice}
    
    for rec in recordings:
        
        idf = re.split('_', rec)[0]
        M,S = load_stateidx(paths[rec], rec)
        sr = get_snr(paths[rec], rec)
        nbin = int(np.round(sr)*2.5)
        dt = (1.0/sr)*nbin

        laser = load_laser(paths[rec], rec)
        rem_trig = so.loadmat(os.path.join(paths[rec], rec, 'rem_trig_%s.mat'%rec), squeeze_me=True)['rem_trig']

        laser = sleepy.downsample_vec(laser, nbin)
        laser[np.where(laser>0)] = 1
        rem_trig = sleepy.downsample_vec(rem_trig, nbin)
        rem_trig[np.where(rem_trig>0)] = 1

        laser_idx = np.where(laser==1)[0]
        rem_idx = np.where(rem_trig==1)[0]

        # REM sequences from offline analysis (assumed to be the
        # "ground truth"
        seq = get_sequences(np.where(M==1)[0])
        for s in seq:
            # check true REM sequences overlapping with online detected sequences
            if len(np.intersect1d(s, rem_idx)) > 0:
                drn = (s[-1]-s[0]+1)*dt
                # does the sequence overlap with laser?
                if len(np.intersect1d(s, laser_idx))>0:
                    dur_exp[idf].append(drn)
                else:
                    dur_ctr[idf].append(drn)
    
    data = {'exp':[], 'ctr':[]}
    # if single_mode put all REM periods together,
    # otherwise average across REM periods for each mouse
    conditions = []
    
    if len(mice) == 1 or single_mode==True:
        for m in mice:
            data['exp'] += dur_exp[m]
            data['ctr'] += dur_ctr[m]                     
    else:
        mouseTotal = []
        mouseAvg = []
        for idf in dur_ctr:
            for x in np.array(dur_ctr[idf]):
                mouseTotal.append(x)
            for y in np.array(dur_exp[idf]):
                mouseTotal.append(y)               
            mouseAvg.append(sum(mouseTotal) / len(mouseTotal))
            mouseTotal = []
            
            dur_ctr[idf] = np.array(dur_ctr[idf]).mean()
            dur_exp[idf] = np.array(dur_exp[idf]).mean()
        data['exp'] = np.array(list(dur_exp.values()))
        data['ctr'] = np.array(list(dur_ctr.values()))
        
        totalMeanDur = sum(mouseAvg) / len(mouseAvg)
                
        for m in mice:
            if m in cellbody:
                conditions.append('cellbody')
            elif m in retro:
                conditions.append('retro')
            else:
                conditions.append('eYFP')
    
    df = pd.DataFrame({'ctr':pd.Series(data['ctr']), 'exp' : pd.Series(data['exp']), 'condition' : pd.Series(conditions)})

    # plot everything
    if not single_mode:
        plt.ion()
        plt.figure()
        ax = plt.axes([0.2, 0.15, 0.3, 0.7])
        df_mean = df.mean()
        plt.bar([1], [df_mean['ctr']], color='grey', label='W/o Laser')
        plt.bar([2], [df_mean['exp']], color='blue', label='With laser')
        plt.xticks([1,2])
        box_off(ax)
        #ax.set_xticklabels(['ctr', 'exp'], rotation=30)
        plt.ylabel('REM duration (s)')
        for (a,b,c) in zip(df['ctr'], df['exp'], df['condition']):
            if c == 'cellbody':
                plt.plot([1,2], [a,b], color='black')
            elif c == 'retro':
                plt.plot([1,2], [a,b], color='red')
            else:
                plt.plot([1,2], [a,b], color='green')
        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=1, frameon=False)
        ax.set_ylim(0,150)
    else:
        plt.figure()
        ax = plt.axes([0.2, 0.15, 0.3, 0.7])
        df_mean = df.mean()
        plt.bar([1], [df_mean['ctr']], color='grey')
        plt.bar([2], [df_mean['exp']], color='blue')
        plt.xticks([1,2])
        box_off(ax)
        #ax.set_xticklabels(['ctr', 'exp'], rotation=30)
        plt.ylabel('REM duration (s)')
        a = df['ctr']
        b = df['exp']
        plt.plot(np.ones((len(a),)), a, '.', color='black', label='W/o Laser')
        plt.plot(2*np.ones((len(b),)), b, '.', color='black', label='With laser')
        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=1, frameon=False)
    plt.show()
    
    # print('No laser avg REM duration ' + str(df_mean['ctr']))
    # print('With aser avg REM duration ' + str(df_mean['exp']))
    # print('Total mean REM duration ' + str(totalMeanDur))

    if len(fig_file) > 0:
        save_figure(fig_file)
        
    #now save the dataframe as a CSV file:
    df.to_csv('closed_loop_df.csv')
    
    # return df, mouseAvg
    return df



### RUN BEFORE CLOSED_LOOP_MIXED_ANOVA FUNCTION BELOW TO GENERATE CSV FILE #######
#basic results for closed loop modified to compare cellbody and retro
def rem_online_analysis_dmm_forANOVA(ppath, recordings, backup='', single_mode=False, fig_file=''):
    """
    analyze results from closed-loop experiments
    :param ppath: base folder
    :param recordings: a list of recordings including BOTH experiment and control CL datasets
    :param backup: string, potential second backup folder with recordings
    :param single_mode: boolean, if True, average across all REM periods (irrespective of mouse)
           and plot each single REM period as dot
    :return: df, pd.DataFrame, with control and experimental REM durations as data columns
    """
    import pandas as pd

    if type(recordings) != list:
        recordings = [recordings]

    paths = dict()
    for rec in recordings:
        if os.path.isdir(os.path.join(ppath, rec)):
            paths[rec] = ppath
        else:
            paths[rec] = backup


    cellbody = ['JS68','JS69','JS156','JS158','JS159','JS325','JS326','JS327','JS328']
    retro = ['JS171','JS172','JS174','JS246','JS247','JS277']
    swichr = ['JS306','JS307','JS336','JS337','JS338','JS339','JS340','JS341','JS342']
    #anything else not in this list is considered eyfp
    
    mice = dict()
    for rec in recordings:
        idf = re.split('_', rec)[0]
        if not idf in mice:
            mice[idf] = 1
    mice = list(mice.keys())
    if len(mice) == 1:
        single_mode=True

    dur_exp = {m:[] for m in mice}
    dur_ctr = {m:[] for m in mice}
    
    for rec in recordings:
        print(rec)
        
        idf = re.split('_', rec)[0]
        M,S = load_stateidx(paths[rec], rec)
        sr = get_snr(paths[rec], rec)
        nbin = int(np.round(sr)*2.5)
        dt = (1.0/sr)*nbin

        laser = load_laser(paths[rec], rec)
        rem_trig = so.loadmat(os.path.join(paths[rec], rec, 'rem_trig_%s.mat'%rec), squeeze_me=True)['rem_trig']

        laser = sleepy.downsample_vec(laser, nbin)
        laser[np.where(laser>0)] = 1
        rem_trig = sleepy.downsample_vec(rem_trig, nbin)
        rem_trig[np.where(rem_trig>0)] = 1

        laser_idx = np.where(laser==1)[0]
        rem_idx = np.where(rem_trig==1)[0]

        # REM sequences from offline analysis (assumed to be the
        # "ground truth"
        seq = get_sequences(np.where(M==1)[0])
        for s in seq:
            # check true REM sequences overlapping with online detected sequences
            if len(np.intersect1d(s, rem_idx)) > 0:
                drn = (s[-1]-s[0]+1)*dt
                # does the sequence overlap with laser?
                if len(np.intersect1d(s, laser_idx))>0:
                    dur_exp[idf].append(drn)
                else:
                    dur_ctr[idf].append(drn)
                  
    
    data = {'exp':[], 'ctr':[]}
    
    # if single_mode put all REM periods together,
    # otherwise average across REM periods for each mouse    
    conditions = []
    if len(mice) == 1 or single_mode==True:
        for m in mice:
            data['exp'] += dur_exp[m]
            data['ctr'] += dur_ctr[m]                       
    else: #Single mode = False, default
        df_For_Anova = pd.DataFrame()
        labels = ['Name','Duration','Virus','Laser']
       
        mouseTotal = []
        mouseAvg = []
        for idf in dur_ctr:
            for x in np.array(dur_ctr[idf]):
                mouseTotal.append(x)
            for y in np.array(dur_exp[idf]):
                mouseTotal.append(y)               
            mouseAvg.append(sum(mouseTotal) / len(mouseTotal))
            mouseTotal = []
            
            dur_ctr[idf] = np.array(dur_ctr[idf]).mean()
            dur_exp[idf] = np.array(dur_exp[idf]).mean()
            
        data['exp'] = np.array(list(dur_exp.values()))
        data['ctr'] = np.array(list(dur_ctr.values()))
        
        #Build the dataframe
        for (key,value) in dur_ctr.items():
            if key in cellbody:                
                tempData = [[key, value, 'AAV1', 'W/o Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            elif key in retro:
                tempData = [[key, value, 'rgAAV', 'W/o Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            elif key in swichr:
                tempData = [[key, value, 'swichr', 'W/o Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            else: 
                tempData = [[key, value, 'eYFP', 'W/o Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
                
        for (key,value) in dur_exp.items():
            if key in cellbody:                
                tempData = [[key, value, 'AAV1', 'Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            elif key in retro:
                tempData = [[key, value, 'rgAAV', 'Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            elif key in swichr:
                tempData = [[key, value, 'swichr', 'Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
            else: 
                tempData = [[key, value, 'eYFP', 'Laser']]
                tempDF = pd.DataFrame(data = tempData, columns = labels)
                df_For_Anova = df_For_Anova.append(tempDF)
        
    #now save the dataframe as a CSV file:
    df_For_Anova.to_csv('closed_loop_df.csv')
    # df_For_Anova.to_csv('closed_loop_retro_df.csv')
    # df_For_Anova.to_csv('closed_loop_swichr_df.csv')
    
    # return df, mouseAvg
    return df_For_Anova


# MIXED ANOVA FOR CLOSED LOOP RESULTS, read from the csv file generated by the function above, rem_online_analysis_dmm_forANOVA:
# THE SAME THING FOR THE INTER-REM HOMEOSTASIS CAN BE FOUND IN closed_loop_homeostasis_and_delta.py
def closed_loop_mixed_ANOVA():
    # mode = 'ChR2'
    # mode = 'retro'
    mode = 'swichr'
    
    df = []
    if mode == 'ChR2':
        df = pd.read_csv('closed_loop_df.csv') #Path is the script folder, which is where rem_online_analysis_dmm_forANOVA saves it.
    elif mode == 'retro':    
        df = pd.read_csv('closed_loop_retro_df.csv') #.csv file for retrograde closed loop
    elif mode == 'swichr':
        df = pd.read_csv('closed_loop_swichr_df.csv')
    
    # pdb.set_trace()
    if mode == 'swichr':
        df = df.loc[df['Name'] != 'JS343'] #eyfp in swichr dataset low duration outlier
    
    df = df.loc[df['Name'] != 'JS340'] #swichr low duration outlier
    # df = df.loc[df['Name'] != 'JS342'] #swichr increase
    
    print('Two-way mixed anova with virus and laser as between and within factors, ' + mode + ':')
    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    
    ANOVA = pg.mixed_anova(data = df, subject = 'Name', dv = 'Duration', within='Laser', between='Virus')    
    print(ANOVA)

    ttests = pg.pairwise_ttests(data = df, subject = 'Name', dv = 'Duration', within='Laser', between='Virus', within_first=True, padjust='bonf')
    print('\n')
    print('Pairwise ttests ' + mode)
    print(ttests)
    
    
    if mode == 'ChR2':
        ChR2LaserDF = df[(df.Virus == 'AAV1') & (df.Laser == 'Laser')]
        ChR2Laser = ChR2LaserDF['Duration'].tolist()
        ChR2NoLaserDF = df[(df.Virus == 'AAV1') & (df.Laser == 'W/o Laser')]
        ChR2NoLaser = ChR2NoLaserDF['Duration'].tolist()
        ChR2Pairedttest = pg.ttest(x = ChR2Laser, y = ChR2NoLaser, paired=True)
        print('\n')
        # pdb.set_trace()
        ChR2Laser = np.asarray(ChR2Laser)
        ChR2NoLaser = np.asarray(ChR2NoLaser)
        ChR2LaserMean = ChR2Laser.mean()
        ChR2NoLaserMean = ChR2NoLaser.mean()
        ChR2MeanDiff = ChR2LaserMean - ChR2NoLaserMean
        cil = ChR2Pairedttest['CI95%'][0][0]   
        cih = ChR2Pairedttest['CI95%'][0][1]   
        print('Difference of the means (ChR2 laser vs no laser) is ' + str(ChR2MeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        print('ChR2 Laser paired ttest:')
        print(ChR2Pairedttest) 
    
    elif mode == 'retro':
        retroLaserDF = df[(df.Virus == 'rgAAV') & (df.Laser == 'Laser')]
        retroLaser = retroLaserDF['Duration'].tolist()
        retroNoLaserDF = df[(df.Virus == 'rgAAV') & (df.Laser == 'W/o Laser')]
        retroNoLaser = retroNoLaserDF['Duration'].tolist()
        retroPairedttest = pg.ttest(x = retroLaser, y = retroNoLaser, paired=True)
        print('\n')
        # pdb.set_trace()
        retroLaser = np.asarray(retroLaser)
        retroNoLaser = np.asarray(retroNoLaser)
        retroLaserMean = retroLaser.mean()
        retroNoLaserMean = retroNoLaser.mean()
        retroMeanDiff = retroLaserMean - retroNoLaserMean
        cil = retroPairedttest['CI95%'][0][0]   
        cih = retroPairedttest['CI95%'][0][1]   
        print('Difference of the means (retro laser vs no laser) is ' + str(retroMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        print('retro Laser paired ttest:')
        print(retroPairedttest) 
    
    elif mode == 'swichr':
        swichrLaserDF = df[(df.Virus == 'swichr') & (df.Laser == 'Laser')]
        swichrLaser = swichrLaserDF['Duration'].tolist()
        swichrNoLaserDF = df[(df.Virus == 'swichr') & (df.Laser == 'W/o Laser')]
        swichrNoLaser = swichrNoLaserDF['Duration'].tolist()
        swichrPairedttest = pg.ttest(x = swichrLaser, y = swichrNoLaser, paired=True)
        swichrLaser = np.asarray(swichrLaser)
        swichrNoLaser = np.asarray(swichrNoLaser)
        swichrLaserMean = swichrLaser.mean()
        swichrNoLaserMean = swichrNoLaser.mean()
        swichrMeanDiff = swichrLaserMean - swichrNoLaserMean 
        cil = swichrPairedttest['CI95%'][0][0]   
        cih = swichrPairedttest['CI95%'][0][1]  
        print('\n')
        print('Difference of the means (swichr laser vs no laser) is ' + str(swichrMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        print('Swichr Laser paired ttest:')
        print(swichrPairedttest)
    
    
    eYFPLaserDF = df[(df.Virus == 'eYFP') & (df.Laser == 'Laser')]
    eYFPLaser = eYFPLaserDF['Duration'].tolist()
    eYFPNoLaserDF = df[(df.Virus == 'eYFP') & (df.Laser == 'W/o Laser')]
    eYFPNoLaser = eYFPNoLaserDF['Duration'].tolist()
    eYFPPairedttest = pg.ttest(x = eYFPLaser, y = eYFPNoLaser, paired=True)
    print('\n')
    eYFPLaser = np.asarray(eYFPLaser)
    eYFPNoLaser = np.asarray(eYFPNoLaser)
    eYFPLaserMean = eYFPLaser.mean()
    eYFPNoLaserMean = eYFPNoLaser.mean()
    eYFPMeanDiff = eYFPLaserMean - eYFPNoLaserMean 
    cil = eYFPPairedttest['CI95%'][0][0]   
    cih = eYFPPairedttest['CI95%'][0][1]    
    print('Difference of the means (eYFP laser vs no laser) is ' + str(eYFPMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    print('eYFP Laser paired ttest:')
    print(eYFPPairedttest) 
    

    #Now get unpaired ttest means and confidence intervals, which will then be paired with the p-values from the pairwise ttests from the ANOVA
    if mode == 'ChR2':    
        print('\n')
        print('95CIs from unpaired ttest to be used with Mixed ANOVA pairwise ttest comparisons above')
        laserMeanDiff = ChR2LaserMean - eYFPLaserMean
        laserUnpairedttest = pg.ttest(x = ChR2Laser, y = eYFPLaser, paired = False)
        cil = laserUnpairedttest['CI95%'][0][0]   
        cih = laserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (With laser, ChR2 vs eYFP) is ' + str(laserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        
        noLaserMeanDiff = ChR2NoLaserMean - eYFPNoLaserMean
        noLaserUnpairedttest = pg.ttest(x = ChR2NoLaser, y = eYFPNoLaser, paired = False)
        cil = noLaserUnpairedttest['CI95%'][0][0]   
        cih = noLaserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (No laser, ChR2 vs eYFP) is ' + str(noLaserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    elif mode == 'retro':
        print('\n')
        print('95CIs from unpaired ttest to be used with Mixed ANOVA pairwise ttest comparisons above')
        laserMeanDiff = retroLaserMean - eYFPLaserMean
        laserUnpairedttest = pg.ttest(x = retroLaser, y = eYFPLaser, paired = False)
        cil = laserUnpairedttest['CI95%'][0][0]   
        cih = laserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (With laser, Retro vs eYFP) is ' + str(laserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        
        noLaserMeanDiff = retroNoLaserMean - eYFPNoLaserMean
        noLaserUnpairedttest = pg.ttest(x = retroNoLaser, y = eYFPNoLaser, paired = False)
        cil = noLaserUnpairedttest['CI95%'][0][0]   
        cih = noLaserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (No laser, Retro vs eYFP) is ' + str(noLaserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    elif mode == 'swichr':
        print('\n')
        print('95CIs from unpaired ttest to be used with Mixed ANOVA pairwise ttest comparisons above')
        laserMeanDiff = swichrLaserMean - eYFPLaserMean
        laserUnpairedttest = pg.ttest(x = swichrLaser, y = eYFPLaser, paired = False)
        cil = laserUnpairedttest['CI95%'][0][0]   
        cih = laserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (With laser, Swichr vs eYFP) is ' + str(laserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
        
        noLaserMeanDiff = swichrNoLaserMean - eYFPNoLaserMean
        noLaserUnpairedttest = pg.ttest(x = swichrNoLaser, y = eYFPNoLaser, paired = False)
        cil = noLaserUnpairedttest['CI95%'][0][0]   
        cih = noLaserUnpairedttest['CI95%'][0][1]
        print('Difference of the means (No laser, Swichr vs eYFP) is ' + str(noLaserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    
    
    ##### PLOT #####
    clrs = ['gray','blue']
    
    fig, axes = plt.subplots(1, 2, sharex = True, sharey = True) #, figsize = figDimensionsBe sure to set figsize to something that looks good      
    numPlot = 0 
    clrs = ['gray','blue']
    linecolor = ('black', 'black','black', 'black','black', 'black','black', 'black','black', 'black','black', 'black')
    
    virus_list = []
    if mode == 'ChR2':
        virus_list = ['AAV1','eYFP']
    elif mode == 'retro':
        virus_list = ['rgAAV','eYFP']
    elif mode == 'swichr':
        virus_list = ['swichr','eYFP']
    
    for virus in virus_list:
    
        dataToUse = df.loc[df['Virus'] == virus]
        sns.barplot(data = dataToUse, y = 'Duration', x = 'Laser', palette = clrs, ax=axes[numPlot])
        
        condition1 = df.loc[df['Virus'] == virus]
        condition1 = condition1.loc[condition1['Laser'] == 'W/o Laser']
        condition1Stat = condition1['Duration'].tolist()
        
        condition2 = df.loc[df['Virus'] == virus]
        condition2 = condition2.loc[condition2['Laser'] == 'Laser']
        condition2Stat = condition2['Duration'].tolist()    
        
        sns.lineplot(data = dataToUse, y='Duration', x='Laser', hue='Name', ax=axes[numPlot], palette=linecolor, legend=False)            
        
        sns.despine()
        axes[numPlot].set(xlabel = virus)
        if numPlot != 0:
            axes[numPlot].get_yaxis().set_visible(False) #get rid of y-axis ticks
            sns.despine(left = True, ax = axes[numPlot])      
        fig.suptitle('Closed loop laser protocol for REM sleep')
        if mode == 'ChR2':
            plt.ylim(0,200)  
        if mode == 'retro':
            plt.ylim(0,150) 
        if mode == 'swichr':
            plt.ylim(0,130)                               
        plt.tight_layout() 
        plt.show()
        numPlot += 1



def closed_loop_mixed_ANOVA_AAV1_rgAAV_comparison():

    dfChR2 = []
    dfRetro = []

    dfChR2 = pd.read_csv('closed_loop_df.csv')
    dfRetro = pd.read_csv('closed_loop_retro_df.csv')
    # pdb.set_trace()
    #put them together without eyfp
    dfTogether = pd.DataFrame()
    dfTogether = dfTogether.append(dfChR2[dfChR2.Virus == 'AAV1'])
    dfTogether = dfTogether.append(dfRetro[dfRetro.Virus == 'rgAAV'])
    
    #calculate differences
    dfToUse = dfChR2[dfChR2.Virus == 'AAV1']
    for mouse in dfToUse.Name.unique():
        noLaserAAV1 = dfToUse[(dfToUse.Name == mouse) & (dfToUse.Virus == 'AAV1') & (dfToUse.Laser == 'W/o Laser')]['Duration'].tolist()[0]
        laserAAV1 = dfToUse[(dfToUse.Name == mouse) & (dfToUse.Virus == 'AAV1') & (dfToUse.Laser == 'Laser')]['Duration'].tolist()[0]
        diffAAV1 = laserAAV1 - noLaserAAV1
        tempDF1 = pd.DataFrame(data = [[mouse, diffAAV1, 'AAV1', 'Diff']], columns = ['Name','Duration','Virus','Laser'])
        dfTogether = dfTogether.append(tempDF1)
        
    dfToUse = dfRetro[dfRetro.Virus == 'rgAAV']
    for mouse in dfToUse.Name.unique():
        noLaserRetro = dfToUse[(dfToUse.Name == mouse) & (dfToUse.Virus == 'rgAAV') & (dfToUse.Laser == 'W/o Laser')]['Duration'].tolist()[0]
        laserRetro = dfToUse[(dfToUse.Name == mouse) & (dfToUse.Virus == 'rgAAV') & (dfToUse.Laser == 'Laser')]['Duration'].tolist()[0]
        diffRetro = laserRetro - noLaserRetro
        tempDF2 = pd.DataFrame(data = [[mouse, diffRetro, 'rgAAV', 'Diff']], columns = ['Name','Duration','Virus','Laser'])
        dfTogether = dfTogether.append(tempDF2)

    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    
        
    ChR2LaserDF = dfChR2[(dfChR2.Virus == 'AAV1') & (dfChR2.Laser == 'Laser')]
    ChR2Laser = ChR2LaserDF['Duration'].tolist()
    ChR2NoLaserDF = dfChR2[(dfChR2.Virus == 'AAV1') & (dfChR2.Laser == 'W/o Laser')]
    ChR2NoLaser = ChR2NoLaserDF['Duration'].tolist()
    ChR2Pairedttest = pg.ttest(x = ChR2Laser, y = ChR2NoLaser, paired=True)
    

    # ChR2MeanDiffDF = pd.DataFrame(data = )
    print('\n')
    # pdb.set_trace()
    ChR2Laser = np.asarray(ChR2Laser)
    ChR2NoLaser = np.asarray(ChR2NoLaser)
    ChR2LaserMean = ChR2Laser.mean()
    ChR2NoLaserMean = ChR2NoLaser.mean()
    ChR2DiffArray = ChR2Laser - ChR2NoLaser
    ChR2MeanDiff = ChR2LaserMean - ChR2NoLaserMean
    cil = ChR2Pairedttest['CI95%'][0][0]   
    cih = ChR2Pairedttest['CI95%'][0][1]   
    print('Difference of the means (ChR2 laser vs no laser) is ' + str(ChR2MeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    print('ChR2 Laser paired ttest:')
    print(ChR2Pairedttest) 
    
    
    retroLaserDF = dfRetro[(dfRetro.Virus == 'rgAAV') & (dfRetro.Laser == 'Laser')]
    retroLaser = retroLaserDF['Duration'].tolist()
    retroNoLaserDF = dfRetro[(dfRetro.Virus == 'rgAAV') & (dfRetro.Laser == 'W/o Laser')]
    retroNoLaser = retroNoLaserDF['Duration'].tolist()
    retroPairedttest = pg.ttest(x = retroLaser, y = retroNoLaser, paired=True)
    print('\n')
    retroLaser = np.asarray(retroLaser)
    retroNoLaser = np.asarray(retroNoLaser)
    retroDiffValues = retroLaser - retroNoLaser    
    retroCIs = pg.compute_bootci(x = retroDiffValues, func = 'mean')
    retroLaserMean = retroLaser.mean()
    retroNoLaserMean = retroNoLaser.mean()
    retroMeanDiff = retroLaserMean - retroNoLaserMean
    cil = retroPairedttest['CI95%'][0][0]   
    cih = retroPairedttest['CI95%'][0][1]   
    print('Difference of the means (retro laser vs no laser) is ' + str(retroMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    print('Retro laser - noLaser compute_bootci is ')
    print(retroCIs)
    print('retro Laser paired ttest:')
    print(retroPairedttest) 
    
    print('ChR2 vs rgAAV Two-way mixed anova with virus and laser as between and within factors:')
    ANOVA = pg.mixed_anova(data = dfTogether, subject = 'Name', dv = 'Duration', within='Laser', between='Virus')    
    print(ANOVA)

    ttests = pg.pairwise_ttests(data = dfTogether, subject = 'Name', dv = 'Duration', within='Laser', between='Virus', within_first=True, padjust='bonf')
    print('\n')
    print('Pairwise ttests')
    print(ttests)
    
    
    print('\n')
    print('95CIs from unpaired ttest to be used with Mixed ANOVA pairwise ttest comparisons above')
    laserMeanDiff = ChR2LaserMean - retroLaserMean
    laserUnpairedttest = pg.ttest(x = ChR2Laser, y = retroLaser, paired = False)
    cil = laserUnpairedttest['CI95%'][0][0]   
    cih = laserUnpairedttest['CI95%'][0][1]
    print('Difference of the means (With laser, ChR2 vs retro) is ' + str(laserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    
    noLaserMeanDiff = ChR2NoLaserMean - retroNoLaserMean
    noLaserUnpairedttest = pg.ttest(x = ChR2NoLaser, y = retroNoLaser, paired = False)
    cil = noLaserUnpairedttest['CI95%'][0][0]   
    cih = noLaserUnpairedttest['CI95%'][0][1]
    print('Difference of the means (No laser, ChR2 vs retro) is ' + str(noLaserMeanDiff) + ', with 95CI = ' + ' (' + str(cil) + ',' + str(cih) + ')')
    

    #paired mean differences:
    print('\n')
    AAV1Diff = dfTogether[(dfTogether.Laser == 'Diff') & (dfTogether.Virus == 'AAV1')]['Duration'].tolist()
    rgAAVDiff = dfTogether[(dfTogether.Laser == 'Diff') & (dfTogether.Virus == 'rgAAV')]['Duration'].tolist()
    # pdb.set_trace()
    AAV1DiffMean = np.asarray(AAV1Diff).mean()
    rgAAVDiffMean = np.asarray(rgAAVDiff).mean()
    print('Paired mean differences (Laser - noLaser) AAV1 - rgAAv unpaired ttest: ' + str(AAV1DiffMean - rgAAVDiffMean))
    unpairedttest = pg.ttest(x = AAV1Diff, y = rgAAVDiff, paired = 'False')
    print(unpairedttest)
    print('\n')
    
    
    
    
    
    
    ##### PLOT #####
    clrs = ['gray','blue']
    
    fig, axes = plt.subplots(1, 2, sharex = True, sharey = True) #, figsize = figDimensionsBe sure to set figsize to something that looks good      
    numPlot = 0 
    clrs = ['gray','blue']
    linecolor = ('black', 'black','black', 'black','black', 'black','black', 'black','black', 'black','black', 'black')
    
    virus_list = ['AAV1', 'rgAAV']
    
    for virus in virus_list:
    
        dataToUse = dfTogether[(dfTogether.Virus == virus) & (dfTogether.Laser != 'Diff')]
        sns.barplot(data = dataToUse, y = 'Duration', x = 'Laser', palette = clrs, ax=axes[numPlot])
        
        condition1 = dfTogether.loc[dfTogether['Virus'] == virus]
        condition1 = condition1.loc[condition1['Laser'] == 'W/o Laser']
        condition1Stat = condition1['Duration'].tolist()
        
        condition2 = dfTogether.loc[dfTogether['Virus'] == virus]
        condition2 = condition2.loc[condition2['Laser'] == 'Laser']
        condition2Stat = condition2['Duration'].tolist()    
        
        sns.lineplot(data = dataToUse, y='Duration', x='Laser', hue='Name', ax=axes[numPlot], palette=linecolor, legend=False)            
        
        sns.despine()
        axes[numPlot].set(xlabel = virus)
        if numPlot != 0:
            axes[numPlot].get_yaxis().set_visible(False) #get rid of y-axis ticks
            sns.despine(left = True, ax = axes[numPlot])      
        fig.suptitle('Closed loop laser protocol for REM sleep')
        plt.ylim(0,200)    
        # plt.ylim(0,150)                                
        plt.tight_layout() 
        plt.show()
        numPlot += 1
        
    #paired mean differences plot
    plt.figure()
    ax = plt.axes([0.15, 0.1, 0.35, 0.4])
    df_sel = dfTogether[dfTogether.Virus.isin(['AAV1', 'rgAAV'])]
    df_sel = df_sel[df_sel.Laser == 'Diff']
    # pdb.set_trace()
    
    sns.barplot(data=df_sel, x = 'Virus', y='Duration', palette={'AAV1':'cornflowerblue', 'rgAAV':'firebrick'}, saturation=0.8)
    g = sns.stripplot(data=df_sel, x = 'Virus', y='Duration', dodge=True, size=4, palette={'AAV1':[0.3]*3, 'rgAAV':[0.3]*3})
    g.legend().remove()
    plt.ylim([0,150])
    plt.ylabel('$\Delta$ Duration (s)')
    plt.plot([-0.5, 2.5], [0,0], 'k', lw=0.8)
    plt.xlabel(None)
    sns.despine()
    plt.xlim([-0.5, 2.5])








#Very basic t-test_ind for eYFP and ChR2 closed loop durations -- USE ANOVA FUNCTIONS INSTEAD
def closed_loop_eYFPvsChR2(x):
    y=x
    ppath1 = r'D:\WeberLabMaterialsJoe\SleepLibraryFunctions\AllMyRecordings\GAD2 dPGi\dPGi_control_eYFP_laser\closed loop' 
    files1 = ['JS68_050819n1','JS68_051219n1','JS69_061319n1','JS156_111119n1','JS156_111419n1','JS156_112019n1','JS158_111319n1','JS158_111519n1','JS159_111319n1','JS159_111419n1','JS159_111519n1']
    ppath2 = r'D:\WeberLabMaterialsJoe\SleepLibraryFunctions\AllMyRecordings\GAD2 dPGi\dPGi_control_eYFP_laser\closed loop' 
    files2 = ['JS250_081320n1','JS250_081420n1','JS251_081320n1','JS251_081420n1','JS252_081320n1','JS252_081420n1','JS253_081320n1','JS253_081420n1','JS272_100120n1','JS272_100220n1','JS273_100120n1','JS273_100220n1']

    print('Running cell body data')
    df1, cellbody = rem_online_analysis(ppath1, files1)
    print('Running eYFP data')
    df2, retro = rem_online_analysis(ppath1, files2)
    
    x = stats.ttest_ind(cellbody, retro)
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return x


#Plots faield transitions and successful transitions as separate
def laser_brainstate_bootstrap_transitions(ppath, recordings, pre, post, edge=0, sf=0, nboots=10000, alpha=0.05, backup='',
                               start_time=0, ma_thr=10, bootstrap_mode=0, fig_file=''):
    """
    Align brain state with laser stimulation and calculate two-sided 1-$alpha confidence intervals using
    bootstrapping.
    :param ppath: base folder
    :param recordings: list of recordings
    :param pre: time before laser
    :param post: time after laser onset
    :param edge: add $edge seconds add beginning and end (that are not shown in the plot) to avoid filtering artifacts
    :param sf: smoothing factor for Gaussian filter; better do not use
    :param nboots: int, how many times the whole data set is resampled for boot-strapping
    :param alpha: plot shows 1-$alpha confidence interval
    :param backup: optional backup folder where recordings are stored
    :param start_time: start time of recordding used for analysis
    :param ma_thr: sleep periods < ma_thr are thrown away
    :param bootstrap_mode: default=0
           bootstrap_mode == 0: Take inter-mouse variance and inter-trial variance (of each mouse) into account.
           That is, bootstrapping re-models the variance expected when re-doing the same
           experimental design (same mouse number and total trial number).
           To account for potentially different number of trials per mouse, resample the data
           during each iteration the following way: Assume that there are n laser trials from m mice;
           randomly select (with replacment) ntrial mice; then select from each mouse randomly one trial.
           bootstrap_mode == 1: Only take inter-trial variance (of each mouse) into account. That is,
           bootstrapping models the variance expected when redoing the experiment with exactly the same mice.
    :param fig_file, if file name is specified, the figure will be saved
    :return: P   - p-values for NREM, REM, Wake
             Mod - by how much the percentage of NREM, REM, Wake is increased compared to baseline
    """

    pre += edge
    post += edge

    rec_paths = dict()
    for rec in recordings:
        if os.path.isdir(os.path.join(ppath, rec)):
            rec_paths[rec] = ppath
        else:
            rec_paths[rec] = backup

    # dict: mouse_id --> laser trials, R W N sequence
    BrainstateDict = {}
    for rec in recordings:
        idf = re.split('_', rec)[0]
        BrainstateDict[idf] = []
    mice = list(BrainstateDict.keys())
    nmice = len(BrainstateDict)

    for rec in recordings:
        ppath = rec_paths[rec]
        SR = get_snr(ppath, rec)
        NBIN = np.round(2.5 * SR)
        dt = NBIN * 1 / SR
        istart_time = int(np.round(start_time / dt))

        M = load_stateidx(ppath, rec)[0]

        seq = get_sequences(np.where(M == 2)[0])
        for s in seq:
            if len(s) * dt <= ma_thr:
                M[s] = 3

        (idxs, idxe) = laser_start_end(load_laser(ppath, rec))
        idf = re.split('_', rec)[0]

        SR = get_snr(ppath, rec)
        NBIN = np.round(2.5 * SR)
        ipre = int(np.round(pre / dt))
        ipost = int(np.round(post / dt))

        idxs = [int(i / NBIN) for i in idxs]
        idxe = [int(i / NBIN) for i in idxe]
        laser_dur = np.mean((np.array(idxe) - np.array(idxs))) * dt

        for (i, j) in zip(idxs, idxe):
            if i >= ipre and j + ipost <= len(M) - 1 and i > istart_time:
                bs = M[i - ipre:i + ipost + 1]
                BrainstateDict[idf].append(bs)

    for mouse in mice:
        BrainstateDict[mouse] = np.array(BrainstateDict[mouse])

    # I assume here that every recording has same dt
    t = np.linspace(-ipre*dt, ipost*dt, ipre+ipost+1)
    Trials = dict()
    for mouse in BrainstateDict:
        Trials[mouse] = np.zeros((BrainstateDict[mouse].shape[0], len(t), 5)) #HERE

    # total number of trials:
    ntrials = 0
    for mouse in BrainstateDict:
        M = np.array(BrainstateDict[mouse])
        for state in range(1, 6): #HERE IMPORTANT CHANGE, went from 1,4 to 1,6
            C = np.zeros(M.shape)
            C[np.where(M == state)] = 100.
            Trials[mouse][:,:,state-1] = C
        ntrials += Trials[mouse].shape[0]

    Prob = np.zeros((nboots, len(t), 5)) #HERE
    if bootstrap_mode == 1:
        for b in range(nboots):
            # average brain state percentage for each mouse during iteration b
            mouse_mean_state = np.zeros((nmice, len(t), 5)) #HERE
            i = 0
            for mouse in mice:
                mmouse = Trials[mouse].shape[0]
                iselect = rand.randint(0, mmouse, (mmouse,))
                for s in [1,2,3,4,5]: #HERE
                    #bBS[s][offset:offset+mmouse,:] = Trials[mouse][iselect,:,s-1]
                    mouse_mean_state[i,:,s-1] = Trials[mouse][iselect,:,s-1].mean(axis=0)
                i += 1

            for s in [1,2,3,4,5]: #HERE
                Prob[b,:,s-1] = mouse_mean_state[:,:,s-1].mean(axis=0)
    else:
        mx_iter = np.zeros((ntrials, len(t), 5)) #HERE
        for b in range(nboots):
            # for each iteration select randomly select ntrials mice
            irand_mice = rand.randint(0, nmice, ntrials)

            # average brain state percentage for each mouse during iteration b
            # mouse_mean_state = np.zeros((nmice, len(t), 3))
            i = 0
            # there are ntrials mice
            for j in irand_mice:
                mouse = mice[irand_mice[j]]
                # we have nmouse trials per mouse
                mmouse = Trials[mouse].shape[0]
                # select one random trial from the current mouse
                iselect = rand.randint(0, mmouse)
                for s in [1, 2, 3,4,5]: #HERE
                    mx_iter[i,:,s-1] = Trials[mouse][iselect,:,s-1]
                i += 1

            # mx_iter is the resampled data set for bootstrap iteration b
            # now we calculate the statistics we're interesting in, which is the mean
            for s in [1, 2, 3,4,5]: #HERE
                Prob[b,:,s-1] = mx_iter[:,:,s-1].mean(axis=0)


    # simple average for each brainstate across mice (w/o) bootstrapping
    Prob_mean = np.zeros((nmice, len(t), 5)) #HERE
    for s in [1,2,3,4,5]: #HERE
        i = 0
        for mouse in mice:
            Prob_mean[i,:,s-1] = Trials[mouse][:,:,s-1].mean(axis=0)
            i += 1

    # pdb.set_trace()
    
    usProb = Prob.copy()
    Prob = np.sort(Prob, axis=0)
    Bounds = np.zeros((2, len(t), 5)) #HERE
    a = int((nboots * alpha) / 2.0)
    for s in [1,2,3,4,5]:#HERE
        Bounds[0,:,s-1] = Prob[a,:,s-1]
        Bounds[1,:,s-1] = Prob[-a,:, s-1]

    # smooth_data
    if sf > 0:
        for s in range(5): #HERE
            Bounds[0, :, s] = smooth_data(Bounds[0, :, s], sf)
            Bounds[1, :, s] = smooth_data(Bounds[1, :, s], sf)
        for i in range(nmice):
            for s in range(5): #HERE
                Prob_mean[i, :, s] = smooth_data(Prob_mean[i,:,s], sf)

    # plot figure
    colors = np.array([[0, 1, 1], [0.5, 0, 1], [0.6, 0.6, 0.6],[0, 0, .5], [1, 0, 0]]) #HERE
    br_states = {1:'REM', 2:'Wake', 3:'NREM', 4:'Trans', 5:'Failed_Trans'} #HERE
    #colors = np.array([[55,255,255], [153,255,153],[153,153,153]])/255.
    it = np.where((t>=-pre+edge) & (t<=post-edge))[0]
    plt.ion()
    plt.figure()
    ax = plt.axes([0.15, 0.15, 0.6, 0.7])
    ax.add_patch(patches.Rectangle((0, 0), laser_dur, 100, facecolor=[0.6, 0.6, 1], edgecolor=[0.6, 0.6, 1]))

    for state in [5,4,3,2,1]: #HERE
        # pdb.set_trace()
        if state == 5:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=5, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=5,)
        elif state == 4:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=4, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=4)
        elif state == 3:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=1, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=1)
        elif state == 2:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=2, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=2)
        elif state == 1:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=3, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=3)
        
        # ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state]) 
        
    
    plt.xlim([-pre+edge, post-edge])
    plt.ylim([0,100])
    plt.xlabel('Time (s)')
    plt.ylabel('Brain state (%)')
    plt.legend(bbox_to_anchor = (1.0, 0.7, 1., .102), loc = 5, mode = 'expand', ncol = 1, frameon = False) #HERE
    box_off(ax)
    plt.draw()

    # statistics
    ibase = np.where((t>=-laser_dur) & (t<0))[0]
    ilsr  = np.where((t>=0) & (t<laser_dur))[0]
    P   = np.zeros((5,)) #HERE
    Mod = np.zeros((5,)) #HERE
    for istate in [1,2,3,4,5]: #HERE
        basel = usProb[:,ibase,istate-1].mean(axis=1)
        laser = usProb[:,ilsr, istate-1].mean(axis=1)
        d = laser - basel
        if np.mean(d) >= 0:
            # now we want all values be larger than 0
            p = len(np.where(d>0)[0]) / (1.0*nboots)
            sig = 1 - p
            if sig == 0:
                sig = 1.0/nboots
            Mod[istate-1] = (np.mean(laser) / np.mean(basel) - 1) * 100
        else:
            p = len(np.where(d<0)[0]) / (1.0*nboots)
            sig = 1 - p
            if sig == 0:
                sig = 1.0/nboots
            Mod[istate-1] = -(1 - np.mean(laser) / np.mean(basel)) * 100
        P[istate-1] = sig

    labels = {1:'REM', 2:'Wake', 3:'NREM', 4:'Trans', 5:'Failed_Trans'} #HERE
    for s in [1,2,3,4,5]: #HERE
        print('%s is changed by %f perc.; P = %f, bootstrap' % (labels[s], Mod[s-1], P[s-1]))
    print("n = %d mice" % len(mice))

    if len(fig_file) > 0:
        plt.savefig(fig_file, bbox_inches="tight")

    return P, Mod


#Plots failed and successful transitions as one state (4+5)
def laser_brainstate_bootstrap_transitions_together(ppath, recordings, pre, post, edge=0, sf=0, nboots=10000, alpha=0.05, backup='',
                               start_time=0, ma_thr=10, bootstrap_mode=0, fig_file=''):
    """
    Align brain state with laser stimulation and calculate two-sided 1-$alpha confidence intervals using
    bootstrapping.
    :param ppath: base folder
    :param recordings: list of recordings
    :param pre: time before laser
    :param post: time after laser onset
    :param edge: add $edge seconds add beginning and end (that are not shown in the plot) to avoid filtering artifacts
    :param sf: smoothing factor for Gaussian filter; better do not use
    :param nboots: int, how many times the whole data set is resampled for boot-strapping
    :param alpha: plot shows 1-$alpha confidence interval
    :param backup: optional backup folder where recordings are stored
    :param start_time: start time of recordding used for analysis
    :param ma_thr: sleep periods < ma_thr are thrown away
    :param bootstrap_mode: default=0
           bootstrap_mode == 0: Take inter-mouse variance and inter-trial variance (of each mouse) into account.
           That is, bootstrapping re-models the variance expected when re-doing the same
           experimental design (same mouse number and total trial number).
           To account for potentially different number of trials per mouse, resample the data
           during each iteration the following way: Assume that there are n laser trials from m mice;
           randomly select (with replacment) ntrial mice; then select from each mouse randomly one trial.
           bootstrap_mode == 1: Only take inter-trial variance (of each mouse) into account. That is,
           bootstrapping models the variance expected when redoing the experiment with exactly the same mice.
    :param fig_file, if file name is specified, the figure will be saved
    :return: P   - p-values for NREM, REM, Wake
             Mod - by how much the percentage of NREM, REM, Wake is increased compared to baseline
    """

    pre += edge
    post += edge

    rec_paths = dict()
    for rec in recordings:
        if os.path.isdir(os.path.join(ppath, rec)):
            rec_paths[rec] = ppath
        else:
            rec_paths[rec] = backup

    # dict: mouse_id --> laser trials, R W N sequence
    BrainstateDict = {}
    for rec in recordings:
        idf = re.split('_', rec)[0]
        BrainstateDict[idf] = []
    mice = list(BrainstateDict.keys())
    nmice = len(BrainstateDict)

    for rec in recordings:
        ppath = rec_paths[rec]
        SR = get_snr(ppath, rec)
        NBIN = np.round(2.5 * SR)
        dt = NBIN * 1 / SR
        istart_time = int(np.round(start_time / dt))

        M = load_stateidx(ppath, rec)[0]

        seq = get_sequences(np.where(M == 2)[0])
        for s in seq:
            if len(s) * dt <= ma_thr:
                M[s] = 3
        
        seq = get_sequences(np.where(M == 5)[0])
        for s in seq:
            M[s] = 4

        (idxs, idxe) = laser_start_end(load_laser(ppath, rec))
        idf = re.split('_', rec)[0]

        SR = get_snr(ppath, rec)
        NBIN = np.round(2.5 * SR)
        ipre = int(np.round(pre / dt))
        ipost = int(np.round(post / dt))

        idxs = [int(i / NBIN) for i in idxs]
        idxe = [int(i / NBIN) for i in idxe]
        laser_dur = np.mean((np.array(idxe) - np.array(idxs))) * dt

        for (i, j) in zip(idxs, idxe):
            if i >= ipre and j + ipost <= len(M) - 1 and i > istart_time:
                bs = M[i - ipre:i + ipost + 1]
                BrainstateDict[idf].append(bs)

    for mouse in mice:
        BrainstateDict[mouse] = np.array(BrainstateDict[mouse])

    # I assume here that every recording has same dt
    t = np.linspace(-ipre*dt, ipost*dt, ipre+ipost+1)
    Trials = dict()
    for mouse in BrainstateDict:
        Trials[mouse] = np.zeros((BrainstateDict[mouse].shape[0], len(t), 4)) #HERE

    # total number of trials:
    ntrials = 0
    for mouse in BrainstateDict:
        M = np.array(BrainstateDict[mouse])
        for state in range(1, 5): #HERE IMPORTANT CHANGE, went from 1,4 to 1,6
            C = np.zeros(M.shape)
            C[np.where(M == state)] = 100.
            Trials[mouse][:,:,state-1] = C
        ntrials += Trials[mouse].shape[0]

    Prob = np.zeros((nboots, len(t), 4)) #HERE
    if bootstrap_mode == 1:
        for b in range(nboots):
            # average brain state percentage for each mouse during iteration b
            mouse_mean_state = np.zeros((nmice, len(t), 4)) #HERE
            i = 0
            for mouse in mice:
                mmouse = Trials[mouse].shape[0]
                iselect = rand.randint(0, mmouse, (mmouse,))
                for s in [1,2,3,4]: #HERE
                    #bBS[s][offset:offset+mmouse,:] = Trials[mouse][iselect,:,s-1]
                    mouse_mean_state[i,:,s-1] = Trials[mouse][iselect,:,s-1].mean(axis=0)
                i += 1

            for s in [1,2,3,4]: #HERE
                Prob[b,:,s-1] = mouse_mean_state[:,:,s-1].mean(axis=0)
    else:
        mx_iter = np.zeros((ntrials, len(t), 4)) #HERE
        for b in range(nboots):
            # for each iteration select randomly select ntrials mice
            irand_mice = rand.randint(0, nmice, ntrials)

            # average brain state percentage for each mouse during iteration b
            # mouse_mean_state = np.zeros((nmice, len(t), 3))
            i = 0
            # there are ntrials mice
            for j in irand_mice:
                mouse = mice[irand_mice[j]]
                # we have nmouse trials per mouse
                mmouse = Trials[mouse].shape[0]
                # select one random trial from the current mouse
                iselect = rand.randint(0, mmouse)
                for s in [1, 2, 3,4]: #HERE
                    mx_iter[i,:,s-1] = Trials[mouse][iselect,:,s-1]
                i += 1

            # mx_iter is the resampled data set for bootstrap iteration b
            # now we calculate the statistics we're interesting in, which is the mean
            for s in [1, 2, 3,4]: #HERE
                Prob[b,:,s-1] = mx_iter[:,:,s-1].mean(axis=0)


    # simple average for each brainstate across mice (w/o) bootstrapping
    Prob_mean = np.zeros((nmice, len(t), 5)) #HERE
    for s in [1,2,3,4]: #HERE
        i = 0
        for mouse in mice:
            Prob_mean[i,:,s-1] = Trials[mouse][:,:,s-1].mean(axis=0)
            i += 1

    # pdb.set_trace()
    
    usProb = Prob.copy()
    Prob = np.sort(Prob, axis=0)
    Bounds = np.zeros((2, len(t), 4)) #HERE
    a = int((nboots * alpha) / 2.0)
    for s in [1,2,3,4]:#HERE
        Bounds[0,:,s-1] = Prob[a,:,s-1]
        Bounds[1,:,s-1] = Prob[-a,:, s-1]

    # smooth_data
    if sf > 0:
        for s in range(4): #HERE
            Bounds[0, :, s] = smooth_data(Bounds[0, :, s], sf)
            Bounds[1, :, s] = smooth_data(Bounds[1, :, s], sf)
        for i in range(nmice):
            for s in range(4): #HERE
                Prob_mean[i, :, s] = smooth_data(Prob_mean[i,:,s], sf)

    # plot figure
    colors = np.array([[0, 1, 1], [0.5, 0, 1], [0.6, 0.6, 0.6],[0, 0, .5]]) #HERE
    br_states = {1:'REM', 2:'Wake', 3:'NREM', 4:'Trans'} #HERE
    #colors = np.array([[55,255,255], [153,255,153],[153,153,153]])/255.
    it = np.where((t>=-pre+edge) & (t<=post-edge))[0]
    plt.ion()
    plt.figure()
    ax = plt.axes([0.15, 0.15, 0.6, 0.7])
    # ax.add_patch(patches.Rectangle((0, 0), laser_dur, 100, facecolor=[0.6, 0.6, 1], edgecolor=[0.6, 0.6, 1]))

    for state in [4,3,2,1]: #HERE
        if state == 4:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=4, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=4)
        elif state == 3:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=1, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=1)
        elif state == 2:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=2, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=2)
        elif state == 1:
            ax.fill_between(t[it], Bounds[0,it,state-1], Bounds[1,it,state-1], color=colors[state-1,:], alpha=0.8, zorder=3, edgecolor=None)
            ax.plot(t[it], Prob_mean[:, it, state-1].mean(axis=0), color=colors[state-1,:], label=br_states[state], zorder=3)
        
        
    ax.add_patch(patches.Rectangle((0, 0), laser_dur, 100, facecolor=[0.6, 0.6, 1], edgecolor=[0.6, 0.6, 1]))
    plt.xlim([-pre+edge, post-edge])
    plt.ylim([0,100])
    plt.xlabel('Time (s)')
    plt.ylabel('Brain state (%)')
    plt.legend(bbox_to_anchor = (1.0, 0.7, 1., .102), loc = 4, mode = 'expand', ncol = 1, frameon = False) #HERE
    box_off(ax)
    plt.draw()

    # statistics
    ibase = np.where((t>=-laser_dur) & (t<0))[0]
    ilsr  = np.where((t>=0) & (t<laser_dur))[0]
    P   = np.zeros((4,)) #HERE
    Mod = np.zeros((4,)) #HERE
    for istate in [1,2,3,4]: #HERE
        basel = usProb[:,ibase,istate-1].mean(axis=1)
        laser = usProb[:,ilsr, istate-1].mean(axis=1)
        d = laser - basel
        if np.mean(d) >= 0:
            # now we want all values be larger than 0
            p = len(np.where(d>0)[0]) / (1.0*nboots)
            sig = 1 - p
            if sig == 0:
                sig = 1.0/nboots
            Mod[istate-1] = (np.mean(laser) / np.mean(basel) - 1) * 100
        else:
            p = len(np.where(d<0)[0]) / (1.0*nboots)
            sig = 1 - p
            if sig == 0:
                sig = 1.0/nboots
            Mod[istate-1] = -(1 - np.mean(laser) / np.mean(basel)) * 100
        P[istate-1] = sig

    labels = {1:'REM', 2:'Wake', 3:'NREM', 4:'Trans'} #HERE
    for s in [1,2,3,4]: #HERE
        print('%s is changed by %f perc.; P = %f, bootstrap' % (labels[s], Mod[s-1], P[s-1]))
    print("n = %d mice" % len(mice))

    if len(fig_file) > 0:
        plt.savefig(fig_file, bbox_inches="tight")

    return P, Mod




    
#simple sleep spectrum function to compare spectrums in recordings without laser
#For more in depth script that can plot multiple conditions (not just 2) use spectrum_conditions.py / spectrum_conditions_transitions.py
def sleep_spectrum_dreadds(ppath, salineRecordings, cnoRecordings, istate, pmode=0, twin=3, ma_thr=20.0, f_max=15, pplot=True, sig_type='EEG', mu=[10, 100],
                   tstart=0, tend=-1, sthres=np.inf, peeg2=False, pnorm=False, single_mode=False, conv=1.0, fig_file='', laser_color='blue', plotSignificance = False):
    """
    calculate power spectrum for brain state i state for the given recordings 
    @Param:
    ppath    -    folder containing all recordings
    recordings -  single recording (string) or list of recordings
    @Optional:
    istate   -    state for which to calculate power spectrum
    twin     -    time window (in seconds) for power spectrum calculation
                  the longer the higher frequency resolution, but the more noisy
    ma_thr   -    short wake periods <= $ma_thr are considered as sleep
    f_max    -    maximal frequency, if f_max==-1: f_max is maximally possible frequency
    pplot    -    if True, plot figure showing result
    pmode    -    mode: 
                  pmode == 1, compare state during laser with baseline outside laser interval
                  pmode == 0, just plot power spectrum for state istate and don't care about laser
    tstart   -    use EEG starting from time point tstart [seconds]
    sig_type -    string, if 'EMG' calculate EMG amplitude (from the EMG spectrum). E.g.,
                  sleepy.sleep_spectrum(ppath, E, istate=2, f_max=30, sig_type='EMG')
    mu       -    tuple, lower and upper range for EMG frequencies used for amplitude calculation
    tend     -    use data up to tend [seconds], if tend == -1, use data till end
    sthres   -    maximum length of bout duration of state $istate used for calculation. If bout duration > $sthres, only
                  use the bout up to $sthres seconds after bout onset.
    peeg2    -    if True, use EEG2 channel for spectrum analysis
    pnorm    -    if True, normalize powerspectrum by dividing each frequency through each average power
                  over the whole EEG recording
    single_mode - if True, plot each single mouse
    fig_file -    if specified save to given file

    errorbars: If it's multiple mice make errorbars over mice; if it's multiple
    recordings of ONE mouse, show errorbars across recordings; 
    if just one recording show now errorbars
                  
    @Return:
    Pow     -    Dict[No loaser = 0|Laser = 1][array], where array: mice x frequencies, if more than one mouse;
                 otherwise, array: recordings x frequencies
    F       -    Frequencies
    """
    both = [salineRecordings, cnoRecordings]
   
    
    for recordingList in both:
      
        recordings = recordingList
        if type(recordings) != list:
            recordings = [recordings]
    
        Mice = {}
        for rec in recordings:
            idf = re.split('_', rec)[0]
            if not(idf in Mice):
                Mice[idf] = Mouse(idf, rec, 'E')
            else:
                Mice[idf].add(rec)
    
        mouse_order = []
        for rec in recordings:
            idf = re.split('_', rec)[0]
            if not idf in mouse_order:
                mouse_order.append(idf)
    
        # Spectra: Dict[mouse_id][laser_on|laser_off][list of powerspectrum_arrays]
        Spectra = {}
        Ids = list(Mice.keys())
        for i in Ids:
            Spectra[i] = {0:[], 1:[]}
            Spectra[i] = {0:[], 1:[]}
    
        for idf in mouse_order:
            for rec in Mice[idf].recordings:
                # load EEG
                if sig_type =='EEG':
                    if not peeg2:
                        EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG.mat'))['EEG']).astype('float')*conv
                    else:
                        EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG2.mat'))['EEG2']).astype('float')*conv
                elif sig_type == 'EMG':
                    if not peeg2:
                        EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG.mat'))['EMG']).astype('float')*conv
                    else:
                        EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG2.mat'))['EMG2']).astype('float')*conv
                else:
                    pass
    
                # load brain state
                M,S = load_stateidx(ppath, rec)
                sr = get_snr(ppath, rec)
                # number of time bins for each time bin in spectrogram
                nbin = int(np.round(sr) * 2.5)
                # duration of time bin in spectrogram / brainstate
                dt = nbin * 1/sr
                nwin = np.round(twin*sr)
    
                istart = int(np.round(tstart/dt))
                if tend==-1:
                    iend = M.shape[0]
                else:
                    iend = int(np.round(tend/dt))
                istart_eeg = istart*nbin
                iend_eeg   = (iend-1)*nbin+1
    
                #Combine NRt states
                M[np.where(M==5)]=4
                
                # flatten out microarousals
                seq = get_sequences(np.where(M==2)[0])
                for s in seq:
                    if len(s)*dt <= ma_thr:
                        M[s] = 3
    
                # get all sequences of state $istate
                M = M[istart:iend]
                seq = get_sequences(np.where(M==istate)[0])
    
                EEG = EEG[istart_eeg:iend_eeg]
    
                if pnorm:
                    pow_norm = power_spectrum(EEG, nwin, 1 / sr)[0]
    
                if pmode == 1:
                    laser = load_laser(ppath, rec)[istart_eeg:iend_eeg]
                    (idxs, idxe) = laser_start_end(laser, SR=sr)
                    # downsample EEG time to spectrogram time    
                    idxs = [int(i/nbin) for i in idxs]
                    idxe = [int(i/nbin) for i in idxe]
                    
                    laser_idx = []
                    for (i,j) in zip(idxs, idxe):
                        laser_idx += list(range(i,j+1))
                    laser_idx = np.array(laser_idx)
                    
                if pmode == 1:
                    print('Your princess is in another castle. (pmode = 1 is for laser and this is a dreadds function now)')
                    return
                            
                # don't care about laser
                if pmode == 0:
                    for s in seq:
                        if len(s)*nbin >= nwin:
                            drn = (s[-1]-s[0])*dt
                            if drn > sthres:
                                b = (s[0] + int(np.round(sthres/dt)))*nbin
                            else:
                                b = int((s[-1]+1)*nbin)
                            sup = list(range(int(s[0]*nbin), b))
                            if sup[-1]>len(EEG):
                                sup = list(range(int(s[0]*nbin), len(EEG)))
    
                            if len(s)*nbin >= nwin:
                                Pow, F = power_spectrum(EEG[sup], nwin, 1/sr) #Gets power spectrum and corresponding frequencies
                                if pnorm:
                                    Pow = np.divide(Pow, pow_norm)
                                Spectra[idf][0].append(Pow)
                
                ##############################   
                if seq[0].size == 0: #Skip on to the next recording if there are no istate bouts in this rec. If there are no rem episodes for instance.
                    # pdb.set_trace()
                    continue

                
                Pow = {0:[], 1:[]}
                if len(Ids)==1:
                    # only one mouse
                    Pow[0] = np.array(Spectra[Ids[0]][0])
                    Pow[1] = np.array(Spectra[Ids[0]][1])
                else:
                    # several mice
                    Pow[0] = np.zeros((len(Ids),len(F)))
                    Pow[1] = np.zeros((len(Ids),len(F)))
                    i = 0
                    for m in Ids:
                        Pow[0][i,:] = np.array(Spectra[m][0]).mean(axis=0)
                        if pmode == 1:
                            Pow[1][i,:] = np.array(Spectra[m][1]).mean(axis=0)
                        i += 1
                        
                        
#        pdb.set_trace()
#        if plotSignificance:
#            x=1
#            #call ttest function
#            
#            #plot bars for delta, theta, sigma power at top of graph axes
#            
#            #plot significance asterisk above bars
        
        
        #now plot stuff
        if f_max > -1:
            ifreq = np.where(F<=f_max)[0]
            F = F[ifreq]
            Pow[0] = Pow[0][:,ifreq]
            if pmode==1:
                Pow[1] = Pow[1][:,ifreq]
        else:
            f_max = F[-1]
        
        if pplot:
            plt.ion()
            if recordings == salineRecordings:
                plt.figure()
            if sig_type == 'EEG':
                ax = plt.axes([0.2, 0.15, 0.6, 0.7])
                n = Pow[0].shape[0]
                clrs = sns.color_palette("husl", len(mouse_order))
                if pmode==1:
                    print('Still pmode1')
                    return
       
                #IF PMODE = 0 
                if not single_mode:                    
                    asdf = np.nanmean(Pow[0], axis=0)-np.nanstd(Pow[0], axis=0)/np.sqrt(n)-np.nanstd(Pow[0], axis=0)/np.sqrt(n)
                    bsdf = np.nanmean(Pow[0], axis=0)+np.nanstd(Pow[0], axis=0)/np.sqrt(n)+np.nanstd(Pow[0], axis=0)/np.sqrt(n)
                    if recordings == salineRecordings:
                        # pdb.set_trace()
                        plt.fill_between(F, asdf, bsdf, alpha=0.5, color='gray')
                        # plt.plot(F, Pow[0].mean(axis=0), color='gray', lw=2, alpha=0.5, label='SALINE')
                        plt.plot(F, np.nanmean(Pow[0],axis=0), color='gray', lw=2, alpha=0.5, label='SALINE')
                    if recordings == cnoRecordings:
                        plt.fill_between(F, asdf, bsdf, alpha=0.5, color='blue')
                        plt.plot(F, np.nanmean(Pow[0], axis=0), color='blue', lw=2, alpha=0.5, label='CNO')
                        
                else:
                    for i in range(len(mouse_order)):
                        plt.plot(F, Pow[0][i, :], label=mouse_order[i], color=clrs[i])
                    plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=len(mouse_order), frameon=False)
    
                if pmode==1 and not single_mode:
                    plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
    
                box_off(ax)
                plt.xlim([0, f_max])
                plt.xlabel('Freq. (Hz)')
                plt.ylabel('Power ($\mathrm{\mu V^2}$)')
                
                plt.show()
    
            #IF PLOTTING EMG
            else:
                # plot EMG amplitude
                nmice = len(mouse_order)
                clrs = sns.color_palette("husl", nmice)
                # plot EMG
                Ampl = {}
                # range of frequencies
                mfreq = np.where((F >= mu[0]) & (F <= mu[1]))[0]
                df = F[1] - F[0]
                if pmode==1:
                    for i in [0, 1]:
                        Ampl[i] = np.sqrt(Pow[i][:,mfreq].sum(axis=1)*df)
                else:
                    Ampl[0] = np.sqrt(Pow[0][:,mfreq].sum(axis=1)*df)
    
                if pmode==1:
                    ax = plt.axes([0.2, 0.15, 0.4, 0.7])
                    ax.bar([0], Ampl[0].mean(), color='gray', label='w/o laser')
                    ax.bar([1], Ampl[1].mean(), color='blue', label='laser')
                    plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
    
                    for i in range(nmice):
                        plt.plot([0,1], [Ampl[0][i], Ampl[1][i]], color=clrs[i], label=mouse_order[i])
    
                    box_off(ax)
                    plt.ylabel('EMG Ampl. ($\mathrm{\mu V}$)')
                    ax.set_xticks([0,1])
                    ax.set_xticklabels(['', ''])
    
                    # some basic stats
                    [tstats, p] = stats.ttest_rel(Ampl[0], Ampl[1])
                    print("Stats for EMG amplitude: t-statistics: %.3f, p-value: %.3f" % (tstats, p))
    
        if len(fig_file) > 0:
            save_figure(fig_file)
    
'''
Plot spectrums for time bins by hours (on the same graph) for the CNO recordings, with saline as 'base case'
'''  
def sleep_spectrum_dreadds_timebins(ppath, salineRecordings, cnoRecordings, istate, pmode=0, twin=3, ma_thr=20.0, f_max=30, pplot=True, sig_type='EEG', mu=[10, 100],
                   tstart=2700, tend=-1, sthres=np.inf, peeg2=False, pnorm=False, single_mode=False, conv=1.0, fig_file='', laser_color='blue'):

    both = [salineRecordings, cnoRecordings] 
    
    for recordingList in both:
        
        if recordingList == both[1]: # if it's the cno recordings do it by timebins
            
            tstartList = [0,3599,7199,10799,14399]
            tendList = [3600,7200,10800,14400,18000]
            colornum = 0
            colors = ['gold','orange', 'darkorange', 'red', 'brown']
            binNum = ['Hour 1','Hour 2','Hour 3','Hour 4','Hour 5']
            for (tstart,tend) in zip(tstartList,tendList):
            
                recordings = recordingList
                if type(recordings) != list:
                    recordings = [recordings]
            
                Mice = {}
                for rec in recordings:
                    idf = re.split('_', rec)[0]
                    if not(idf in Mice):
                        Mice[idf] = Mouse(idf, rec, 'E')
                    else:
                        Mice[idf].add(rec)
            
                mouse_order = []
                for rec in recordings:
                    idf = re.split('_', rec)[0]
                    if not idf in mouse_order:
                        mouse_order.append(idf)
            
                # Spectra: Dict[mouse_id][laser_on|laser_off][list of powerspectrum_arrays]
                Spectra = {}
                Ids = list(Mice.keys())
                for i in Ids:
                    Spectra[i] = {0:[], 1:[]}
                    Spectra[i] = {0:[], 1:[]}
            
                for idf in mouse_order:
                    for rec in Mice[idf].recordings:
                        # load EEG
                        if sig_type =='EEG':
                            if not peeg2:
                                EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG.mat'))['EEG']).astype('float')*conv
                            else:
                                EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG2.mat'))['EEG2']).astype('float')*conv
                        elif sig_type == 'EMG':
                            if not peeg2:
                                EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG.mat'))['EMG']).astype('float')*conv
                            else:
                                EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG2.mat'))['EMG2']).astype('float')*conv
                        else:
                            pass
            
                        # load brain state
                        M,S = load_stateidx(ppath, rec)
                        sr = get_snr(ppath, rec)
                        # number of time bins for each time bin in spectrogram
                        nbin = int(np.round(sr) * 2.5)
                        # duration of time bin in spectrogram / brainstate
                        dt = nbin * 1/sr
                        nwin = np.round(twin*sr)
            
                        istart = int(np.round(tstart/dt))
                        if tend==-1:
                            iend = M.shape[0]
                        else:
                            iend = int(np.round(tend/dt))
                        istart_eeg = istart*nbin
                        iend_eeg   = (iend-1)*nbin+1
            
                        M[np.where(M==5)]=2
                        # flatten out microarousals
                        seq = get_sequences(np.where(M==2)[0])
                        for s in seq:
                            if len(s)*dt <= ma_thr:
                                M[s] = 3
            
                        # get all sequences of state $istate
                        M = M[istart:iend]
                        seq = get_sequences(np.where(M==istate)[0])
            
                        EEG = EEG[istart_eeg:iend_eeg]
            
                        if pnorm:
                            pow_norm = power_spectrum(EEG, nwin, 1 / sr)[0]
            
                        if pmode == 1:
                            laser = load_laser(ppath, rec)[istart_eeg:iend_eeg]
                            (idxs, idxe) = laser_start_end(laser, SR=sr)
                            # downsample EEG time to spectrogram time    
                            idxs = [int(i/nbin) for i in idxs]
                            idxe = [int(i/nbin) for i in idxe]
                            
                            laser_idx = []
                            for (i,j) in zip(idxs, idxe):
                                laser_idx += list(range(i,j+1))
                            laser_idx = np.array(laser_idx)
                            
                        if pmode == 1:
                            # first analyze frequencies not overlapping with laser
                            seq_nolsr = []
                            for s in seq:
                                s = np.setdiff1d(s, laser_idx)
                                if len(s) > 0:
                                    q = get_sequences(s)
                                    seq_nolsr += q
            
                            for s in seq_nolsr:
                                #s = np.setdiff1d(s, laser_idx)
                                if len(s)*nbin >= nwin:
                                    drn = (s[-1]-s[0])*dt
                                    if drn > sthres:
                                        # b is the end of segment used for power spectrum calculation;
                                        # that is, the last index (in raw EEG) of the segment
                                        b = (s[0] + int(np.round(sthres/dt)))*nbin
                                    else:
                                        b = int((s[-1]+1)*nbin)
            
                                    sup = list(range(int(s[0]*nbin), b))
            
                                    if sup[-1]>len(EEG):
                                        sup = list(range(int(s[0]*nbin), len(EEG)))
                                    if len(s)*nbin >= nwin:
                                        Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                        if pnorm:
                                            Pow = np.divide(Pow, pow_norm)
                                        Spectra[idf][0].append(Pow)
                                    
                            # now analyze sequences overlapping with laser
                            seq_lsr = []
                            for s in seq:
                                s = np.intersect1d(s, laser_idx)
                                if len(s) > 0:
                                    q = get_sequences(s)
                                    seq_lsr += q
            
                            for s in seq_lsr:
                                s = np.intersect1d(s, laser_idx)
                                
                                if len(s)*nbin >= nwin:
                                    # calculate power spectrum
                                    # upsample indices
                                    # brain state time 0     1         2
                                    # EEG time         0-999 1000-1999 2000-2999
                                    drn = (s[-1]-s[0])*dt
                                    if drn > sthres:
                                        b = (s[0] + int(np.round(sthres/dt)))*nbin
                                        #pdb.set_trace()
                                    else:
                                        b = int((s[-1]+1)*nbin)
            
                                    sup = list(range(int(s[0]*nbin), b))
            
                                    if sup[-1]>len(EEG):
                                        sup = list(range(int(s[0]*nbin), len(EEG)))
                                    if len(s)*nbin >= nwin:
                                        Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                        if pnorm:
                                            Pow = np.divide(Pow, pow_norm)
                                        Spectra[idf][1].append(Pow)
                                    
                        # don't care about laser
                        if pmode == 0:
                            for s in seq:
                                if len(s)*nbin >= nwin:
                                    drn = (s[-1]-s[0])*dt
                                    if drn > sthres:
                                        b = (s[0] + int(np.round(sthres/dt)))*nbin
                                    else:
                                        b = int((s[-1]+1)*nbin)
                                    sup = list(range(int(s[0]*nbin), b))
                                    if sup[-1]>len(EEG):
                                        sup = list(range(int(s[0]*nbin), len(EEG)))
            
                                    if len(s)*nbin >= nwin:
                                        Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                        if pnorm:
                                            Pow = np.divide(Pow, pow_norm)
                                        Spectra[idf][0].append(Pow)
                            
                        Pow = {0:[], 1:[]}
                        if len(Ids)==1:
                            # only one mouse
                            Pow[0] = np.array(Spectra[Ids[0]][0])
                            Pow[1] = np.array(Spectra[Ids[0]][1])
                        else:
                            # several mice
                            Pow[0] = np.zeros((len(Ids),len(F)))
                            Pow[1] = np.zeros((len(Ids),len(F)))
                            i = 0
                            for m in Ids:
                                Pow[0][i,:] = np.array(Spectra[m][0]).mean(axis=0)
                                if pmode == 1:
                                    Pow[1][i,:] = np.array(Spectra[m][1]).mean(axis=0)
                                i += 1
            
                if f_max > -1:
                    ifreq = np.where(F<=f_max)[0]
                    F = F[ifreq]
                    Pow[0] = Pow[0][:,ifreq]
                    if pmode==1:
                        Pow[1] = Pow[1][:,ifreq]
                else:
                    f_max = F[-1]
                
                if pplot:
                    plt.ion()
                    if recordings == salineRecordings:
                        plt.figure()
                    if sig_type == 'EEG':
                        ax = plt.axes([0.2, 0.15, 0.6, 0.7])
                        n = Pow[0].shape[0]
                        clrs = sns.color_palette("husl", len(mouse_order))
                        if pmode==1:
                            if not single_mode:
                                a = Pow[1].mean(axis=0) - Pow[1].std(axis=0) / np.sqrt(n)
                                b = Pow[1].mean(axis=0) + Pow[1].std(axis=0) / np.sqrt(n)
                                plt.fill_between(F, a, b, alpha=0.5, color=laser_color)
                                plt.plot(F, Pow[1].mean(axis=0), color=laser_color, lw=2, label='With laser')
                            else:
                                for i in range(len(mouse_order)):                                                        
                                    plt.plot(F, Pow[1][i,:], '--', color=clrs[i])                            
                                #plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=len(mouse_order),
                                #           frameon=False)
            
            
                        #IF PMODE = 0 
                        if not single_mode:
                            a = Pow[0].mean(axis=0)-Pow[0].std(axis=0)/np.sqrt(n)
                            b = Pow[0].mean(axis=0)+Pow[0].std(axis=0)/np.sqrt(n)
                            if recordings == salineRecordings:
                                plt.fill_between(F, a, b, alpha=0.5, color='gray')
                                plt.plot(F, Pow[0].mean(axis=0), color='gray', lw=2, alpha=0.5, label='SALINE')
                            if recordings == cnoRecordings:
                                plt.fill_between(F, a, b, alpha=0, color=colors[colornum])
                                plt.plot(F, Pow[0].mean(axis=0), color=colors[colornum], lw=2, alpha=0.5,label=binNum[colornum]) #label=binNum[colornum]
                                plt.legend(loc='best')
                                
                        else:
                            for i in range(len(mouse_order)):
                                plt.plot(F, Pow[0][i, :], label=mouse_order[i], color=clrs[i])
                            plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=len(mouse_order), frameon=False)
            
                        if pmode==1 and not single_mode:
                            plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
            
                        box_off(ax)
                        plt.xlim([0, f_max])
                        plt.xlabel('Freq. (Hz)')
                        plt.ylabel('Power ($\mathrm{\mu V^2}$)')
                        
                        plt.show()
            
                    #IF PLOTTING EMG
                    else:
                        # plot EMG amplitude
                        nmice = len(mouse_order)
                        clrs = sns.color_palette("husl", nmice)
                        # plot EMG
                        Ampl = {}
                        # range of frequencies
                        mfreq = np.where((F >= mu[0]) & (F <= mu[1]))[0]
                        df = F[1] - F[0]
                        if pmode==1:
                            for i in [0, 1]:
                                Ampl[i] = np.sqrt(Pow[i][:,mfreq].sum(axis=1)*df)
                        else:
                            Ampl[0] = np.sqrt(Pow[0][:,mfreq].sum(axis=1)*df)
            
                        if pmode==1:
                            ax = plt.axes([0.2, 0.15, 0.4, 0.7])
                            ax.bar([0], Ampl[0].mean(), color='gray', label='w/o laser')
                            ax.bar([1], Ampl[1].mean(), color='blue', label='laser')
                            plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
            
                            for i in range(nmice):
                                plt.plot([0,1], [Ampl[0][i], Ampl[1][i]], color=clrs[i], label=mouse_order[i])
            
                            box_off(ax)
                            plt.ylabel('EMG Ampl. ($\mathrm{\mu V}$)')
                            ax.set_xticks([0,1])
                            ax.set_xticklabels(['', ''])
            
                            # some basic stats
                            [tstats, p] = stats.ttest_rel(Ampl[0], Ampl[1])
                            print("Stats for EMG amplitude: t-statistics: %.3f, p-value: %.3f" % (tstats, p))
            
                if len(fig_file) > 0:
                    save_figure(fig_file) 
                    
                colornum+=1     
                    
                    
                    
        else:

            recordings = recordingList
            if type(recordings) != list:
                recordings = [recordings]
        
            Mice = {}
            for rec in recordings:
                idf = re.split('_', rec)[0]
                if not(idf in Mice):
                    Mice[idf] = Mouse(idf, rec, 'E')
                else:
                    Mice[idf].add(rec)
        
            mouse_order = []
            for rec in recordings:
                idf = re.split('_', rec)[0]
                if not idf in mouse_order:
                    mouse_order.append(idf)
        
            # Spectra: Dict[mouse_id][laser_on|laser_off][list of powerspectrum_arrays]
            Spectra = {}
            Ids = list(Mice.keys())
            for i in Ids:
                Spectra[i] = {0:[], 1:[]}
                Spectra[i] = {0:[], 1:[]}
        
            for idf in mouse_order:
                for rec in Mice[idf].recordings:
                    # load EEG
                    if sig_type =='EEG':
                        if not peeg2:
                            EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG.mat'))['EEG']).astype('float')*conv
                        else:
                            EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EEG2.mat'))['EEG2']).astype('float')*conv
                    elif sig_type == 'EMG':
                        if not peeg2:
                            EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG.mat'))['EMG']).astype('float')*conv
                        else:
                            EEG = np.squeeze(so.loadmat(os.path.join(ppath, rec, 'EMG2.mat'))['EMG2']).astype('float')*conv
                    else:
                        pass
        
                    # load brain state
                    M,S = load_stateidx(ppath, rec)
                    sr = get_snr(ppath, rec)
                    # number of time bins for each time bin in spectrogram
                    nbin = int(np.round(sr) * 2.5)
                    # duration of time bin in spectrogram / brainstate
                    dt = nbin * 1/sr
                    nwin = np.round(twin*sr)
        
                    istart = int(np.round(tstart/dt))
                    if tend==-1:
                        iend = M.shape[0]
                    else:
                        iend = int(np.round(tend/dt))
                    istart_eeg = istart*nbin
                    iend_eeg   = (iend-1)*nbin+1
        
                    M[np.where(M==5)]=2
                    # flatten out microarousals
                    seq = get_sequences(np.where(M==2)[0])
                    for s in seq:
                        if len(s)*dt <= ma_thr:
                            M[s] = 3
        
                    # get all sequences of state $istate
                    M = M[istart:iend]
                    seq = get_sequences(np.where(M==istate)[0])
        
                    EEG = EEG[istart_eeg:iend_eeg]
        
                    if pnorm:
                        pow_norm = power_spectrum(EEG, nwin, 1 / sr)[0]
        
                    if pmode == 1:
                        laser = load_laser(ppath, rec)[istart_eeg:iend_eeg]
                        (idxs, idxe) = laser_start_end(laser, SR=sr)
                        # downsample EEG time to spectrogram time    
                        idxs = [int(i/nbin) for i in idxs]
                        idxe = [int(i/nbin) for i in idxe]
                        
                        laser_idx = []
                        for (i,j) in zip(idxs, idxe):
                            laser_idx += list(range(i,j+1))
                        laser_idx = np.array(laser_idx)
                        
                    if pmode == 1:
                        # first analyze frequencies not overlapping with laser
                        seq_nolsr = []
                        for s in seq:
                            s = np.setdiff1d(s, laser_idx)
                            if len(s) > 0:
                                q = get_sequences(s)
                                seq_nolsr += q
        
                        for s in seq_nolsr:
                            #s = np.setdiff1d(s, laser_idx)
                            if len(s)*nbin >= nwin:
                                drn = (s[-1]-s[0])*dt
                                if drn > sthres:
                                    # b is the end of segment used for power spectrum calculation;
                                    # that is, the last index (in raw EEG) of the segment
                                    b = (s[0] + int(np.round(sthres/dt)))*nbin
                                else:
                                    b = int((s[-1]+1)*nbin)
        
                                sup = list(range(int(s[0]*nbin), b))
        
                                if sup[-1]>len(EEG):
                                    sup = list(range(int(s[0]*nbin), len(EEG)))
                                if len(s)*nbin >= nwin:
                                    Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                    if pnorm:
                                        Pow = np.divide(Pow, pow_norm)
                                    Spectra[idf][0].append(Pow)
                                
                        # now analyze sequences overlapping with laser
                        seq_lsr = []
                        for s in seq:
                            s = np.intersect1d(s, laser_idx)
                            if len(s) > 0:
                                q = get_sequences(s)
                                seq_lsr += q
        
                        for s in seq_lsr:
                            s = np.intersect1d(s, laser_idx)
                            
                            if len(s)*nbin >= nwin:
                                # calculate power spectrum
                                # upsample indices
                                # brain state time 0     1         2
                                # EEG time         0-999 1000-1999 2000-2999
                                drn = (s[-1]-s[0])*dt
                                if drn > sthres:
                                    b = (s[0] + int(np.round(sthres/dt)))*nbin
                                    #pdb.set_trace()
                                else:
                                    b = int((s[-1]+1)*nbin)
        
                                sup = list(range(int(s[0]*nbin), b))
        
                                if sup[-1]>len(EEG):
                                    sup = list(range(int(s[0]*nbin), len(EEG)))
                                if len(s)*nbin >= nwin:
                                    Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                    if pnorm:
                                        Pow = np.divide(Pow, pow_norm)
                                    Spectra[idf][1].append(Pow)
                                
                    # don't care about laser
                    if pmode == 0:
                        for s in seq:
                            if len(s)*nbin >= nwin:
                                drn = (s[-1]-s[0])*dt
                                if drn > sthres:
                                    b = (s[0] + int(np.round(sthres/dt)))*nbin
                                else:
                                    b = int((s[-1]+1)*nbin)
                                sup = list(range(int(s[0]*nbin), b))
                                if sup[-1]>len(EEG):
                                    sup = list(range(int(s[0]*nbin), len(EEG)))
        
                                if len(s)*nbin >= nwin:
                                    Pow, F = power_spectrum(EEG[sup], nwin, 1/sr)
                                    if pnorm:
                                        Pow = np.divide(Pow, pow_norm)
                                    Spectra[idf][0].append(Pow)
                        
                    Pow = {0:[], 1:[]}
                    if len(Ids)==1:
                        # only one mouse
                        Pow[0] = np.array(Spectra[Ids[0]][0])
                        Pow[1] = np.array(Spectra[Ids[0]][1])
                    else:
                        # several mice
                        Pow[0] = np.zeros((len(Ids),len(F)))
                        Pow[1] = np.zeros((len(Ids),len(F)))
                        i = 0
                        for m in Ids:
                            Pow[0][i,:] = np.array(Spectra[m][0]).mean(axis=0)
                            if pmode == 1:
                                Pow[1][i,:] = np.array(Spectra[m][1]).mean(axis=0)
                            i += 1
        
            if f_max > -1:
                ifreq = np.where(F<=f_max)[0]
                F = F[ifreq]
                Pow[0] = Pow[0][:,ifreq]
                if pmode==1:
                    Pow[1] = Pow[1][:,ifreq]
            else:
                f_max = F[-1]
            
            if pplot:
                plt.ion()
                if recordings == salineRecordings:
                    plt.figure()
                if sig_type == 'EEG':
                    ax = plt.axes([0.2, 0.15, 0.6, 0.7])
                    n = Pow[0].shape[0]
                    clrs = sns.color_palette("husl", len(mouse_order))
                    if pmode==1:
                        if not single_mode:
                            a = Pow[1].mean(axis=0) - Pow[1].std(axis=0) / np.sqrt(n)
                            b = Pow[1].mean(axis=0) + Pow[1].std(axis=0) / np.sqrt(n)
                            plt.fill_between(F, a, b, alpha=0.5, color=laser_color)
                            plt.plot(F, Pow[1].mean(axis=0), color=laser_color, lw=2, label='With laser')
                        else:
                            for i in range(len(mouse_order)):                                                        
                                plt.plot(F, Pow[1][i,:], '--', color=clrs[i])                            
                            #plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=len(mouse_order),
                            #           frameon=False)
        
        
                    #IF PMODE = 0 
                    if not single_mode:
                        a = Pow[0].mean(axis=0)-Pow[0].std(axis=0)/np.sqrt(n)
                        b = Pow[0].mean(axis=0)+Pow[0].std(axis=0)/np.sqrt(n)
                        if recordings == salineRecordings:
                            plt.fill_between(F, a, b, alpha=0, color='gray')
                            plt.plot(F, Pow[0].mean(axis=0), color='black', lw=2, alpha=0.5, label='SALINE')
                            plt.legend(loc='best')
                        if recordings == cnoRecordings:
                            plt.fill_between(F, a, b, alpha=0.5, color='blue')
                            plt.plot(F, Pow[0].mean(axis=0), color='blue', lw=2, alpha=0.5, label='CNO')
                            
                    else:
                        for i in range(len(mouse_order)):
                            plt.plot(F, Pow[0][i, :], label=mouse_order[i], color=clrs[i])
                        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', ncol=len(mouse_order), frameon=False)
        
                    if pmode==1 and not single_mode:
                        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
        
                    box_off(ax)
                    plt.xlim([0, f_max])
                    plt.xlabel('Freq. (Hz)')
                    plt.ylabel('Power ($\mathrm{\mu V^2}$)')
                    
                    plt.show()
        
                #IF PLOTTING EMG
                else:
                    # plot EMG amplitude
                    nmice = len(mouse_order)
                    clrs = sns.color_palette("husl", nmice)
                    # plot EMG
                    Ampl = {}
                    # range of frequencies
                    mfreq = np.where((F >= mu[0]) & (F <= mu[1]))[0]
                    df = F[1] - F[0]
                    if pmode==1:
                        for i in [0, 1]:
                            Ampl[i] = np.sqrt(Pow[i][:,mfreq].sum(axis=1)*df)
                    else:
                        Ampl[0] = np.sqrt(Pow[0][:,mfreq].sum(axis=1)*df)
        
                    if pmode==1:
                        ax = plt.axes([0.2, 0.15, 0.4, 0.7])
                        ax.bar([0], Ampl[0].mean(), color='gray', label='w/o laser')
                        ax.bar([1], Ampl[1].mean(), color='blue', label='laser')
                        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, mode='expand', frameon=False)
        
                        for i in range(nmice):
                            plt.plot([0,1], [Ampl[0][i], Ampl[1][i]], color=clrs[i], label=mouse_order[i])
        
                        box_off(ax)
                        plt.ylabel('EMG Ampl. ($\mathrm{\mu V}$)')
                        ax.set_xticks([0,1])
                        ax.set_xticklabels(['', ''])
        
                        # some basic stats
                        [tstats, p] = stats.ttest_rel(Ampl[0], Ampl[1])
                        print("Stats for EMG amplitude: t-statistics: %.3f, p-value: %.3f" % (tstats, p))
        
            if len(fig_file) > 0:
                save_figure(fig_file)    
    
    
    


#Compute the significance of each time point from sleepy.sleep_timecourse. Needs post-hoc correction as well (divide computed p-value by # of time bins)
def timecourse_ttest(ppath, recordingFile):

    condition = 'dreaddCNO'
#    condition = 'mcherryCNO'
    
    statsList = ['perc']
    brainStateList = ['REM', 'Wake', 'NREM']    
    timebinPValuesDF = pd.DataFrame()
    labels = ['Stat', 'Brain State', 'time', 'pvalue']
    
    for eachStat in statsList:
        stat = eachStat
        TimeMxCtr, TimeMxExp, df = sleepy.sleep_timecourse(ppath, trace_file = recordingFile, tbin = 1800, n = 10, stats = stat)
        
        #start reorganizing the dataframe by mouse
        
#        pdb.set_trace()
        
        salineDF = df.loc[df['dose'] == '0']
        time0Saline = salineDF.loc[salineDF['time'] == 't0']
        time1Saline =  salineDF.loc[salineDF['time'] == 't1']
        time2Saline =  salineDF.loc[salineDF['time'] == 't2']
        time3Saline =  salineDF.loc[salineDF['time'] == 't3']
        time4Saline =  salineDF.loc[salineDF['time'] == 't4']
        time5Saline =  salineDF.loc[salineDF['time'] == 't5']
        time6Saline =  salineDF.loc[salineDF['time'] == 't6']
        time7Saline =  salineDF.loc[salineDF['time'] == 't7']
        time8Saline =  salineDF.loc[salineDF['time'] == 't8']
        time9Saline =  salineDF.loc[salineDF['time'] == 't9']
        
    #    pdb.set_trace()
        
        cnoDF = df.loc[df['dose'] == condition]
        time0CNO = cnoDF.loc[cnoDF['time'] == 't0']
        time1CNO = cnoDF.loc[cnoDF['time'] == 't1']
        time2CNO = cnoDF.loc[cnoDF['time'] == 't2']
        time3CNO = cnoDF.loc[cnoDF['time'] == 't3']
        time4CNO = cnoDF.loc[cnoDF['time'] == 't4']
        time5CNO = cnoDF.loc[cnoDF['time'] == 't5']
        time6CNO = cnoDF.loc[cnoDF['time'] == 't6']
        time7CNO = cnoDF.loc[cnoDF['time'] == 't7']
        time8CNO = cnoDF.loc[cnoDF['time'] == 't8']
        time9CNO = cnoDF.loc[cnoDF['time'] == 't9']
        
    #    pdb.set_trace()
    
        
        for eachBrainState in brainStateList:
            #SALINE
            time0SalineBrainState = time0Saline.loc[time0Saline['state'] == eachBrainState]
            t0SalineStatList = time0SalineBrainState[stat].tolist()
            
            time1SalineBrainState = time1Saline.loc[time1Saline['state'] == eachBrainState]
            t1SalineStatList = time1SalineBrainState[stat].tolist()
            
            time2SalineBrainState = time2Saline.loc[time2Saline['state'] == eachBrainState]
            t2SalineStatList = time2SalineBrainState[stat].tolist()
            
            time2SalineBrainState = time2Saline.loc[time2Saline['state'] == eachBrainState]
            t2SalineStatList = time2SalineBrainState[stat].tolist()
            
            time3SalineBrainState = time3Saline.loc[time3Saline['state'] == eachBrainState]
            t3SalineStatList = time3SalineBrainState[stat].tolist()
            
            time4SalineBrainState = time4Saline.loc[time4Saline['state'] == eachBrainState]
            t4SalineStatList = time4SalineBrainState[stat].tolist()
            
            time5SalineBrainState = time5Saline.loc[time5Saline['state'] == eachBrainState]
            t5SalineStatList = time5SalineBrainState[stat].tolist()
            
            time6SalineBrainState = time6Saline.loc[time6Saline['state'] == eachBrainState]
            t6SalineStatList = time6SalineBrainState[stat].tolist()
            
            time7SalineBrainState = time7Saline.loc[time7Saline['state'] == eachBrainState]
            t7SalineStatList = time7SalineBrainState[stat].tolist()
            
            time8SalineBrainState = time8Saline.loc[time8Saline['state'] == eachBrainState]
            t8SalineStatList = time8SalineBrainState[stat].tolist()
            
            time9SalineBrainState = time9Saline.loc[time9Saline['state'] == eachBrainState]
            t9SalineStatList = time9SalineBrainState[stat].tolist()
            
            
            #CNO
            time0CNOBrainState = time0CNO.loc[time0CNO['state'] == eachBrainState]
            t0CNOStatList = time0CNOBrainState[stat].tolist()
            
            time1CNOBrainState = time1CNO.loc[time1CNO['state'] == eachBrainState]
            t1CNOStatList = time1CNOBrainState[stat].tolist()
            
            time2CNOBrainState = time2CNO.loc[time2CNO['state'] == eachBrainState]
            t2CNOStatList = time2CNOBrainState[stat].tolist()
            
            time3CNOBrainState = time3CNO.loc[time3CNO['state'] == eachBrainState]
            t3CNOStatList = time3CNOBrainState[stat].tolist()
            
            time4CNOBrainState = time4CNO.loc[time4CNO['state'] == eachBrainState]
            t4CNOStatList = time4CNOBrainState[stat].tolist()
            
            time5CNOBrainState = time5CNO.loc[time5CNO['state'] == eachBrainState]
            t5CNOStatList = time5CNOBrainState[stat].tolist()
            
            time6CNOBrainState = time6CNO.loc[time6CNO['state'] == eachBrainState]
            t6CNOStatList = time6CNOBrainState[stat].tolist()
            
            time7CNOBrainState = time7CNO.loc[time7CNO['state'] == eachBrainState]
            t7CNOStatList = time7CNOBrainState[stat].tolist()
            
            time8CNOBrainState = time8CNO.loc[time8CNO['state'] == eachBrainState]
            t8CNOStatList = time8CNOBrainState[stat].tolist()
            
            time9CNOBrainState = time9CNO.loc[time9CNO['state'] == eachBrainState]
            t9CNOStatList = time9CNOBrainState[stat].tolist()
            
            
            #PVALUES
            ttest_indResultT0 = stats.ttest_ind(t0SalineStatList, t0CNOStatList)
            pvalueT0 = ttest_indResultT0[1]        
            pvalueActualT0 = pvalueT0 / 10 #Number of bins
            
            ttest_indResultT1 = stats.ttest_ind(t1SalineStatList, t1CNOStatList)
            pvalueT1 = ttest_indResultT1[1]        
            pvalueActualT1 = pvalueT1 / 10 #Number of bins
            
            ttest_indResultT2 = stats.ttest_ind(t2SalineStatList, t2CNOStatList)
            pvalueT2 = ttest_indResultT2[1]        
            pvalueActualT2 = pvalueT2 / 10 #Number of bins
            
            ttest_indResultT3 = stats.ttest_ind(t3SalineStatList, t3CNOStatList)
            pvalueT3 = ttest_indResultT3[1]        
            pvalueActualT3 = pvalueT3 / 10 #Number of bins
            
            ttest_indResultT4 = stats.ttest_ind(t4SalineStatList, t4CNOStatList)
            pvalueT4 = ttest_indResultT4[1]        
            pvalueActualT4 = pvalueT4 / 10 #Number of bins
            
            ttest_indResultT5 = stats.ttest_ind(t5SalineStatList, t5CNOStatList)
            pvalueT5 = ttest_indResultT5[1]        
            pvalueActualT5 = pvalueT5 / 10 #Number of bins
            
            ttest_indResultT6 = stats.ttest_ind(t6SalineStatList, t6CNOStatList)
            pvalueT6 = ttest_indResultT6[1]        
            pvalueActualT6 = pvalueT6 / 10 #Number of bins
            
            ttest_indResultT7 = stats.ttest_ind(t7SalineStatList, t7CNOStatList)
            pvalueT7 = ttest_indResultT7[1]        
            pvalueActualT7 = pvalueT7 / 10 #Number of bins
            
            ttest_indResultT8 = stats.ttest_ind(t8SalineStatList, t8CNOStatList)
            pvalueT8 = ttest_indResultT8[1]        
            pvalueActualT8 = pvalueT8 / 10 #Number of bins
            
            ttest_indResultT9 = stats.ttest_rel(t9SalineStatList, t9CNOStatList)
            pvalueT9 = ttest_indResultT9[1]        
            pvalueActualT9 = pvalueT9 / 10 #Number of bins
        
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T0', pvalueActualT0]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T1', pvalueActualT1]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T2', pvalueActualT2]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T3', pvalueActualT3]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T4', pvalueActualT4]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T5', pvalueActualT5]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T6', pvalueActualT6]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T7', pvalueActualT7]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T8', pvalueActualT8]], columns=labels))
            timebinPValuesDF = timebinPValuesDF.append(pd.DataFrame(data = [[stat, eachBrainState, 'T9', pvalueActualT9]], columns=labels))
            
#            pdb.set_trace()

#    pdb.set_trace()
    return timebinPValuesDF




#############Sleep_stats Helper Functions###############################
def swap_eeg_emg(ppath, rec):
    """
    swap EEG and EEG2 or EMG with EMG2 if $ch='EMG'
    """
    
    
    EEG = so.loadmat(os.path.join(ppath, rec, 'EEG.mat'))
    EMG = so.loadmat(os.path.join(ppath, rec,'EMG.mat'))
        
    tmp = EEG
    EEG = EMG
    EMG = tmp
    
    file_eeg = os.path.join(ppath, rec, 'EEG.mat')
    file_emg = os.path.join(ppath, rec, 'EMG.mat')
    so.savemat(file_eeg, {'EEG' : EEG})        
    so.savemat(file_emg, {'EMG' : EMG})



class Mouse :    
    def __init__(self, idf, list=None, typ='') :
        self.recordings = []
        self.recordings.append(list)
        self.typ = typ
        self.idf = idf
    
    def add(self, rec) :
        self.recordings.append(rec)

    def __len__(self) :
        return len(self.recordings)

    def __repr__(self) :
        return ", ".join(self.recordings)


def load_recordings(ppath, rec_file) :
    """
    load_recordings(ppath, rec_file)
    
    load recording listing with syntax:
    [E|C] \s+ recording_name
    
    #COMMENT
    
    @RETURN:
        (list of controls, lis of experiments)
    """
    exp_list = []
    ctr_list = []    

    rfile = os.path.join(ppath, rec_file)
    f = open(rfile, newline=None)
    lines = f.readlines()
    f.close()

    for l in lines :
        if re.search('^\s+$', l) :
            continue
        if re.search('^\s*#', l) :
            continue
        
        a = re.split('\s+', l)
        
        if re.search('E', a[0]) :
            exp_list.append(a[1])
            
        if re.search('C', a[0]) :
            ctr_list.append(a[1])
            
    return ctr_list, exp_list



def load_dose_recordings(ppath, rec_file):
    """
    load recording list with following syntax:
    A line is either control or experiments; Control recordings look like:
    C \s recording_name
    Experimental recordings also come with an additional dose parameter 
    (allowing for comparison of multiple doses with controls)
    
    E \s recording_name \s dose_1
    E \s recording_name \s dose_2
    """
    
    rfile = os.path.join(ppath, rec_file)
    f = open(rfile, newline=None)
    lines = f.readlines()
    f.close()

    # first get all potential doses
    doses = {}
    ctr_list = []
    for l in lines :
        if re.search('^\s+$', l):
            continue
        if re.search('^\s*#', l):
            continue        
        a = re.split('\s+', l)
        
        if re.search('E', a[0]):
            if a[2] in doses:
                doses[a[2]].append(a[1])
            else:
                doses[a[2]] = [a[1]]

        if re.search('C', a[0]):
            ctr_list.append(a[1])

    return ctr_list, doses
    
    

def load_stateidx(ppath, name, ann_name=''):
    """ load the sleep state file of recording (folder) $ppath/$name
    @Return:
        M,K         sequence of sleep states, sequence of 
                    0'1 and 1's indicating non- and annotated states
    """   
    
    if ann_name == '':
        ann_name = name
        
    file = os.path.join(ppath, name, 'remidx_' + ann_name + '.txt')
    
    f = open(file, 'r')    
    lines = f.readlines()
    f.close()
    
    n = 0
    for l in lines:
        if re.match('\d', l):
            n = n+1
            
    M = np.zeros(n)
    K = np.zeros(n)
    
    i = 0
    for l in lines :
        
        if re.search('^\s+$', l) :
            continue
        if re.search('\s*#', l) :
            continue
        
        if re.match('\d+\s+\d+', l) :
            a = re.split('\s+', l)
            M[i] = int(a[0])
            K[i] = int(a[1])
            i = i+1
            
    return M,K



def get_snr(ppath, name) :
    """
    read and return SR from file $ppath/$name/info.txt 
    """
    fid = open(os.path.join(ppath, name, 'info.txt'), 'rU')
    lines = fid.readlines()
    fid.close()
    values = []
    for l in lines :
        a = re.search("^" + 'SR' + ":" + "\s+(.*)", l)
        if a :
            values.append(a.group(1))            
    return float(values[0])




def get_sequences(idx, ibreak=1): #default get_sequences for use with sleepy
    #def get_sequences(idx, ibreak):
    """
    get_sequences(idx, ibreak=1)
    idx     -    np.vector of indices
    @RETURN:
    seq     -    list of np.vectors
    """
    diff = idx[1:] - idx[0:-1]
    breaks = np.nonzero(diff>ibreak)[0]
    breaks = np.append(breaks, len(idx)-1)
    
    seq = []    
    iold = 0
    for i in breaks:
        r = list(range(iold, i+1))
        seq.append(idx[r])
        iold = i+1
        
    return seq

### MANIPULATING FIGURES ##############################################################
def set_fontsize(fs):
    import matplotlib
    matplotlib.rcParams.update({'font.size': fs})



def save_figure(fig_file):
    import matplotlib
    # alternative way of setting nice fonts:
    
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42
    #matplotlib.pylab.savefig(fig_file, dpi=300)

    matplotlib.rcParams['text.usetex'] = False 
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.pylab.savefig(fig_file, bbox_inches="tight", dpi=300)
    matplotlib.rcParams['text.usetex'] = False   


def box_off(ax):
    ax.spines["top"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()
     

def power_spectrum(data, length, dt):
    """
    scipy's implementation of Welch's method using hanning window to estimate
    the power spectrum
    @Parameters
        data    -   time series; float vector!
        length  -   length of hanning window, even integer!
    
    @Return:
        power density, frequencies

        The function return power density in units V^2 / Hz
        Note that
            np.var(data) ~ np.sum(power density) * (frequencies[1]-frequencies[0])
    """
    f, pxx = scipy.signal.welch(data, fs=1/dt, window='hann', nperseg=int(length), noverlap=int(length/2))
    return pxx, f

def load_laser(ppath, name):
    """
    load laser from recording ppath/name ...
    @RETURN: 
    @laser, vector of 0's and 1's 
    """ 
    # laser might be .mat or h5py file
    # perhaps we could find a better way of testing that
    file = os.path.join(ppath, name, 'laser_'+name+'.mat')
    try:
        laser = np.array(h5py.File(file,'r').get('laser'))
    except:
        laser = so.loadmat(file)['laser']
    return np.squeeze(laser)
    
def laser_start_end(laser, SR=1525.88, intval=5):
    """laser_start_end(ppath, name) ...
    print start and end index of laser stimulation trains: For example,
    if you was stimulated for 2min every 20 min with 20 Hz, return the
    start and end index of the each 2min stimulation period (train)

    returns the tuple (istart, iend), both indices are inclusive,
    i.e. part of the sequence
    @Param:
    laser    -    laser, vector of 0s and 1s
    intval   -    minimum time separation [s] between two laser trains
    @Return:
    (istart, iend) - tuple of two np.arrays with laser start and end indices
    """
    idx = np.where(laser > 0.5)[0]
    if len(idx) == 0 :
        #return (None, None)
        return ([], [])
    
    idx2 = np.nonzero(np.diff(idx)*(1./SR) > intval)[0]
    istart = np.hstack([idx[0], idx[idx2+1]])
    iend   = np.hstack([idx[idx2], idx[-1]])    

    return (istart, iend)
    

