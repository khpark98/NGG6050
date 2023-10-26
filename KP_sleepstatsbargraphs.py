# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:24:44 2022

Bar graphs of sleep stats -- duration, freq, % for NREM, REM, wake

@author: ChungWeberPC_19
"""

import sleepy
import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from statsmodels.stats.multicomp import (pairwise_tukeyhsd, MultiComparison)
import pingouin as pg



""" ppath to your text file (and recordings)"""

ppath = r'C:\\Users\\ChungWeberPC_23\\Desktop\\Kristen\\KP01_warmingexp'
textfile = 'KP01_filelistingEEselected.txt'
"""your recording and group list as a text file"""
data=sleepy.load_dose_recordings(ppath, textfile)[1]

"""if your recordings are in different folders, write that folders=True so you will be asked to input the ppath to each folder separately, else folders=False"""
folders=False
groups=['BL', 'warming']
states=['REM', 'Wake', 'NREM']
sleepstat=['Perc', 'Freq', 'Dur']
mathr=20 #or 20 #0 for NREM and 20 for wake
tstart=0*3600
tend=1.99*3600
xlabels=['Baseline', 'Warming']
timecourse=False

if timecourse==False:
    """creates dataframes for the sleep stats of each group"""
   
    dataframelist=[]
    
    for g in groups:
        groupfile=str(g)+'.csv'
        ##if folders==True:
            ##print(g)
            ##if g=='saline_c' or g=='cno1_c':
                ##ppathg=r'E:\LCretrodreadd\LCtoPOAretro-idreadd_DBH'
           
          
           
            ##elif g =='Esaline' or g=='cno1':
                ##ppathg=r'E:\LCretrodreadd\LCtoPOAretro-edreadd_DBH'
    
         
            ##ppath=ppathg
            
        sleepy.sleep_stats(ppath, data[g],mathr, tstart, tend, pplot=False, csv_file=os.path.join(ppath, groupfile))
        sleepy.plot_hypnograms(ppath, data[g], tstart=tstart, tend=tend, ma_thr=20)
       
      
        dataframe=pd.DataFrame(pd.read_csv(os.path.join(ppath, groupfile)))
   
        
        dataframe['group'] = g
     
      #  if g=='control':
      #      dataframe['genotype']='gad'
       
      #  elif  g=='' or g=='':
      #      dataframe['genotype']='16p'
          
        #if g=='saline_c' or g=='Esaline':
            #dataframe['treatment']=''
        #elif g=='cno1_c' or g=='cno1':
            #dataframe['treatment']='cno'

        dataframelist.append(dataframe)
        dataframe.to_csv(os.path.join(ppath, groupfile))
       
           
            
    alldata=pd.concat(dataframelist, ignore_index=True)

    print(alldata)
        
   #creating different combinations  of what sleep stats to plot
   
    for a in states:
        for b in sleepstat:
            tocompare=[a,b]
            print('Plotting next: ' + tocompare[0] + ' '+ tocompare[1])
    
   
        
            #melted=pd.melt(alldata[alldata.state==tocompare[0]], id_vars=['group', 'mouse', 'genotype', 'treatment'], value_vars=[tocompare[1]])
            melted=pd.melt(alldata[alldata.state==tocompare[0]], id_vars=['group', 'mouse'], value_vars=[tocompare[1]])
            plt.ion()
            sns.set_style('white')
            sns.set_style("ticks")
            dotcolors=['black','black','black','black','black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black']
            indmousecolor=sns.set_palette(sns.color_palette(dotcolors))
            #ax=sns.catplot(x='group', y='value', kind='swarm',hue='mouse', palette=indmousecolor, data=melted)   #plots each mouse separately
            #ax=sns.catplot(x='group', y='value', kind='swarm',hue='mouse', palette=[muted], data=melted)   #plots each mouse separately #[muted] each data point is a different color w palette muted
            ax=sns.relplot(x='group', y='value', hue='mouse', kind='line',sort=False,palette=indmousecolor,linewidth=0.5, aspect=0.7,data=melted)  #PLOTTING THE SAME MOUSE IN MULTIPLE GROUPS
                  
            colors=['lightgray', 'lightblue', 'rosybrown','rosybrown', 'lightblue', 'lightblue']
            mycolors=sns.set_palette(colors)
            
            ax=sns.barplot(x='group', y='value',palette=mycolors, ci=68, errwidth=2, capsize=0.05, data=melted) #ci=68 is SEM
            #ax=sns.swarmplot(x='group', y='value', data=melted, color='black', size=7) #plot individual data points
           
        
            ax.tick_params(labelsize=18) #changes the size of x and y tick labels
            if len(xlabels) == 0:
                ax.set_xticklabels([groups[0], groups[1]])
            if len(xlabels) ==2:
                ax.set_xticklabels([xlabels[0], xlabels[1]])
            if len(xlabels) ==3:
                ax.set_xticklabels([xlabels[0], xlabels[1], xlabels[2]])
            if len(xlabels) ==4:
                ax.set_xticklabels([xlabels[0], xlabels[1], xlabels[2], xlabels[3]],rotation=45)
            if len(xlabels) ==6:
                ax.set_xticklabels([xlabels[0], xlabels[1], xlabels[2], xlabels[3], xlabels[4], xlabels[5]],rotation=45)
            
            ax.set_xlabel(' ')
            ax.set_title(tocompare[0], fontsize=20)
            if tocompare[1]=='Freq':
                ax.set_ylabel('Frequency /h',fontsize=20)
            elif tocompare[1]=='Perc':
                ax.set_ylabel('%',fontsize=20)
            elif tocompare[1] == 'Dur':
                ax.set_ylabel('Duration (s)', fontsize=20)
            sns.set_style("ticks")
            plt.show()
           
            #statistics=pg.pairwise_ttests(data=melted, dv='value', within='group', subject='mouse') #paired 
            #pg.print_table(statistics)
            statistics=pg.pairwise_ttests(data=melted, dv='value', between='group', subject='mouse') #independent
            pg.print_table(statistics)
            #statistics=pg.mixed_anova(data=melted, dv='value', between='genotype',within='treatment', subject='mouse', correction=True)
            #description=('Statistics '+ tocompare[0] + ' '+ tocompare[1] + ': ' + 
            # '\n' + 'Info: mathr '+ str(mathr) +' tstart '+ str(tstart) + ' tend '+ str(tend))
        
            #print(description)
            #pg.print_table(statistics) 

            #ttesti=pg.pairwise_ttests(data=melted, dv='value', between='genotype',within='treatment', subject='mouse', padjust='bonf', within_first=False)
            #pg.print_table(ttesti)
            #ttesti=pg.pairwise_ttests(data=melted, dv='value', between='genotype',within='treatment', subject='mouse', padjust='bonf', within_first=True)
            #pg.print_table(ttesti)
           # statsit=pg.pairwise_ttests(data=melted,dv='value', subject='mouse')
           # pg.print_table(statsit)
            
          
           