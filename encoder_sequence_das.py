# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 11:09:36 2017

@author: kisha_000
"""
import pandas as pd 
import numpy as np

def encode(df):

    
    #read .csv file
    #df = pd.read_csv('C:/Users/kisha_000/Desktop/tese/Temporal Needleman Wunsch/das28_sequence.csv',sep = ';',decimal=',')
    
    #remove duplicates [id_doente,dt_consulta]
    df.drop_duplicates(['id_doente','dt_consulta'],inplace = True)
    
    #remove rows with NaN values
    df.dropna(inplace = True)
    
    #reset index to make from 0 to len(dataframe)
    df = df.reset_index(drop=True)
    
    #convert dt_consulta column into datetime column
    df['dt_consulta'] = pd.to_datetime(df['dt_consulta'],dayfirst=True)
    
    #create auxily column with time intervals between appointments
    df['time_intervals'] = df['dt_consulta'].sub(df['dt_consulta'].shift())
    
    #convert previous result time_intervals in floats
    df['time_intervals'] = df['time_intervals'] / np.timedelta64(1, 'D')

    #set the time interval of the first appointment of each patient to zero
    df.loc[df.groupby('id_doente',as_index=False).head(1).index,'time_intervals'] = 0

    #create new column with das28_4v values replaced into 4 different events A,B,C,D
    #remission das28<2.6
    df.loc[df['das28_4v'] < 2.6, 'event'] = 'A'
    #low disease activity 2.6<=das28<3.2
    df.loc[(df['das28_4v'] <3.2) & (df['das28_4v']>=2.6), 'event'] = 'B'    
    #medium disease activity 3.2<=das28<=5.1
    df.loc[(df['das28_4v'] <=5.1) & (df['das28_4v']>=3.2), 'event'] = 'C'
    #high disease activity
    df.loc[df['das28_4v'] > 5.1, 'event'] = 'D'
    
    #create new column that adds column time_intervals and event 
    #this is an auxiliar step to be able to encode our sequences
    df['aux_encode'] = df['time_intervals'].astype(int).astype(str) + '.' + df['event']

    #groupby id_doente and create our sequences - now we have a dataframe with id_doente and encode (encoded sequence)
    df= df.groupby('id_doente')['aux_encode'].apply(lambda x: ','.join(x)).reset_index()
    
    #Add event Z at the end for all sequences
    #df['aux_encode'] = df['aux_encode'] + ',100.Z'
    return df