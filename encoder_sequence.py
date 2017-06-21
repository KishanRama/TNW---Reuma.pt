# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 17:05:25 2017

@author: kisha_000
"""

#import numpy as np
#import pandas as pd

def encode(df):

    #dictionary used to replace n_bio_corrente by alphabetical order letters
    di = {1: "A", 2: "B", 3:"C", 4:"D", 5:"E", 6:"F", 7:"G", 8:"H", 9:"I", 10:"J"}
    
    #read .csv file
    #df = pd.read_csv('C:/Users/kisha_000/Desktop/tese/Temporal Needleman Wunsch/bio_corrente_sequence.csv',sep = ';')
    
    #remove duplicates [id_doente,n_bio_corrente]
    df.drop_duplicates(['id_doente','n_bio_corrente'],inplace = True)
    
    #remove rows with NaN values
    df.dropna(inplace = True)
    
    #reset index to make from 0 to len(dataframe)
    df = df.reset_index(drop=True)
    
    #replace n_bio_corrente int values by alphatebical letters
    df.replace({"n_bio_corrente": di}, inplace = True)
    
    #create new column that adds column n_bio_corrente and n_dias_duracao
    #this is an auxiliar step to be able to encode our sequences
    df['aux_encode'] = df['n_bio_corrente'] + ',' + df['n_dias_duracao'].astype(int).astype(str)

    #groupby id_doente and create our sequences - now we have a dataframe with id_doente and encode (encoded sequence)
    df = df.groupby('id_doente')['aux_encode'].apply(lambda x: '.'.join(x)).reset_index()
    
    df['aux_encode'] = '0.' + df['aux_encode'] + '.Z'
    
    return df