import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from initialise_db import Patient, Sample, session, initialise_db

# def found_in_control():

sample = list(session.query(Sample).filter_by(id='KCL717_Pool10_TCRB'))[0]
df = sample.get_df()
expanded = df[df['proportion']>0.5]

# print(expanded['aminoacid'][0])

CTR = list(session.query(Sample).filter_by(id='KCL717_CEF_TCRB'))[0].pool_id
print(CTR)

# print(expanded['aminoacid'][0] in CTR['aminoacid'].to_list())
