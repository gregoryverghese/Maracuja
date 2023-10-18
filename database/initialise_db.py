from sqlalchemy import create_engine, Column, String, ForeignKey
from sqlalchemy.orm import declarative_base, relationship, sessionmaker
import pandas as pd

Base = declarative_base()
engine = create_engine('sqlite:///analysis.db')


class Sample(Base):
    __tablename__ = 'sample'
    id = Column('id', String, primary_key=True)
    pool = Column('pool', String)
    protocol = Column('protocol', String)
    patient_id = Column('patient_id', String, ForeignKey('patient.id'))
    data_path = Column('data_path', String)
    HLA_A = Column('HLA_A', String)
    
    
    def __init__(self, id):
        self.id = id        
        self.pool = metadata['pool'][metadata['Sample']==id].values[0]
        self.protocol = metadata['protocol'][metadata['Sample']==id].values[0]
        self.patient_id = metadata['patient_id'][metadata['Sample']==id].values[0]
        self.antigens = list(pool_peptides['Sequence'][pool_peptides['Pool']==self.pool]) 
        self.data_path = data_path
        self.HLA_A = metadata['HLA A'][metadata['Sample']==id].values[0]

    def get_df(self):
        df = pd.read_csv(self.data_path + '/' + self.id + '.tsv', delimiter='\t')
        df = df.rename(columns={'aminoAcid': 'b_aminoacid', 'count (templates/reads)': 'count'})
        df = df [['b_aminoacid', 'count']]
        df = df[df['b_aminoacid'].notna()]
        df['proportion'] = df['count']/sum(df['count'])
        df=df.sort_values('count', ascending = False).head(100).reset_index().drop('index', axis=1)

        antigens = list(pool_peptides['Sequence'][pool_peptides['Pool']==self.pool]) 
        df['pMTnet'] = [1]*len(df)
        for i in range(len(df)):
            all_pred = preds[preds['CDR3']==df['b_aminoacid'][i]]
            antigen_pred = all_pred[all_pred['Antigen'].isin(antigens)]
            pred = (antigen_pred['Rank'].min()+0.0001)**0.2
            df.at[i, 'pMTnet'] = pred

        return df
    
    def __repr__(self):
        return self.id
    
class Patient(Base):
    __tablename__ = 'patient'
    id = Column('id', String, primary_key=True)
    hormad = Column('hormad', String)
    genomics = Column('genomics', String)
    samples = relationship('Sample', backref='patient', lazy=True)
    sc_samples = relationship('SC_Sample', backref='patient', lazy=True)

    def __init__(self, id):
        self.id = id
        self.hormad = metadata['HORMAD1'][metadata['patient_id']==id].values[0]
        self.genomics = metadata['Genomics'][metadata['patient_id']==id].values[0]
    
    def __repr__(self):
        return self.id

class SC_Sample(Base):
    __tablename__ = 'sc_sample'
    id = Column('id', String, primary_key=True)
    pool = Column('pool', String)
    protocol = Column('protocol', String)
    patient_id = Column('patient_id', String, ForeignKey('patient.id'))
    data_path = Column('data_path', String)
    
    def __init__(self, id):
        self.id = id        
        self.pool = metadata['pool'][metadata['Sample']==id].values[0]
        self.patient_id = metadata['patient_id'][metadata['Sample']==id].values[0]
        self.protocol = metadata['protocol'][metadata['Sample']==id].values[0]
        self.data_path = data_path

    def get_df(self):
        df = pd.read_csv(self.data_path + '/' + self.id + '.csv', delimiter=',')
        df = df.rename(columns={'TCR_Alpha_Gamma_CDR3_Translation_Dominant': 'a_aminoacid', 'TCR_Beta_Delta_CDR3_Translation_Dominant': 'b_aminoacid'})
        
        # Fill missing alphas and betas
        alphas = df['a_aminoacid'][df['a_aminoacid'].notna()][df['TCR_Paired_Chains']==False]
        betas = df['b_aminoacid'][df['b_aminoacid'].notna()][df['TCR_Paired_Chains']==False]
        for i in alphas.index:
            partners = df['b_aminoacid'][df['a_aminoacid']==alphas[i]]
            missing_beta = list(set(partners[partners.notna()]))
            if missing_beta:
                df.at[i,'b_aminoacid'] = missing_beta[0]
        for i in betas.index:
            partners = df['a_aminoacid'][df['b_aminoacid']==betas[i]]
            missing_alphas = list(set(partners[partners.notna()]))
            if missing_alphas:
                df.at[i,'a_aminoacid'] = missing_alphas[0]
        
        df = df.groupby(['a_aminoacid','b_aminoacid'], dropna=0).size().reset_index()
        df.rename(columns={ df.columns[2]: "count" }, inplace = True)
        df = df[df['a_aminoacid'].notna()]
        df=df.sort_values('count', ascending = False).head(100).reset_index().drop('index', axis=1)
        df['proportion'] = df['count']/sum(df['count'])
        return df
    
    def __repr__(self):
        return self.id

def initialise_db(path):
    global metadata; global pool_peptides; global preds; global data_path
    data_path = path
    metadata = pd.read_csv(data_path + '/metadata.csv', delimiter='\t')
    pool_peptides = pd.read_csv(data_path + '/HORMAD1pool_withNetMHC.csv', delimiter=',')
    preds = pd.read_csv(data_path + '/pMTnet_prediction.csv', delimiter=',')

    patients = list(set(metadata['patient_id']))
    session.query(Patient).delete()
    for patient in patients:
        patient_entry = Patient(patient)
        session.add(patient_entry)


    samples = list(set(metadata['Sample']))
    session.query(Sample).delete()
    for sample in samples:
        sample_entry = Sample(sample)
        session.add(sample_entry)
    session.commit()


def initialise_SC_db(path):
    global metadata; global pool_peptides; global preds; global data_path
    data_path = path
    metadata = pd.read_csv(data_path + '/metadata.csv', delimiter='\t')
    patients = list(set(metadata['patient_id']))
    session.query(Patient).delete()
    for patient in patients:
        patient_entry = Patient(patient)
        session.add(patient_entry)

    samples = list(set(metadata['Sample']))
    session.query(SC_Sample).delete()
    for sample in samples:
        sample_entry = SC_Sample(sample)
        session.add(sample_entry)
    session.commit()



Base.metadata.create_all(bind=engine)
Session = sessionmaker(bind=engine)
session = Session()

