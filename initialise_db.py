from sqlalchemy import create_engine, Column, String, ForeignKey
from sqlalchemy.orm import declarative_base, relationship, sessionmaker
import pandas as pd

Base = declarative_base()
engine = create_engine('sqlite:///analysis.db')

class Sample(Base):
    __tablename__ = 'sample'
    id = Column('id', String, primary_key=True)
    pool_id = Column('pool_id', String)
    protocol = Column('protocol', String)
    patient_id = Column('patient_id', String, ForeignKey('patient.id'))
    
    def __init__(self, id):
        self.id = id
        metadata = pd.read_csv('../KCL-Raw-data/metadata.csv', delimiter='\t')
        self.pool_id = metadata['pool_id'][metadata['Sample']==id].values[0]
        self.protocol = metadata['protocol'][metadata['Sample']==id].values[0]
        self.patient_id = metadata['patient_id'][metadata['Sample']==id].values[0]

    def get_df(self):
        df = pd.read_csv('../KCL-Raw-data/' + self.id + '.tsv', delimiter='\t')
        df = df.rename(columns={'aminoAcid': 'aminoacid', 'count (templates/reads)': 'count'})
        df = df [['aminoacid', 'count']]
        df = df[df['aminoacid'].notna()]
        df['proportion'] = df['count']/sum(df['count'])
        return df
    
    def __repr__(self):
        return self.id
    
class Patient(Base):
    __tablename__ = 'patient'
    id = Column('id', String, primary_key=True)
    hormad = Column('hormad', String)
    genomics = Column('genomics', String)
    samples = relationship('Sample', backref='patient', lazy=True)

    def __init__(self, id):
        self.id = id
        metadata = pd.read_csv('../KCL-Raw-data/metadata.csv', delimiter='\t')
        self.hormad = metadata['HORMAD1'][metadata['patient_id']==id].values[0]
        self.genomics = metadata['Genomics'][metadata['patient_id']==id].values[0]
    
    def __repr__(self):
        return self.id
    
Base.metadata.create_all(bind=engine)
Session = sessionmaker(bind=engine)
session = Session()


def initialise_db():
    metadata = pd.read_csv('../KCL-Raw-data/metadata.csv', delimiter='\t')
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