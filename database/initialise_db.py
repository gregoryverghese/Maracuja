from sqlalchemy import create_engine, Column, String, ForeignKey
from sqlalchemy.orm import declarative_base, relationship, sessionmaker
import pandas as pd

Base = declarative_base()
engine = create_engine('sqlite:///analysis.db')

class HO1_protein:
        def __init__(self):
                self.whole = 'MATAQLQRTPMSALVFPNKISTEHQSLVLVKRLLAVSVSCITYLRGIFPECAYGTRYLDDLCVKILREDKNCPGSTQLVKWMLGCYDALQKKYLRMVVLAVYTNPEDPQTISECYQFKFKYTNNGPLMDFISKNQSNESSMLSTDTKKASILLIRKIYILMQNLGPLPNDVCLTMKLFYYDEVTPPDYQPPGFKDGDCEGVIFEGEPMYLNVGEVSTPFHIFKVKVTTERERMENIDSTILSPKQIKTPFQKILRDKDVEDEQEHYTSDDLDIETKMEEQEKNPASSELEEPSLVCEEDEIMRSKESPDLSISHSQVEQLVNKTSELDMSESKTRSGKVFQNKMANGNQPVKSSKENRKRSQHESGRIVLHHFDSSSQESVPKRRKKFSEPKEHI'
                self.Pool1 = 'MATAQLQRTPMSALVFPNKISTEHQSLVLVKRLLAVSVSCITYLRGI'
                self.Pool2 = 'ITYLRGIFPECAYGTRYLDDLCVKILREDKNCPGSTQLVKWMLGCYD'
                self.Pool3 = 'WMLGCYDALQKKYLRMVVLAVYTNPEDPQTISECYQFKFKYTNNGPL'
                self.Pool4 = 'YTNNGPLMDFISKNQSNESSMLSTDTKKASILLIRKIYILMQNLGPL'
                self.Pool5 = 'MQNLGPLPNDVCLTMKLFYYDEVTPPDYQPPGFKDGDCEGVIFEGEP'
                self.Pool6 = 'VIFEGEPMYLNVGEVSTPFHIFKVKVTTERERMENIDSTILSPKQIK'
                self.Pool7 = 'LSPKQIKTPFQKILRDKDVEDEQEHYTSDDLDIETKMEEQEKNPASS'
                self.Pool8 = 'EKNPASSELEEPSLVCEEDEIMRSKESPDLSISHSQVEQLVNKTSEL'
                self.Pool9 = 'VNKTSELDMSESKTRSGKVFQNKMANGNQPVKSSKENRKRSQHESGR'
                self.Pool10 = 'SQHESGRIVLHHFDSSSQESVPKRRKKFSEPKEHI'
HO1 = HO1_protein()

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
        #self.antigens = list(pool_peptides['Sequence'][pool_peptides['Pool']==self.pool]) 
        self.data_path = data_path
        self.HLA_A = metadata['HLAA'][metadata['Sample']==id].values[0]

    def get_df(self):
        df = pd.read_csv(self.data_path + '/' + self.id + '.tsv', delimiter='\t')
        # df=df.sort_values('count', ascending = False).reset_index().drop('index', axis=1).head(100)

        #antigens = list(pool_peptides['Sequence'][pool_peptides['Pool']==self.pool]) 
        #df['pMTnet'] = [1]*len(df)
        #for i in range(len(df)):
            #all_pred = preds[preds['CDR3']==df['aminoacid'][i]]
            #antigen_pred = all_pred[all_pred['Antigen'].isin(antigens)]
            #pred = (antigen_pred['Rank'].min()+0.0001)**0.2
            #df.at[i, 'pMTnet'] = pred

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
        self.hormad = metadata['HORMAD1'][metadata['patient_id']==id].values[0]
        self.genomics = metadata['Genomics'][metadata['patient_id']==id].values[0]
    
    def __repr__(self):
        return self.id


def initialise_db(path):
    global metadata; global pool_peptides; global preds; global data_path
    data_path = path
    metadata = pd.read_csv(data_path + '/metadata.csv')
    print(metadata.columns)
    #pool_peptides = pd.read_csv(data_path + '/HORMAD1pool_withNetMHC.csv', delimiter=',')
    #preds = pd.read_csv(data_path + '/pMTnet_prediction.csv', delimiter=',')

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




Base.metadata.create_all(bind=engine)
Session = sessionmaker(bind=engine)
session = Session()

initialise_db('/Users/w2030634/CancerHub/TCR/Maracuja/KCL-Clean-data')


