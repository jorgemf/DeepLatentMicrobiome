import numpy as np
import pandas as pd
import random
import tensorflow.keras as keras
from sklearn.model_selection import train_test_split


def read_data(random_state=42):
    otu = pd.read_csv('data/otu_table_all_80.csv', index_col=0, header=None, sep='\t').T
    otu = otu.set_index('otuids')
    otu = otu.astype('int32')
    metadata = pd.read_csv('data/metadata_table_all_80.csv', sep='\t')
    metadata = metadata.set_index('X.SampleID')
    domain = metadata[['age',
                       'Temperature',
                       'Precipitation3Days',
                       'INBREDS',
                       'Maize_Line']]
    domain = pd.concat([domain, pd.get_dummies(domain['INBREDS'], prefix='INBREDS')], axis=1)
    domain = pd.concat([domain, pd.get_dummies(domain['Maize_Line'], prefix='Maize_Line')], axis=1)
    domain = domain.drop(['INBREDS', 'Maize_Line'], axis=1)
    df = pd.concat([otu, domain], axis=1, sort=True, join='outer')
    data_microbioma = df[otu.columns].to_numpy(dtype=np.float32)
    data_domain = df[domain.columns].to_numpy(dtype=np.float32)
    data_microbioma_train, data_microbioma_test, data_domain_train, data_domain_test = \
        train_test_split(data_microbioma, data_domain, test_size=0.1, random_state=random_state)
    return data_microbioma_train, data_microbioma_test, data_domain_train, data_domain_test, otu.columns, domain.columns


class DatasetSequence(keras.utils.Sequence):

    def __init__(self, data_microbioma, data_domain, idx, latent_space,
                 batch_size, shuffle, random_seed,
                 encoder_domain, encoder_bioma):
        self.idx = idx
        self.data_microbioma = data_microbioma
        self.data_domain = data_domain
        self.zeros = np.zeros((batch_size, latent_space), dtype=data_domain.dtype)
        self.batch_size = batch_size
        self.shuffle = shuffle
        random.seed(random_seed)
        self.encoder_domain = encoder_domain
        self.encoder_bioma = encoder_bioma
        self.on_epoch_end()

    def __len__(self):
        return int(np.ceil(len(self.idx) / float(self.batch_size)))

    def __getitem__(self, idx):
        idx_init = idx * self.batch_size
        idx_end = (idx + 1) * self.batch_size
        m = self.data_microbioma[self.idx[idx_init:idx_end]]
        d = self.data_domain[self.idx[idx_init:idx_end]]

        if self.encoder_bioma is not None and self.encoder_domain is not None:
            x = (m, d)
            y = (m, m, self.zeros)
        elif self.encoder_bioma is not None:
            x = m
            y = m
        elif self.encoder_domain is not None:
            x = d
            y = m
        return x, y

    def on_epoch_end(self):
        if self.shuffle:
            random.shuffle(self.idx)