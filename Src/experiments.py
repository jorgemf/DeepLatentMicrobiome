import os
import gc

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import warnings

warnings.filterwarnings('ignore')

import tensorflow as tf

physical_devices = tf.config.list_physical_devices('GPU')
for gpu in physical_devices:
    tf.config.experimental.set_memory_growth(gpu, True)

from utils import *
from layers import *
from loss import *
from results import save_summary
from train import perform_experiment
from tensorflow.keras import optimizers
import multiprocessing
import tqdm
import numpy as np

np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)


def imap_perform_experiment(p):
    experiment_metrics, models, result = perform_experiment(*p)
    del models
    del result
    gc.collect()
    return experiment_metrics

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

if __name__ == '__main__':

    all_metrics = []

    cv_folds = 5
    epochs = 100
    optimizer = optimizers.Adam

    params = [
        (None, layers.ReLU, losses.MeanSquaredError),
        (CenterLogRatio, None, MakeLoss(LossMeanSquaredErrorWrapper, CenterLogRatio, None)),
        (
            Percentage, layers.Softmax,
            MakeLoss(LossCategoricalCrossentropyWrapper, Percentage, None)),
        (CenterLogRatio, layers.Softmax,
         MakeLoss(LossCategoricalCrossentropyWrapper, Percentage, None)),
        (Percentage, layers.Softmax, MakeLoss(LossBrayCurtis, Percentage, None)),
    ]
    all_params = []

    for input_transform, output_transform, reconstruction_loss in params:
        for batch_size, learning_rate, learning_rate_scheduler in [
            [128, 0.01, None],
            [64, 0.001, None],
            [128, 0.01, [ExpDecayScheluder(initial_lr=0.01, k=0.92), ('lr=0.01,k=0.92')]]]:
            for activation in ['tanh', 'relu', 'sigmoid']:
                activation_latent = activation
                for latent_space in [10, 50, 100]:
                    for autoencoder_layers in [[512, 256], [256], []]:
                        v = [cv_folds, epochs, batch_size, learning_rate, optimizer,
                             learning_rate_scheduler, input_transform, output_transform,
                             reconstruction_loss, latent_space,
                             autoencoder_layers, activation, activation_latent]
                        all_params.append(v)
    all_params = all_params
    print('Number of experiments: ', len(all_params))

    PROCESS_PER_GPU = 4

    physical_devices = tf.config.list_physical_devices('GPU')
    for gpu in physical_devices:
        try:
            tf.config.experimental.set_memory_growth(gpu, True)
        except:
            # Invalid device or cannot modify virtual devices once initialized.
            pass

    if len(physical_devices) == 0:
        physical_devices = ['/CPU:0']
    else:
        physical_devices = ['/GPU:{}'.format(i) for i in range(len(physical_devices))]

    print(physical_devices)

    GPUS = len(physical_devices)

    all_params_gpu = [x + [False, physical_devices[i % GPUS]] for i, x in enumerate(all_params)]

    multiprocessing.set_start_method('spawn', force=True)

    all_metrics = []

    with tqdm.tqdm(total=len(all_params_gpu)) as pbar:
        for b in batch(all_params_gpu, (PROCESS_PER_GPU-1) * 10):
            with multiprocessing.Pool(processes=PROCESS_PER_GPU * GPUS) as pool:
                for experiment_metrics in pool.imap(imap_perform_experiment, b):
                    all_metrics.append(experiment_metrics)
                    gc.collect()
                    pbar.update()
    save_summary(all_metrics, show=False)
