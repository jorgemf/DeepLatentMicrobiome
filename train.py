import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import warnings

warnings.filterwarnings('ignore')

from metric import *
from utils import ExpDecayScheluder
from layers import *
from model import autoencoder
from results import print_results, plot_models
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, losses, metrics, optimizers, callbacks
from sklearn.model_selection import KFold
from IPython.display import Markdown, display
from tqdm.keras import TqdmCallback
from data import read_data

data_microbioma_train, data_microbioma_test, data_domain_train, data_domain_test, \
otu_columns, domain_columns = read_data()


def compile_train(model, encoder_bioma=None, encoder_domain=None,
                  reconstruction_error=losses.MeanSquaredError(),
                  encoded_comparison_error=losses.MeanAbsoluteError(),
                  metrics=[
                      [
                          metrics.MeanSquaredError(),
                          metrics.MeanAbsoluteError(),
                          metrics.MeanAbsolutePercentageError(),
                      ], [
                          metrics.MeanSquaredError(),
                          metrics.MeanAbsoluteError(),
                          metrics.MeanAbsolutePercentageError(),
                      ], [
                          metrics.MeanAbsoluteError(),
                      ]
                  ],
                  optimizer=optimizers.SGD(lr=0.01)):
    if encoder_domain is not None and encoder_bioma is not None:
        model.compile(optimizer=optimizer,
                      loss=[reconstruction_error, reconstruction_error, encoded_comparison_error],
                      metrics=metrics)
    elif encoder_bioma is not None:
        model.compile(optimizer=optimizer,
                      loss=reconstruction_error,
                      metrics=metrics[0])
    elif encoder_domain is not None:
        model.compile(optimizer=optimizer,
                      loss=reconstruction_error,
                      metrics=metrics[1])
    else:
        raise Exception('Not domain nor bioma models')


def train(model_fn,
          data_microbioma,
          data_domain,
          latent_space=10,
          folds=5,
          epochs=20,
          batch_size=128,
          learning_rate_scheduler=ExpDecayScheluder(),
          random_seed=347,
          verbose=0):
    data_zeros_latent = np.zeros((data_domain.shape[0], latent_space), dtype=data_domain.dtype)
    kf = KFold(n_splits=folds, random_state=random_seed, shuffle=True)
    results = []
    models = []
    train_callbacks = [
        callbacks.EarlyStopping(monitor='val_loss', patience=epochs + 1, restore_best_weights=True)]
    if verbose >= 0:
        train_callbacks += [TqdmCallback(verbose=verbose)]
    if learning_rate_scheduler is not None:
        train_callbacks += [learning_rate_scheduler.make()]

    tf.random.set_seed(random_seed)

    for train_index, test_index in kf.split(data_microbioma):
        m_train, m_test = data_microbioma[train_index], data_microbioma[test_index]
        d_train, d_test = data_domain[train_index], data_domain[test_index]
        z_train, z_test = data_zeros_latent[train_index], data_zeros_latent[test_index]
        all_models = model_fn()
        model, encoder_bioma, encoder_domain, decoder_bioma = all_models

        metrics_prefix = None
        if encoder_bioma is not None and encoder_domain is not None:
            x_train = (m_train, d_train)
            y_train = (m_train, m_train, z_train)
            x_test = (m_test, d_test)
            y_test = (m_test, m_test, z_test)
        elif encoder_bioma is not None:
            x_train = m_train
            y_train = m_train
            x_test = m_test
            y_test = m_test
            metrics_prefix = 'bioma'
        elif encoder_domain is not None:
            x_train = d_train
            y_train = m_train
            x_test = d_test
            y_test = m_test
            metrics_prefix = 'domain'

        train_dataset = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(5000).batch(
            batch_size)
        train_dataset = train_dataset.prefetch(tf.data.experimental.AUTOTUNE)
        val_dataset = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(batch_size)
        val_dataset = val_dataset.prefetch(tf.data.experimental.AUTOTUNE)
        r = model.fit(train_dataset,
                      epochs=epochs,
                      validation_data=val_dataset,
                      callbacks=train_callbacks,
                      verbose=0)
        if metrics_prefix is not None:
            old_keys = r.history
            r.history = {}
            for k, v in old_keys.items():
                if k == 'loss' or k == 'val_loss':
                    new_key = k
                elif k.startswith('val_'):
                    new_key = 'val_{}_{}'.format(metrics_prefix, k[4:])
                else:
                    new_key = '{}_{}'.format(metrics_prefix, k)
                r.history[new_key] = v
        del val_dataset
        del train_dataset
        del x_train
        del y_train
        del x_test
        del y_test
        results.append(r)
        models.append(all_models)
    return results, models


def get_experiment_metrics(input_transform, output_transform):
    name_in = input_transform.__class__.__name__ if input_transform is not None else ""
    name_out = output_transform.__class__.__name__ if output_transform is not None else ""

    relative_in_transform = None
    relative_out_transform = None

    relative_in_transform = Percentage()
    if name_in == "":
        relative_out_transform = Percentage()
    elif name_in == CenterLogRatio.__name__ and name_out != layers.Softmax.__name__:
        relative_out_transform = layers.Softmax()

    return [
        [
            MeanSquaredErrorWrapper(y_true_transformer=input_transform,
                                    y_pred_transformer=None),
            MeanAbsoluteErrorWrapper(y_true_transformer=input_transform,
                                     y_pred_transformer=None),
            MeanAbsolutePercentageErrorWrapper(y_true_transformer=relative_in_transform,
                                               y_pred_transformer=relative_out_transform),
            BrayCurtisDissimilarity(y_true_transformer=relative_in_transform,
                                    y_pred_transformer=relative_out_transform),
            PearsonCorrelation(y_true_transformer=relative_in_transform,
                               y_pred_transformer=relative_out_transform),
            # SpearmanCorrelation(y_true_transformer=relative_in_transform,
            #                    y_pred_transformer=relative_out_transform),
            JensenShannonDivergence(y_true_transformer=relative_in_transform,
                                    y_pred_transformer=relative_out_transform),
            # CrossEntropy(y_true_transformer=relative_in_transform,
            #             y_pred_transformer=relative_out_transform),
        ], [
            MeanSquaredErrorWrapper(y_true_transformer=input_transform,
                                    y_pred_transformer=None),
            MeanAbsoluteErrorWrapper(y_true_transformer=input_transform,
                                     y_pred_transformer=None),
            MeanAbsolutePercentageErrorWrapper(y_true_transformer=relative_in_transform,
                                               y_pred_transformer=relative_out_transform),
            BrayCurtisDissimilarity(y_true_transformer=relative_in_transform,
                                    y_pred_transformer=relative_out_transform),
            PearsonCorrelation(y_true_transformer=relative_in_transform,
                               y_pred_transformer=relative_out_transform),
            # SpearmanCorrelation(y_true_transformer=relative_in_transform,
            #                    y_pred_transformer=relative_out_transform),
            JensenShannonDivergence(y_true_transformer=relative_in_transform,
                                    y_pred_transformer=relative_out_transform),
            # CrossEntropy(y_true_transformer=relative_in_transform,
            #             y_pred_transformer=relative_out_transform),
        ], [
            metrics.MeanAbsoluteError(name='mae'),
        ]
    ]


def perform_experiment(cv_folds, epochs, batch_size, learning_rate, optimizer,
                       learning_rate_scheduler, input_transform, output_transform,
                       reconstruction_loss, latent_space, layers,
                       activation, activation_latent, show_results=False, device='/CPU:0'):
    if input_transform is not None:
        input_transform = input_transform()
    if output_transform is not None:
        output_transform = output_transform()
    if reconstruction_loss.__class__.__name__ == 'MakeLoss':
        reconstruction_loss = reconstruction_loss.make()
    else:
        reconstruction_loss = reconstruction_loss()
    domain_layers = [l // 16 for l in layers]
    bioma_autoencoder = " -> ".join(["b"] +
                                    [str(l) for l in layers] +
                                    [str(latent_space)] +
                                    [str(l) for l in reversed(layers)] +
                                    ["b"])
    domain_autoencoder = " -> ".join(["d"] +
                                     [str(l) for l in domain_layers] +
                                     [str(latent_space)] +
                                     [str(l) for l in reversed(layers)] +
                                     ["b"])
    in_transform_name = input_transform.__class__.__name__ if input_transform else "none"
    out_transform_name = output_transform.__class__.__name__ if output_transform else "none"
    lr_scheduler_text = learning_rate_scheduler[1] if learning_rate_scheduler is not None else "none"
    lr_text = learning_rate if learning_rate_scheduler is not None else "constant = {}".format(
        learning_rate)
    learning_rate_scheduler = learning_rate_scheduler[
        0] if learning_rate_scheduler is not None else None
    optimizer = optimizer(learning_rate=learning_rate)
    experiment_parameters = [
        ("Input transform", in_transform_name),
        ("Output transform", out_transform_name),
        ("Reconstruction Loss", reconstruction_loss.__class__.__name__),
        ("Latent Space", latent_space),
        ("Bioma Autoencoder", bioma_autoencoder),
        ("Domain Autoencoder", domain_autoencoder),
        ("Activation Encoder", activation),
        ("Activation Decoder", activation),
        ("Activation Latent", activation_latent),
        ("CV folds", cv_folds),
        ("Epochs", epochs),
        ("Batch Size", batch_size),
        ("Learning Rate Scheduler", lr_scheduler_text),
        ("Learning Rate", lr_text),
        ("Optimizer", optimizer.__class__.__name__),
    ]

    if show_results:
        md_text = ""
        md_text += "| Parameter             | Value         |\n"
        md_text += "|:----------------------|:--------------|\n"
        for n, v in experiment_parameters:
            md_text += "| {} | {} |\n".format(n, v)

        display(Markdown(md_text))

    def create_model(print_data=False):

        models = autoencoder(bioma_shape=717,
                             domain_shape=36,
                             output_shape=717,
                             latent_space=latent_space,
                             bioma_layers=layers,
                             domain_layers=domain_layers,
                             input_transform=input_transform,
                             output_transform=output_transform,
                             activation_function_encoder=activation,
                             activation_function_decoder=activation,
                             activation_function_latent=activation_latent)
        model, encoder_bioma, encoder_domain, decoder_bioma = models

        if print_data:
            plot_models(model, encoder_bioma, encoder_domain, decoder_bioma)
        compile_train(model,
                      encoder_bioma=encoder_bioma,
                      encoder_domain=encoder_domain,
                      reconstruction_error=reconstruction_loss,
                      encoded_comparison_error=losses.MeanAbsoluteError(),
                      metrics=get_experiment_metrics(input_transform, output_transform),
                      optimizer=optimizer)

        return model, encoder_bioma, encoder_domain, decoder_bioma

    create_model(print_data=False)

    with tf.device(device):
        results, models = train(create_model,
                                data_microbioma_train,
                                data_domain_train,
                                latent_space=latent_space,
                                folds=cv_folds,
                                epochs=epochs,
                                batch_size=batch_size,
                                learning_rate_scheduler=learning_rate_scheduler,
                                verbose=-1)

    validation_results = print_results(results, show_results=show_results)
    if show_results:
        display(Markdown("*************"))

    return experiment_parameters + validation_results, models, results
