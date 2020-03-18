from layers import *
from metris import *
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import layers, losses, metrics, optimizers


def autoencoder(bioma_shape=717,
                domain_shape=36,
                output_shape=717,
                latent_space=10,
                bioma_layers=[128, 64],
                domain_layers=[32, 16],
                input_transform=CenterLogRatio(),
                output_transform=None,
                activation_function_encoder=layers.ReLU(),
                activation_function_decoder=layers.ReLU(),
                activation_function_latent='tanh',
                ):
    has_domain = domain_shape is not None
    has_bioma = bioma_shape is not None

    if not has_bioma and not has_domain:
        raise Exception('Either bioma or domain has to be expecified.')

    # encoder bioma

    if has_bioma:
        in_bioma = layers.Input(shape=(bioma_shape,), name='bioma_input_{}'.format(bioma_shape))
        net = in_bioma
        if input_transform is not None:
            net = input_transform(net)
        for s in bioma_layers:
            net = layers.Dense(s, activation=activation_function_encoder,
                               name="encoder_bioma_dense_{}".format(s))(net)
        encoded_bioma = layers.Dense(latent_space, activation=activation_function_latent,
                                     name='encoded_bioma_{}'.format(latent_space))(net)
        encoder_bioma = keras.Model(inputs=in_bioma, outputs=encoded_bioma, name='EncoderBioma')
    else:
        encoded_bioma = None
        encoder_bioma = None

    # encoder domain

    if has_domain:
        in_domain = layers.Input(shape=(domain_shape,), name='domain_input_{}'.format(domain_shape))
        net = in_domain
        for s in domain_layers:
            net = layers.Dense(s, activation=activation_function_encoder,
                               name="encoder_domain_dense_{}".format(s))(net)
        encoded_domain = layers.Dense(latent_space, activation=activation_function_latent,
                                      name='encoded_domain_{}'.format(latent_space))(net)
        encoder_domain = keras.Model(inputs=in_domain, outputs=encoded_domain, name='EncoderDomain')
    else:
        encoded_domain = None
        encoder_domain = None

    # decoder bioma for both autoencoders

    in_latent_space = layers.Input(shape=(latent_space,), name='latent_space_input')
    net = in_latent_space
    net_bioma = encoded_bioma
    net_domain = encoded_domain
    for s in reversed(bioma_layers):
        layer = layers.Dense(s, activation=activation_function_decoder,
                             name="decoder_dense_{}".format(s))
        net = layer(net)
        if has_bioma:
            net_bioma = layer(net_bioma)
        if has_domain:
            net_domain = layer(net_domain)

    layer = layers.Dense(output_shape, activation=None, name='decoded_bioma')
    decoded_bioma = layer(net)
    if has_bioma:
        net_bioma = layer(net_bioma)
    if has_domain:
        net_domain = layer(net_domain)

    if output_transform is not None:
        decoded_bioma = output_transform(decoded_bioma)
        if has_bioma:
            net_bioma = output_transform(net_bioma)
        if has_domain:
            net_domain = output_transform(net_domain)

    decoder_bioma = keras.Model(inputs=in_latent_space, outputs=decoded_bioma, name='DecoderBioma')

    # combined model for training

    if has_domain and has_bioma:
        diff_encoders = tf.math.abs(encoded_domain - encoded_bioma, name='diff_encoded')
        diff_encoders = Identity(name='latent')(diff_encoders)
        net_bioma = Identity(name='bioma')(net_bioma)
        net_domain = Identity(name='domain')(net_domain)

        model = keras.Model(inputs=[in_bioma, in_domain],
                            outputs=[net_bioma, net_domain, diff_encoders],
                            name='model')
    else:
        if has_bioma:
            net_bioma = Identity(name='bioma')(net_bioma)
            model = keras.Model(inputs=[in_bioma],
                                outputs=[net_bioma],
                                name='model')
        if has_domain:
            net_domain = Identity(name='domain')(net_domain)
            model = keras.Model(inputs=[in_domain],
                                outputs=[net_domain],
                                name='model')

    return model, encoder_bioma, encoder_domain, decoder_bioma
