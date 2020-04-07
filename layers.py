import tensorflow as tf
from tensorflow.keras import layers


class CenterLogRatio(layers.Layer):

    def __init__(self, calculate_percentage=True, epsilon=1E-6):
        super(CenterLogRatio, self).__init__(name='center_log_ratio')
        self.calculate_percentage = calculate_percentage
        self.epsilon = epsilon

    def build(self, input_shape):
        super(CenterLogRatio, self).build(input_shape)

    def call(self, x, **kwargs):
        epsilon = tf.constant(self.epsilon)
        if self.calculate_percentage:
            x_sum = tf.math.reduce_sum(x, axis=-1, keepdims=True)
            elements = tf.shape(x)[-1]
            x_sum += epsilon * tf.cast(elements, dtype=tf.float32)
            x = tf.math.divide(x, x_sum, name='percentage')
        x += epsilon
        logx = tf.math.log(x)
        log_mean = tf.math.reduce_mean(logx, axis=-1, keepdims=True)
        return tf.subtract(logx, log_mean, name='center_log_ratio')

    def compute_output_shape(self, input_shape):
        return input_shape

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'calculate_percentage': self.calculate_percentage,
            'epsilon': self.epsilon,
        })
        return config


class Percentage(layers.Layer):

    def __init__(self, epsilon=1E-6):
        super(Percentage, self).__init__(name='percentage')
        self.epsilon = epsilon

    def build(self, input_shape):
        super(Percentage, self).build(input_shape)

    def call(self, x, **kwargs):
        epsilon = tf.constant(self.epsilon)
        x_sum = tf.math.reduce_sum(x, axis=-1, keepdims=True)
        elements = tf.shape(x)[-1]
        x_sum += epsilon * tf.cast(elements, dtype=tf.float32)
        percent = tf.math.divide(x, x_sum, name='percentage')
        percent += epsilon
        return percent

    def compute_output_shape(self, input_shape):
        return input_shape

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'epsilon': self.epsilon,
        })
        return config


class Identity(layers.Layer):

    def __init__(self, name):
        super(Identity, self).__init__(name=name)

    def build(self, input_shape):
        super(Identity, self).build(input_shape)

    def call(self, x, **kwargs):
        return x

    def compute_output_shape(self, input_shape):
        return input_shape
