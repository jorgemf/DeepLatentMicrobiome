import tensorflow as tf
from tensorflow.keras import metrics
import scipy

class MeanSquaredErrorWrapper(metrics.MeanSquaredError):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(MeanSquaredErrorWrapper, self).__init__(name='mse')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def __call__(self, y_true, y_pred, *args, **kwargs):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        return super(MeanSquaredErrorWrapper, self).__call__(y_true, y_pred, *args, **kwargs)


class MeanAbsoluteErrorWrapper(metrics.MeanAbsoluteError):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(MeanAbsoluteErrorWrapper, self).__init__(name='mae')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def __call__(self, y_true, y_pred, *args, **kwargs):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        return super(MeanAbsoluteErrorWrapper, self).__call__(y_true, y_pred, *args, **kwargs)


class MeanAbsolutePercentageErrorWrapper(metrics.MeanAbsolutePercentageError):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(MeanAbsolutePercentageErrorWrapper, self).__init__(name='mape')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def __call__(self, y_true, y_pred, *args, **kwargs):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        return super(MeanAbsolutePercentageErrorWrapper, self).__call__(y_true, y_pred, *args,
                                                                        **kwargs)


class BrayCurtisDissimilarity(metrics.Mean):
    def __init__(self, y_true_transformer, y_pred_transformer):
        super(BrayCurtisDissimilarity, self).__init__(name='BrayCurtis')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def update_state(self, y_true, y_pred, sample_weight=None):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        diff = tf.math.abs(y_true - y_pred)
        sum = tf.math.abs(y_true + y_pred)
        value = tf.reduce_sum(diff, axis=-1) / tf.reduce_sum(sum, axis=-1)
        return super(BrayCurtisDissimilarity, self).update_state(
            value, sample_weight=sample_weight)


class PearsonCorrelation(metrics.Mean):
    def __init__(self, y_true_transformer, y_pred_transformer):
        super(PearsonCorrelation, self).__init__(name='pearson_corr')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def update_state(self, y_true, y_pred, sample_weight=None):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        mean_y_true = tf.reduce_mean(y_true, axis=-1, keepdims=True)
        mean_y_pred = tf.reduce_mean(y_pred, axis=-1, keepdims=True)
        dev_y_true = y_true - mean_y_true
        dev_y_pred = y_pred - mean_y_pred
        l2_norm_y_true = dev_y_true * dev_y_true
        l2_norm_y_true = tf.sqrt(tf.reduce_sum(l2_norm_y_true, axis=-1))
        l2_norm_y_pred = dev_y_pred * dev_y_pred
        l2_norm_y_pred = tf.sqrt(tf.reduce_sum(l2_norm_y_pred, axis=-1))
        r = tf.reduce_sum(dev_y_true * dev_y_pred, axis=-1) / (l2_norm_y_true * l2_norm_y_pred)
        r = tf.where(tf.math.is_nan(r), tf.zeros_like(r), r) # Avoid nan
        return super(PearsonCorrelation, self).update_state(
            r, sample_weight=sample_weight)


class SpearmanCorrelation(metrics.Mean):
    def __init__(self, y_true_transformer, y_pred_transformer):
        super(SpearmanCorrelation, self).__init__(name='spearman_corr')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def spearmanr(self, y_true, y_pred):
        values = []
        for x, y in zip(y_true, y_pred):
            value, _ = scipy.stats.spearmanr(x, y)
            values.append(value)
        return value

    def update_state(self, y_true, y_pred, sample_weight=None):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        value = tf.py_function(self.spearmanr, [y_pred, y_true], Tout=tf.float32)
        return super(SpearmanCorrelation, self).update_state(
            value, sample_weight=sample_weight)


class JensenShannonDivergence(metrics.Mean):
    def __init__(self, y_true_transformer, y_pred_transformer):
        super(JensenShannonDivergence, self).__init__(name='jensen_shannon_divergence')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def kl_divergence(self, p, q):
        r = p * tf.math.log(p / q)
        return tf.reduce_sum(r, axis=-1)

    def update_state(self, y_true, y_pred, sample_weight=None):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        kl_1 = self.kl_divergence(y_true, y_pred)
        kl_2 = self.kl_divergence(y_pred, y_true)
        value = (kl_1 + kl_2) / 2
        return super(JensenShannonDivergence, self).update_state(
            value, sample_weight=sample_weight)


class CrossEntropy(metrics.Mean):
    def __init__(self, y_true_transformer, y_pred_transformer):
        super(CrossEntropy, self).__init__(name='crossentropy')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def update_state(self, y_true, y_pred, sample_weight=None):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        r = y_true * tf.math.log(y_pred)
        value = -tf.reduce_sum(r, axis=-1)
        return super(CrossEntropy, self).update_state(
            value, sample_weight=sample_weight)
