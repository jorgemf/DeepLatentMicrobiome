from tensorflow.keras import losses


class LossMeanSquaredErrorWrapper(losses.MeanSquaredError):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(LossMeanSquaredErrorWrapper, self).__init__()
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def call(self, y_true, y_pred):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        return super(LossMeanSquaredErrorWrapper, self).call(y_true, y_pred)


class LossCategoricalCrossentropyWrapper(losses.CategoricalCrossentropy):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(LossCategoricalCrossentropyWrapper, self).__init__()
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def call(self, y_true, y_pred):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        return super(LossCategoricalCrossentropyWrapper, self).call(y_true, y_pred)


class LossBrayCurtis(losses.Loss):

    def __init__(self, y_true_transformer, y_pred_transformer):
        super(LossBrayCurtis, self).__init__(name='bray_curtis')
        self.y_true_transformer = y_true_transformer
        self.y_pred_transformer = y_pred_transformer

    def call(self, y_true, y_pred):
        if self.y_true_transformer is not None:
            y_true = self.y_true_transformer(y_true)
        if self.y_pred_transformer is not None:
            y_pred = self.y_pred_transformer(y_pred)
        diff = tf.math.abs(y_true - y_pred)
        sum = tf.math.abs(y_true + y_pred)
        value = tf.reduce_sum(diff, axis=-1) / tf.reduce_sum(sum, axis=-1)
        return value
