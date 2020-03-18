from tensorflow.keras import callbacks


class ExpDecayScheluder():

    def __init__(self, initial_lr=0.01, k=0.9):
        self.initial_lr = initial_lr
        self.k = k

    def make(self):
        exp_decay_fn = lambda epoch: self.initial_lr * self.k ** epoch
        return callbacks.LearningRateScheduler(exp_decay_fn)


class MakeLoss():
    def __init__(self, loss_fn, true_fn, pred_fn):
        self.loss_fn = loss_fn
        self.true_fn = true_fn
        self.pred_fn = pred_fn

    def make(self):
        true_fn = self.true_fn() if self.true_fn is not None else None
        pred_fn = self.pred_fn() if self.pred_fn is not None else None
        return self.loss_fn(true_fn, pred_fn)
