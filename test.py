import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import warnings

warnings.filterwarnings('ignore')

from train import *

def test_experiment(models, input_transform, output_transform, show=True):
    data_bioma_test_transformed = data_microbioma_test
    if input_transform is not None:
        input_transform = input_transform()
        data_bioma_test_transformed = input_transform(data_microbioma_test)
    if output_transform is not None:
        output_transform = output_transform()
    metrics_results = {}
    metrics = get_experiment_metrics(input_transform, output_transform)[1][3:]
    otus_errors = []
    for cv_models in models:
        _, _, encoder_domain, decoder_bioma = cv_models
        encoded = encoder_domain.predict(data_domain_test)
        decoded = decoder_bioma.predict(encoded)
        for m in metrics:
            if m.name not in metrics_results:
                metrics_results[m.name] = []
            result = m(data_microbioma_test, decoded)
            m.reset_states()
            metrics_results[m.name].append(result.numpy())
        # otus error
        se = tf.math.squared_difference(decoded, data_bioma_test_transformed)
        mse = tf.reduce_mean(se, axis=0)
        otus_errors.append(mse)
    mse_otus = tf.reduce_mean(tf.stack(otus_errors, axis=0), axis=0)
    mse_otus_keys = sorted(zip(mse_otus.numpy(), otu_columns), key=lambda x: x[0])
    for k, v in list(metrics_results.items()):
        v = np.asarray(v)
        metrics_results[k] = (v.mean(), v.min(), v.max())
    if show:
        md_text = "## Test results \n"
        md_text += "| Metric           | Mean    | Min     | Max     |\n"
        md_text += "|:-----------------|--------:|--------:|--------:|\n"
        for k, v in metrics_results.items():
            md_text += "| {} | {} | {} | {} |\n".format(k, v[0], v[1], v[2])

        display(Markdown(md_text))

        md_text ="### Best Otus\n"
        md_text += "| OTU | mse |\n"
        md_text += "|:----|----:|\n"
        for v, k in mse_otus_keys[:10]:
            md_text += "| {} | {} |\n".format(k, v)
        md_text += "\n\n"
        md_text +="### Worst Otus\n"
        md_text += "| OTU | mse |\n"
        md_text += "|:----|----:|\n"
        for v, k in reversed(mse_otus_keys[-10:]):
            md_text += "| {} | {} |\n".format(k, v)

        display(Markdown(md_text))

    return metrics_results


def perform_test_experiment(cv_folds, epochs, batch_size, learning_rate, optimizer,
                           learning_rate_scheduler, input_transform, output_transform,
                           reconstruction_loss, latent_space, layers,
                           activation, activation_latent, show_results=False, device='/CPU:0'):
    summary, models, results = perform_experiment(cv_folds, epochs, batch_size, learning_rate,
                                                  optimizer, learning_rate_scheduler,
                                                  input_transform, output_transform,
                                                  reconstruction_loss, latent_space, layers,
                                                  activation, activation_latent,
                                                  show_results, device)
    metrics_test = test_experiment(models, input_transform=CenterLogRatio, output_transform=None)
    return models, metrics_test