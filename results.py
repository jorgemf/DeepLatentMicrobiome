import pandas as pd
import numpy as np
import tensorflow.keras as keras
import matplotlib.pyplot as plt
from IPython.display import SVG, Markdown, display, display_pretty


def save_plot_model(model, name='model'):
    filename = '{}.png'.format(name)
    keras.utils.plot_model(model, show_shapes=True, show_layer_names=True, to_file=filename)


def plot_models(model, encoder_bioma, encoder_domain, decoder_bioma):
    def _print(t):
        display_pretty(t, raw=True)

    display(Markdown("### Encoder Bioma"))
    encoder_bioma.summary(line_length=120, print_fn=_print)
    display(SVG(
        keras.utils.model_to_dot(encoder_bioma, show_shapes=True).create(prog='dot', format='svg')))
    if encoder_domain is not None:
        display(Markdown("### Encoder Domain"))
        encoder_domain.summary(line_length=120, print_fn=_print)
        display(SVG(keras.utils.model_to_dot(encoder_domain, show_shapes=True).create(prog='dot',
                                                                                      format='svg')))
    display(Markdown("### Decoder Bioma"))
    decoder_bioma.summary(line_length=120, print_fn=_print)
    display(SVG(
        keras.utils.model_to_dot(decoder_bioma, show_shapes=True).create(prog='dot', format='svg')))
    display(Markdown("### Model for training"))
    model.summary(line_length=120, print_fn=_print)
    display(SVG(keras.utils.model_to_dot(model, show_shapes=True).create(prog='dot', format='svg')))


def print_results(results, show_results=True):
    sequences = {}
    for k in results[0].history.keys():
        arrays = np.array([r.history[k] for r in results])
        sequences[k] = [arrays.mean(axis=0), arrays.min(axis=0), arrays.max(axis=0)]

    # best loss
    min_loss_i = np.argmin(sequences['val_loss'][0])
    best_iteration = [('best_lost_epoch', min_loss_i)]
    for k in results[0].history.keys():
        if k.startswith('val_'):
            best_iteration.append(
                (k, [sequences[k][x][min_loss_i] for x in range(len(sequences[k]))]))
    if not show_results:
        return best_iteration
    display(Markdown("<p>Best iteration: <b>{}</b></p>".format(min_loss_i)))
    keys = [k for k in results[0].history.keys() if not k.startswith('val_')]
    x = results[0].epoch

    for k in keys:
        mean_seq, min_seq, max_seq = sequences[k]
        display(
            Markdown("<b>{}</b>: {:.4f} (min: {:.4f}, max: {:.4f})".format(k, mean_seq[min_loss_i],
                                                                           min_seq[min_loss_i],
                                                                           max_seq[min_loss_i])))
        plt.fill_between(x, min_seq, max_seq, color="blue", alpha=0.2)
        plt.plot(x, mean_seq, color="blue", alpha=0.8)
        k_val = 'val_{}'.format(k)
        if k_val in sequences:
            mean_seq, min_seq, max_seq = sequences[k_val]
            plt.fill_between(x, min_seq, max_seq, color="red", alpha=0.2)
            plt.plot(x, mean_seq, color="red", alpha=0.8)

        plt.title(k)
        plt.ylabel(k)
        plt.xlabel('Epoch')
        if k_val in sequences:
            plt.legend(['train', 'validation'])
        plt.show()

    md_text = "| Metric           | Mean    | Min     | Max     |\n"
    md_text += "|:-----------------|--------:|--------:|--------:|\n"
    for k in keys:
        mean_seq, min_seq, max_seq = sequences[k]
        md_text += "| {} | {} | {} | {} |\n".format(k, mean_seq[min_loss_i],
                                                    min_seq[min_loss_i],
                                                    max_seq[min_loss_i])

    display(Markdown(md_text))

    return best_iteration

def print_results_noEnsemble(results, show_results=True):
    sequences = {}
    for k in results.history.keys():
        arrays = np.array(results.history[k])
        sequences[k] = [arrays.mean(), arrays.min(), arrays.max()]
    # best loss
    min_loss_i = np.argmin(np.array(results.history['val_loss']))
    best_iteration = [('best_lost_epoch', min_loss_i)]

    display(Markdown("<p>Best iteration: <b>{}</b></p>".format(min_loss_i)))
    keys = [k for k in results.history.keys() if not k.startswith('val_')]
    x = results.epoch

    md_text = "| Metric           | Mean    | Min     | Max     |\n"
    md_text += "|:-----------------|--------:|--------:|--------:|\n"
    for k in keys:
        mean_seq, min_seq, max_seq = sequences[k]
        md_text += "| {} | {} | {} | {} |\n".format(k, mean_seq,
                                                    min_seq,
                                                    max_seq)
    display(Markdown(md_text))


def save_summary(all_metrics, filename='results.csv', show=True):
    keys_set = set()
    keys_list = []
    for m in all_metrics:
        for k, v in m:
            if k not in keys_set:
                keys_set.add(k)
                keys_list.append(k)

    data_results = []
    data_columns = ['experiment n.']

    md_text = "| experiment n. "
    for k in keys_list:
        if k.startswith('val_'):
            md_text += "| {} ".format(k[4:])
            data_columns.append(k[4:])
            data_columns.append('{}_min'.format(k[4:]))
            data_columns.append('{}_max'.format(k[4:]))
        else:
            md_text += "| {} ".format(k)
            data_columns.append(k)
    md_text += "|\n"
    md_text += "|:--"
    for k in keys_list:
        if k.startswith('val_'):
            md_text += "|--:"
        else:
            md_text += "|:--"
    md_text += "|\n"

    for idx, m in enumerate(all_metrics):
        data_row = [idx]
        md_text += "| {} ".format(idx)
        i = 0
        for k, v in m:
            while k != keys_list[i]:
                md_text += "| "
                i += 1
                data_row.append(None)
            i += 1
            if k.startswith('val_'):
                md_text += "| **{:.4f}** [{:.4f},{:.4f}] ".format(*v)
                data_row.extend(v)
            else:
                md_text += "| {} ".format(v)
                data_row.append(v)
        md_text += "|\n"
        data_results.append(data_row)

    if show:
        display(Markdown(md_text))

    df_results = pd.DataFrame(data_results, columns=data_columns)
    df_results = df_results.set_index(data_columns[0])
    df_results.to_csv(filename)
