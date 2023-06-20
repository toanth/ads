import argparse
import json
from argparse import ArgumentParser
from itertools import groupby
from os.path import realpath, dirname
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

basepath = dirname(realpath(__file__)) + "/"


def color_y_axis(ax, color):
    for t in ax.get_yticklabels():
        t.set_color(color)


def is_quick_run(all_benchmarks):
    return all_benchmarks[0][0]['name'] == all_benchmarks[0][0]['run_name']


def check_same_config(data, bms, baseline_data, baseline_bms):
    config = data['context']
    baseline_config = baseline_data['context']
    for value in ['host_name', 'num_cpus', 'mhz_per_cpu', 'cpu_scaling_enabled', 'library_build_type']:
        if config[value] != baseline_config[value]:
            print(
                f"Warning: {value} differs between current ({config[value]}) and baseline ({baseline_config[value]}),"
                "comparisons may not be meaningful")
    if is_quick_run(bms):
        print("Current results are from a quick run, measurements may be more noisy than usual")
    if is_quick_run(baseline_bms):
        print("Baseline was a quick run, measurements may be more noisy than usual")


def print_avg_variation(benchmarks, quick_run):
    if quick_run:
        print("results are from a quick run, so there is no variation info available")
    else:
        cvs = [float(bm['cpu_time']) for _, runs in benchmarks.items() for bm in runs if bm['aggregate_name'] == 'cv']
        variation = np.mean(cvs)
        print(f"average coefficient of variation (ie standard deviation / mean): {variation * 100:.2f}%")


def get_benchmarks(data):
    benchmarks = data['benchmarks']
    res = {}
    for idx, it in groupby(benchmarks, key=lambda bm: bm['family_index']):
        res[idx] = list(it)
    return res


def get_runs(all_benchmarks, family_index, quick_run):
    if quick_run:
        if family_index >= len(all_benchmarks):
            return None, None
        res = [bm for bm in all_benchmarks[family_index] if bm['run_type'] == 'iteration']
        return res, np.array([[bm['cpu_time'] for bm in res], [0] * len(res)])
    else:
        if not family_index in all_benchmarks:
            return None, None
    info = all_benchmarks[family_index]
    runs = [bm for bm in info if bm['aggregate_name'] == 'median']
    stddevs = np.array([bm['cpu_time'] for bm in info if bm['aggregate_name'] == 'stddev'])
    means = np.array([bm['cpu_time'] for bm in info if bm['aggregate_name'] == 'mean'])
    error_bars = np.array([means, stddevs])
    max_variation = max([float(bm['cpu_time']) for bm in info if bm['aggregate_name'] == 'cv'])
    if max_variation > 0.05:
        print("Warning: High maximum measurement variation of "
              + str(max_variation * 100) + " percent for " + get_name(runs))
    return runs, error_bars


def get_name(runs):
    if runs is None:
        return "<None>"
    name = runs[0]['run_name'].split('/')[0].split('<')[0]
    if name[:3] == "BM_":
        name = name[3:]
    return name


def get_repetitions(runs):
    return int(runs[0]['repetitions'])


def get_performance_measurements(runs, errors, force_relative=None):
    assert len(runs) > 1
    divide_by_n = 'perN' in runs[0]
    if force_relative is not None:
        divide_by_n = force_relative
    num_iters = [int(benchmark['run_name'].split('/')[-1]) for benchmark in runs]
    times = np.array([benchmark['cpu_time'] for benchmark in runs])
    if divide_by_n:
        times = times / num_iters
        errors = errors / num_iters  # mean and stddev scale with constant factors
    bits = np.array([benchmark['bits'] for benchmark in runs])
    bits = bits / num_iters
    assert (times.shape == bits.shape == errors[0].shape == errors[1].shape == (len(runs),))
    return times, bits, num_iters, errors, divide_by_n


def generate_plots(data, baseline_data):
    num_skipped = 0
    baseline_benchmarks = None
    baseline_runs = None
    baseline_errors = None
    all_benchmarks = get_benchmarks(data)
    if not all_benchmarks:
        raise ValueError("No benchmarks were run")
    num_families = len(all_benchmarks)
    quick_run = is_quick_run(all_benchmarks)
    print_avg_variation(all_benchmarks, quick_run)
    if baseline_data is not None:
        baseline_benchmarks = get_benchmarks(baseline_data)
        check_same_config(data, all_benchmarks, baseline_data, baseline_benchmarks)

    for i in range(0, num_families):
        runs, errors = get_runs(all_benchmarks, i, quick_run)
        name = get_name(runs)
        if baseline_data is not None:
            baseline_runs, baseline_errors = get_runs(baseline_benchmarks, i - num_skipped,
                                                      is_quick_run(baseline_benchmarks))
            if get_name(baseline_runs) != name:
                print("Warning: No matching baseline results found for " + name)
                num_skipped = num_skipped + 1
                baseline_runs = None
        times, bits, num_iters, errors, per_n = get_performance_measurements(runs, errors)
        if baseline_runs is not None:
            baseline_times, baseline_bits, baseline_num_iters, baseline_errors, _ \
                = get_performance_measurements(baseline_runs, baseline_errors, per_n)
            assert baseline_num_iters == num_iters
            fig, ax = plt.subplots(ncols=2, figsize=(10, 5), layout='constrained')
            ax[0].plot(num_iters, times, 'b-s', label='current (median)')
            if get_repetitions(runs) > 1:
                ax[0].errorbar(num_iters, errors[0], errors[1], fmt='b.', label='current µ ± σ')
            ax[0].plot(baseline_num_iters, baseline_times, 'g-o', label='baseline (median)')
            if get_repetitions(baseline_runs) > 1:
                ax[0].errorbar(baseline_num_iters, baseline_errors[0], baseline_errors[1], fmt='g.',
                               label='baseline µ ± σ')
            ax[0].set_xlabel('n')
            ax[0].set_xscale('log')
            if per_n:
                ax[0].set_ylabel('time in ' + runs[0]['time_unit'] + ' divided by n')
            else:
                ax[0].set_ylabel('time in ' + runs[0]['time_unit'])
            ax[0].set_yscale('log')
            ax[0].grid()
            title = "time"
            if per_n:
                title += " (divided by n)"
            ax[0].set_title(title)
            rel_time = ax[0].twinx()
            time_percent = [100 * times[i] / baseline_times[i] for i in range(len(times))]
            color_y_axis(rel_time, 'c')
            rel_time.set_ylabel('relative time in percent')
            rel_time.plot(num_iters, time_percent, 'c-p', label='percent of baseline')
            rel_time.axhline(y=100, color='c', linestyle='--')

            ax[1].plot(num_iters, bits, 'b-s')
            ax[1].plot(baseline_num_iters, baseline_bits, 'g-o')
            ax[1].set_xlabel('n')
            ax[1].set_xscale('log')
            ax[1].set_ylabel('number of bits divided by n')
            ax[1].set_yscale('log')
            ax[1].grid()
            ax[1].set_title("space (divided by n)")
            ax[1].set_ylim([None, 1000])
            rel_space = ax[1].twinx()
            time_percent = [100 * bits[i] / baseline_bits[i] for i in range(len(bits))]
            rel_space.plot(num_iters, time_percent, 'c-p')
            rel_space.axhline(y=100, color='c', linestyle='--')
            color_y_axis(rel_space, 'c')
            rel_space.set_ylabel('relative space in percent')
            plt.suptitle(name + f"\n#repetitions: {get_repetitions(runs)}, {get_repetitions(baseline_runs)}")
            fig.legend(loc='outside upper right')
        else:
            fig, ax = plt.subplots()
            ax2 = ax.twinx()
            time_name = "time (median)"
            if per_n:
                time_name += " / n"
            ax.plot(num_iters, times, 'b-s', label=time_name)
            if get_repetitions(runs) > 1:
                ax.errorbar(num_iters, errors[0], errors[1], fmt='b.', label='time µ ± σ')
            ax2.plot(num_iters, bits, 'c-o', label='space / n')
            color_y_axis(ax, 'b')
            ax.set_xlabel('n')
            ax.set_xscale('log')
            if per_n:
                ax.set_ylabel('time divided by n in ' + runs[0]['time_unit'])
            else:
                ax.set_ylabel('time in ' + runs[0]['time_unit'])
            ax.set_yscale('log')
            ax.grid()
            color_y_axis(ax2, 'c')
            ax2.set_ylabel('number of bits divided by n')
            ax2.set_yscale('log')
            ax2.set_ylim([None, 1000])
            fig.legend(loc='outside upper right')
            plt.title(name + f"\n#repetitions: {get_repetitions(runs)}")
        yield fig, name


def parse_args():
    parser = ArgumentParser(description="create plots from one or two benchmarking runs")
    parser.add_argument('-i', '--input-file', dest='input_file', default='benchmark_results.json',
                        type=argparse.FileType('r'), nargs=1,
                        help="The input file from which the benchmarking results are read. Must be a JSON file,"
                             " defaults to 'benchmark_results.json'.")
    parser.add_argument('-c', '--compare', dest='compare_file', type=argparse.FileType('r'), nargs=1,
                        help="Baseline to compare the benchmarking results against. Must be a JSON file."
                             "If this option is omitted, the results aren't compared against other results.")
    parser.add_argument('-s', '--show', dest='show_figures', default=False, action="store_true",
                        help="Interactively show the plots in addition to saving them as png images")
    parser.add_argument('-o', '--output-dir', dest='output_dir', default='plots', nargs=1,
                        help="Save plots in this directory.")

    args = parser.parse_args()
    return args


def main():
    baseline_data = None
    args = parse_args()
    plots_folder = args.output_dir
    Path(plots_folder).mkdir(parents=True, exist_ok=True)

    if args.compare_file is not None:
        with open(args.compare_file[0].name, mode='r') as reference_file:
            baseline_data = json.load(reference_file)
    with open(args.input_file.name, mode='r') as json_file:
        data = json.load(json_file)

    for fig, name in generate_plots(data, baseline_data):
        if args.show_figures:
            plt.show()
        fig.savefig(plots_folder + '/' + name + '.png')


if __name__ == '__main__':
    main()
