import argparse
import json
from argparse import ArgumentParser
from dataclasses import dataclass
from itertools import groupby
from os.path import realpath, dirname, commonprefix
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy import isnan

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


def print_avg_variation(benchmarks):
    if is_quick_run(benchmarks):
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


@dataclass
class FamilyInfo:
    median_data: []
    error_bars: np.array
    name: str
    group: float
    quick: bool
    comp_name: str
    has_error_data: bool


def get_family(benchmark_results: [], family_index: int, comp_name: str) -> FamilyInfo:
    group = float('nan')
    if is_quick_run(benchmark_results):
        if family_index >= len(benchmark_results):
            return FamilyInfo(None, None, "<error>", group, True, comp_name, False)
        res = [bm for bm in benchmark_results[family_index] if bm['run_type'] == 'iteration']
        assert len(res) > 0
        if "group" in res[0]:
            group = res[0]["group"]
        return FamilyInfo(res, np.array([[bm['cpu_time'] for bm in res], [0] * len(res)]), get_name(res), group, True,
                          comp_name, False)
    else:
        if family_index not in benchmark_results:
            return FamilyInfo(None, None, "<error>", group, False, comp_name)
        info = benchmark_results[family_index]
        median_data = [bm for bm in info if bm['aggregate_name'] == 'median']
        max_variation = max([float(bm['cpu_time']) for bm in info if bm['aggregate_name'] == 'cv'])
        if max_variation > 0.05:
            print(f"Warning: High maximum measurement variation of {max_variation * 100:.1f} percent for " + get_name(
                median_data))
        stddevs = np.array([bm['cpu_time'] for bm in info if bm['aggregate_name'] == 'stddev'])
        means = np.array([bm['cpu_time'] for bm in info if bm['aggregate_name'] == 'mean'])
        assert len(median_data) == len(stddevs) == len(means) > 0
        error_bars = np.array([means, stddevs])
        if "group" in median_data[0]:
            group = median_data[0]["group"]
        return FamilyInfo(median_data, error_bars, get_name(median_data), group, False, comp_name, max_variation > 0.05)


def get_name(runs):
    if runs is None:
        return "<None>"
    name = runs[0]['run_name'].split('/')[0].split('<')[0]
    if name[:3] == "BM_":
        name = name[3:]
    return name


def get_repetitions(data):
    return int(data[0]['repetitions'])


@dataclass
class FamilyResults:
    times: np.array
    bits: np.array
    input_sizes: np.array
    errors: np.array
    per_n: bool
    subtract_n_from_bits: bool
    repetitions: int
    name: str
    comp_name: str
    has_errors: bool


def get_performance_measurements(family_info: FamilyInfo) -> FamilyResults:
    data = family_info.median_data
    assert len(data) > 1
    per_n_val = int(data[0]['perN']) if 'perN' in data[0] else 0
    subtract_n_from_bits = bool(per_n_val & 2)
    divide_by_n = bool(per_n_val & 1)
    input_sizes = [int(benchmark['run_name'].split('/')[-1]) for benchmark in data]
    times = np.array([benchmark['cpu_time'] for benchmark in data])
    errors = family_info.error_bars
    bits = np.array([benchmark['bits'] for benchmark in data])
    if subtract_n_from_bits:
        bits -= input_sizes
    bits /= input_sizes
    assert (times.shape == bits.shape == errors[0].shape == errors[1].shape == (len(data),))
    return FamilyResults(times, bits, input_sizes, errors, divide_by_n, subtract_n_from_bits, get_repetitions(data),
                         get_name(data), family_info.comp_name, family_info.has_error_data)


def generate_single_plot(family_infos: [FamilyInfo]):
    family_infos = [f for f in family_infos if f is not None]
    assert len(family_infos) > 0
    results: [FamilyResults] = [get_performance_measurements(family) for family in family_infos]
    styles = ['b-s', 'g-o', 'r-*', 'y-3', 'm-x', 'c-d', 'k-+', 'tab:orange', 'tab:purple', 'tab:gray',
              'tab:brown', 'tab:olive', 'tab:pink']
    per_n = all([r.per_n for r in results])
    times = [r.times for r in results]
    has_errors = [r.has_errors for r in results]
    errors = [np.array(r.errors) for r in results]
    bits = [r.bits for r in results]
    if per_n:
        for res, t, e, in zip(results, times, errors):
            t /= res.input_sizes
            e /= res.input_sizes  # median, mean and stddev scale with constant factors
    time_unit = family_infos[0].median_data[0]['time_unit']
    if len(results) > 1:
        # compare several results against each other in two separate graphs for time and space
        for r in results:
            assert len(r.input_sizes) >= len(results[0].input_sizes)
        for i in family_infos:
            assert i.median_data[0]['time_unit'] == time_unit
        fig, ax = plt.subplots(ncols=2, figsize=(10, 5), layout='constrained')
        for res, t, e, show_e, style in zip(results, times, errors, has_errors, styles):
            ax[0].plot(res.input_sizes, t, style, label=res.comp_name + ' (median)')
            if show_e:
                ax[0].errorbar(res.input_sizes, e[0], e[1], fmt=style[0] + '.',
                               label=res.comp_name + ' µ ± σ')
            ax[1].plot(res.input_sizes, res.bits, style)

        ax[0].set_xlabel('n')
        ax[0].set_xscale('log')
        ax[1].set_xlabel('n')
        ax[1].set_xscale('log')
        if results[0].subtract_n_from_bits:
            ax[1].set_ylabel('number of additional bits divided by n')
        else:
            ax[1].set_ylabel('number of bits divided by n')
        ax[1].set_yscale('log')
        ax[1].grid(which='major', alpha=0.7)
        ax[1].grid(which='minor', alpha=0.2)
        if results[0].subtract_n_from_bits:
            ax[1].set_title("additional space (divided by n)")
        else:
            ax[1].set_title("space (divided by n)")
        if max([max(res.bits) for res in results]) > 1000:
            ax[1].set_ylim([None, 1000])
        if per_n:
            ax[0].set_ylabel('time in ' + time_unit + ' divided by n')
        else:
            ax[0].set_ylabel('time in ' + time_unit)
        ax[0].set_yscale('log')
        ax[0].grid(which='major', alpha=0.7)
        ax[0].grid(which='minor', alpha=0.2)
        title = "time"
        if per_n:
            title += " (divided by n)"
        ax[0].set_title(title)
        if len(family_infos) == 2:
            rel_time = ax[0].twinx()
            n0 = results[0].input_sizes
            n1 = results[1].input_sizes
            while n0[0] < n1[0]:
                times[0] = times[0][1:]
                n0 = n0[1:]
                bits[0] = bits[0][1:]
            while n0[0] > n1[0]:
                times[1] = times[1][1:]
                n1 = n1[1:]
                bits[1] = bits[1][1:]
            times[0] = times[0][:min(len(times[0]), len(times[1]))]
            if len(n0) > 1 and n0[1] == n1[1]:
                time_percent = [100 * times[0][i] / times[1][i] for i in range(len(times[0]))]
                color_y_axis(rel_time, 'c')
                rel_time.set_ylabel('relative time in percent')
                rel_time.plot(n0, time_percent, 'c-p', label='percent of baseline')
                rel_time.axhline(y=100, color='c', linestyle='--')

                rel_space = ax[1].twinx()
                time_percent = [100 * bits[0][i] / bits[1][i] for i in range(len(bits[0]))]
                rel_space.plot(n0, time_percent, 'c-p')
                rel_space.axhline(y=100, color='c', linestyle='--')
                color_y_axis(rel_space, 'c')
                rel_space.set_ylabel('relative space in percent')
        names = list(set(r.name for r in results))
        prefix = commonprefix(names)
        name = prefix + "\n" + names[0][len(prefix):]
        for n in names[1:]:
            name += ' vs ' + n[len(prefix):]
        repetitions = str(results[0].repetitions)
        for i in range(1, len(results)):
            repetitions += ', ' + str(results[i].repetitions)
        plt.suptitle(name + f"\n#repetitions: {repetitions}")
        fig.legend(loc='outside upper right')
    else:  # len results == 1
        res = results[0]
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        time_name = "time (median)"
        if per_n:
            time_name += " / n"
        ax.plot(res.input_sizes, res.times, 'b-s', label=time_name)
        if res.repetitions > 1:
            ax.errorbar(res.input_sizes, errors[0][0], errors[0][1], fmt="b.", label='time µ ± σ')
        ax2.plot(res.input_sizes, res.bits, styles[1], label='space')
        color_y_axis(ax, styles[0][0])
        ax.set_xlabel('n')
        ax.set_xscale('log')
        if per_n:
            ax.set_ylabel('time divided by n in ' + time_unit)
        else:
            ax.set_ylabel('time in ' + time_unit)
        ax.set_yscale('log')
        ax.grid(which='major', alpha=0.7)
        ax.grid(which='minor', alpha=0.2)
        color_y_axis(ax2, styles[1][0])
        if res.subtract_n_from_bits:
            ax2.set_ylabel('number of additional bits divided by n')
        else:
            ax2.set_ylabel('number of bits divided by n')
        ax2.set_yscale('log')
        if max(res.bits) > 1000:
            ax2.set_ylim([None, 1000])
        fig.legend(loc='outside upper right')
        plt.title(res.name + f"\n#repetitions: {res.repetitions}")
    return fig


def generate_plots(data, baseline_data):
    baseline_benchmarks = None
    current_benchmark_results = get_benchmarks(data)
    if not current_benchmark_results:
        raise ValueError("No benchmarks were run")
    num_families = len(current_benchmark_results)
    print_avg_variation(current_benchmark_results)
    if baseline_data is not None:
        baseline_benchmarks = get_benchmarks(baseline_data)
        check_same_config(data, current_benchmark_results, baseline_data, baseline_benchmarks)

    group_results = {}
    baseline_families = {}
    if baseline_data is not None:
        for i in range(len(baseline_benchmarks)):
            tmp = get_family(baseline_benchmarks, i, 'baseline')
            baseline_families[tmp.name] = tmp
    for i in range(0, num_families):
        current_family = get_family(current_benchmark_results, i, 'current')
        baseline_family = None
        name = current_family.name
        if not isnan(current_family.group):
            if current_family.group not in group_results:
                group_results[current_family.group] = []
            if -int(current_family.group) not in group_results:
                group_results[-int(current_family.group)] = []
            group_results[current_family.group].append(current_family)
            if current_family.group != int(current_family.group):
                group_results[-int(current_family.group)].append(current_family)
        if baseline_data is not None:
            if name not in baseline_families:
                print("Warning: No matching baseline results found for " + name)
            else:
                baseline_family = baseline_families[name]
        fig = generate_single_plot([current_family, baseline_family])
        yield fig, name
    for group in group_results.values():
        if len(group) <= 1:
            continue
        for family in group:
            family.comp_name = family.name
        fig = generate_single_plot(group)
        prefix = commonprefix([g.name for g in group])
        name = "compare_" + str(len(group)) + '_' + prefix + 's' + '_' + group[0].name[len(prefix):]
        yield fig, name


def parse_args():
    parser = ArgumentParser(description="create plots from one or two benchmarking runs")
    parser.add_argument('-i', '--input-file', dest='input_file', default='results/benchmark_results.json',
                        type=argparse.FileType('r'),
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
