import json
from os.path import realpath, dirname
from pathlib import Path

import matplotlib.pyplot as plt


def color_y_axis(ax, color):
    for t in ax.get_yticklabels():
        t.set_color(color)


def get_runs(all_benchmarks, family_index, is_quick_run):
    if is_quick_run:
        return [bm for bm in all_benchmarks if bm['family_index'] == family_index and bm['run_type'] == 'iteration']
    else:
        max_variation = max([float(bm['cpu_time']) for bm in all_benchmarks if
                             bm['family_index'] == family_index and bm['aggregate_name'] == 'cv'])
        if max_variation > 0.02:
            print("Warning: high maximum measurement variation of " + str(max_variation * 100) + " percent")
        return [bm for bm in all_benchmarks if bm['family_index'] == family_index and bm['aggregate_name'] == 'median']


def generate_plots(data):
    all_benchmarks = data['benchmarks']
    if not all_benchmarks:
        raise ValueError("No benchmarks were run")
    num_families = all_benchmarks[-1]['family_index']
    is_quick_run = all_benchmarks[0]['name'] == all_benchmarks[0]['run_name']
    for i in range(0, num_families + 1):
        runs = get_runs(all_benchmarks, i, is_quick_run)
        name = runs[0]['run_name'].split('/')[0].split('<')[0]
        if name[:3] == "BM_":
            name = name[3:]
        assert len(runs) > 1
        times = [benchmark['cpu_time'] for benchmark in runs]
        bits = [benchmark['bits'] for benchmark in runs]
        num_iters = [int(benchmark['run_name'].split('/')[-1]) for benchmark in runs]
        print(num_iters)
        print(times)
        fig, ax = plt.subplots()
        ax.plot(num_iters, times, 'b')
        color_y_axis(ax, 'b')
        ax.set_xlabel('n')
        ax.set_xscale('log')
        ax.set_ylabel('time in ' + runs[0]['time_unit'])
        ax.set_yscale('log')
        plt.grid()
        ax2 = ax.twinx()
        ax2.plot(num_iters, bits, 'c')
        color_y_axis(ax2, 'c')
        ax2.set_ylabel('number of bits')
        ax2.set_yscale('log')
        plt.title(name)
        yield fig, name


def main():
    basepath = dirname(realpath(__file__)) + "/"
    with open(basepath + "benchmark_results.json", mode='r') as json_data:
        data = json.load(json_data)
    plots_folder = basepath + "plots"
    Path(plots_folder).mkdir(parents=True, exist_ok=True)
    for fig, name in generate_plots(data):
        fig.savefig(plots_folder + '/' + name + '.png')


if __name__ == '__main__':
    main()
