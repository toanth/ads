#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

usage() {
    cat << EOF >&2
Usage: $0 [options] | $0 [options] -- [binary options]
Run benchmarks with the given options.
Options include:

-o | --output-file: The name of the JSON output file. Defaults to benchmark_results.json.
-b | --binary:      Path to the binary to execute. Defaults to "../build/Release/benchmarks/benchmarks".
-c | --change-environment:  Temporarily change system settings to decrease measurement
                    variance similarly to ./setupTestEnv, which requires root privileges. To be used with caution.
                    If the script is aborted (such as by pressing ctrl c), the system may be left in the benchmarking state,
                    which can be left by running ./setupTestEnv --unset
-q | --quick:       Only repeat each test once (the default behavior of google benchmark; otherwise 3 repetitions are used)
--help:             Print this message.
EOF
}

if [[ $1 == "--help" ]]; then
    usage
    exit 0
fi


# getopt usage inspired by
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash/29754866#29754866
getopt --test
if [[ $? -ne 4 ]]; then
    echo "getopt didn't work"
fi

LONGOPTS=output-file:,binary:,change-environment,quick
OPTIONS=o:,b:,c,q

PARSED="$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")"
if [[ $? -ne 0 ]]; then
    echo "couldn't parse options"
    usage
    exit 2
fi

eval set -- "$PARSED"

sudoChanges=0
quick=0
outputFile=benchmark_results.json
binPath="../build/Release/benchmarks/benchmarks"

moreThan4Cores=0 # if there are more than 4 logical CPUs, reserve 2 for benchmarking and turn off their SMT siblings
if [[ $(nproc) -gt 4 ]]; then
    moreThan4Cores=1 # set that now because cset shield reduces the number of available cores
fi

while true; do
    case "$1" in
        --output-file|-o)
          outputFile="$2"
          shift 2
          ;;
        --binary|-b)
          binPath="$2"
          shift 2
          ;;
        --quick|-q)
            quick=1
            shift
            ;;
        --change-environment|-c)
            sudoChanges=1
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unrecognized option"
            usage
            exit 3
            ;;
    esac
done

oldGovernor="schedutil"

runBenchmark() {
    if [[ $sudoChanges == 1 ]]; then
        user=$(whoami)
        oldGovernor=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
        echo "cpufreq performance governor was ${oldGovernor} and will now be set to 'performance'"
        sudo cpupower frequency-set -g performance
        if [[ $moreThan4Cores == 1 ]]; then
            # reserve 2 CPUs so that we can also run perf, which should run on a different CPU
            # `cset shield -c 0,2 -k on` reserves the logical cpus 0 and 1 and moves all threads, including kernel threads, out of them
            sudo cset shield -c 0,1,2,3 -k on
            echo "allowing user ${user} to run in the 'user' cpuset with cset exec without root privileges"
            sudo chown -R ${user} /cpusets/user # the benchmark shouldn't be executed with root privileges
            # typically, cpus 0 and 1 form an SMT pair, as well as cpus 2 and 3
            echo "disabling CPUs 1 and 3, current status: $(sudo cat /sys/devices/system/cpu/cpu1/online) $(sudo cat /sys/devices/system/cpu/cpu3/online)"
            echo 0 | sudo tee /sys/devices/system/cpu/cpu1/online
            echo 0 | sudo tee /sys/devices/system/cpu/cpu3/online
        fi

        # `cset shield --exec --`` runs the given program in the user cpuset, on which no other tasks are being run
        # `setarch $(uname -m) -R` disables ASLR to reduce variance (but obviously weakens security)
        cset shield --exec -- setarch "$(uname -m)" -R "${binPath}" ${benchmarkOpts}

        if [[ $moreThan4Cores == 1 ]]; then
            echo "Enabling cpus 1 and 3"
            echo 1 | sudo tee /sys/devices/system/cpu/cpu1/online
            echo 1 | sudo tee /sys/devices/system/cpu/cpu3/online
            echo "resetting cset shield"
            sudo cset shield --reset
        fi
        echo "Resetting performance governor to ${oldGovernor}"
        sudo cpupower frequency-set -g ${oldGovernor}
    else
        if [[ $(cset shield | grep '"user" cpuset of CPUSPEC(0,2) with 0 tasks running' -c) == 1 ]]; then
            echo "Executing benchmarks in the benchmarking system state"
            # in benchmarking state (user ran ./setupTestEnv --set)
            # `cset shield --exec --`` runs the given program in the user cpuset, on which no other tasks are being run
            # `setarch $(uname -m) -R` disables ASLR to reduce variance (but obviously weakens security)
            # because the user should have been set as the owner of the user cpuset, this should work without root privileges
            cset shield --exec -- setarch "$(uname -m)" -R "${binPath}" ${benchmarkOpts}
        else
            echo "Executing benchmarks in the normal system state"
            # `taskset -c 0` sets the task affinity to cpu 0, but doesn't prevent the scheduler from scheduling other programs there
            taskset -c 0 setarch "$(uname -m)" -R "${binPath}" ${benchmarkOpts}
        fi
    fi
}

export CLICOLOR_FORCE=1
benchmarkOpts="$@ $binFlags --benchmark_out=$outputFile --benchmark_out_format=json --benchmark_counters_tabular=true"
if [[ $quick == 0 ]]; then
    benchmarkOpts="$benchmarkOpts --benchmark_repetitions=3 --benchmark_report_aggregates_only=true"
fi

echo "executing command ${binPath} $benchmarkOpts"

(cd "${SCRIPTPATH}" &&
time runBenchmark "$@")
