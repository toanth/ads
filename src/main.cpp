
#include "../include/elias_fano.hpp"
#include "../include/linear_space_rmq.hpp"
#include "../include/succinct_rmq.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace ads;

using Rmq = LinearSpaceRMQ<U64>;

struct PredecessorInput {
    std::vector<U64> sequence;
    std::vector<U64> queries;
};

struct RmqInput {
    std::vector<U64> sequence;
    std::vector<std::pair<U64, U64>> queries;
};


PredecessorInput readPredecessorInput(std::string_view filename) {
    PredecessorInput result;

    std::ifstream file(filename.data());
    std::ifstream::sync_with_stdio(false);
    Index n;
    file >> n;
    result.sequence.resize(n);
    for (Index i = 0; i < n; ++i) {
        file >> result.sequence[i];
    }
    U64 query;
    while (file >> query) {
        result.queries.push_back(query);
    }
    return result;
}

RmqInput readRmqInput(std::string_view filename) {
    RmqInput input;


    std::ifstream file(filename.data());
    std::ifstream::sync_with_stdio(false);
    Index n;
    file >> n;
    input.sequence.resize(n);
    for (Index i = 0; i < n; ++i) {
        file >> input.sequence[i];
    }
    U64 query_min, query_max;
    char comma;
    while (file >> query_min >> comma >> query_max) {
        ++query_max; // use [min, max) instead of [min, max]
        input.queries.emplace_back(query_min, query_max);
    }
    return input;
}

void writePredecessorOutput(std::string_view filename, Span<U64> output) {
    std::ofstream file(filename.data());
    std::ofstream::sync_with_stdio(false);
    for (U64 val : output) {
        file << val << '\n';
    }
}

void writeRmqOutput(std::string_view filename, Span<std::pair<U64, U64>> output) {
    std::ofstream file(filename.data());
    std::ofstream::sync_with_stdio(false);
    for (auto val : output) {
        file << val.first << '\n';
    }
}

/// \brief Builds the Elias-Fano datastructure from `input.sequence` and uses it to compute the predecessors to
/// `input.queries`.
/// \param input A list of Integers that will make up the sequence stored in Elias-Fano and a list of
/// integers for which the predecessors will be computed, which will be overwritten with the answers to those queries.
/// \return The total space in bits needed by the Elias fano datastructure (on the heap and on the stack), without
/// taking into account the space needed to store the input. The answers will be stored in the input.queries vector.
Index eliasFano(PredecessorInput& input) {
    EliasFano<> ef(input.sequence);
    // TODO: Overwrite queries vector instead of creating new vector
    for (Index i = 0; i < input.queries.size(); ++i) {
        input.queries[i] = ef.predecessor(input.queries[i]);
    }
    std::cout << ef.getUpper().allocatedSizeInBits() << " " << ef.numUpperBitsPerNumber() << " "
              << ef.numLowerBitsPerNumber() << std::endl;
    return ef.numAllocatedBits();
}

Index rmq(RmqInput& input) {
    Rmq rmq(input.sequence);
    for (auto& query : input.queries) {
        query.first = rmq(query.first, query.second);
    }
    return rmq.allocatedSizeInBits();
}


void printUsage() {
    std::cout << "Usage: `<program> [pd|rmq] input_file output_file";
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printUsage();
        return 1;
    }
    Index spaceInBits = -1;
    Index time = -1;
    if (argv[1] == std::string("pd")) {
        PredecessorInput input = readPredecessorInput(argv[2]);
        auto before = std::chrono::system_clock::now();
        spaceInBits = eliasFano(input);
        auto after = std::chrono::system_clock::now();
        time = std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
        writePredecessorOutput(argv[3], input.queries);

    } else if (argv[1] == std::string("rmq")) {
        RmqInput input = readRmqInput(argv[2]);
        auto before = std::chrono::system_clock::now();
        spaceInBits = rmq(input);
        auto after = std::chrono::system_clock::now();
        time = std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
        writeRmqOutput(argv[3], input.queries);
    } else {
        printUsage();
        return 2;
    }
    std::cout << "RESULT algo=rmq name=tobias_theuer time=" << time << " space=" << spaceInBits << std::endl;
    return 0;
}
