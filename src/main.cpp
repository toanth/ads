#include "../include/bitvector.hpp"
#include "../include/elias_fano.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace ads;

struct Input {
    std::vector<Elem> sequence;
    std::vector<Elem> queries;
};


Input readInput(std::string_view filename) {
    Input result;

    std::ifstream file(filename.data());
    std::ifstream::sync_with_stdio(false);
    Index n;
    file >> n;
    result.sequence.resize(n);
    for (Index i = 0; i < n; ++i) {
        file >> result.sequence[i];
    }
    Elem query;
    while (file >> query) {
        result.queries.push_back(query);
    }
    return result;
}

void writeOutput(std::string_view filename, Span<Elem> output) {
    std::ofstream file(filename.data());
    std::ofstream::sync_with_stdio(false);
    for (Elem val : output) {
        file << val << '\n';
    }
}

/// \brief Builds the Elias-Fano datastructure from `input.sequence` and uses it to compute the predecessors to
/// `input.queries`. \param input A list of Integers that will make up the sequence stored in Elias-Fano and a list of
/// integers for which the predecessors will be computed. \param spaceInBits An out param that stores the total space in
/// bits needed by the Elias fano datastructure (on the heap and on the stack), without taking into account the space
/// needed to store the input. \return The results of the predecessor queries.
std::vector<Elem> eliasFano(Input input, Index& spaceInBits) {
    EliasFano<> ef(input.sequence);
    std::vector<Elem> result(input.queries.size());
    // TODO: Overwrite queries vector instead of creating new vector
    for (Index i = 0; i < result.size(); ++i) {
        result[i] = ef.predecessor(input.queries[i]);
    }
    spaceInBits = ef.numAllocatedBits();
    return result;
}


void printUsage() {
    std::cout << "Usage: `<program> [pd|rmq] input_file output_file";
}

int main(int argc, char** argv) {
    // TODO: Modulo and division by powers of two is less efficient for signed integer types than for unsigned
    if (argc != 4) {
        printUsage();
        return 1;
    }
    if (argv[1] == std::string("pd")) {
        Input input = readInput(argv[2]);
        Index spaceInBits = 0;
        auto before = std::chrono::system_clock::now();
        std::vector<Elem> results = eliasFano(std::move(input), spaceInBits);
        auto after = std::chrono::system_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
        writeOutput(argv[3], results);
        std::cout << "RESULT algo=md nametobias_theuer time=" << time.count() << " space=" << spaceInBits << std::endl;

    } else if (argv[1] == std::string("rmq")) {
        throw "not implemented!";
    } else {
        printUsage();
        return 2;
    }
    return 0;
}
