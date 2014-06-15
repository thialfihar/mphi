#include <cstdio>
#include <cstring>
#include <boost/preprocessor/repetition/repeat.hpp>

#include <fstream>

#include "Matrix.h"

#define N 3
#define K 3

int main(int argc, char *argv[]) {
    FILE *json;
    char buf[100];
    ExponentMapResult result = Matrix<N, K>::get_exponent_map();
    mpz_class sample_size = (result.confidence == 0 ? result.phi : result.confidence);
    sprintf(buf, "map_%d_%d.json", N, K);
    json = fopen(buf, "w");
    printf("phi: %s\ncon: %d\n", result.phi.get_str().c_str(), result.confidence);
    fprintf(json, "{\"n\": %d, \"k\": %d, \"phi\": %s, \"confidence\": %d, \"map\": {\n",
        N, K, result.phi.get_str().c_str(), result.confidence);
    unsigned int i = 0;
    for (const auto &pair : result.map) {
        mpz_class percentage = 10000 * pair.second / sample_size;
        printf("%s %s %.2f%% %s\n",
            pair.first.get_str().c_str(), pair.second.get_str().c_str(),
            percentage.get_ui() / 100.0, result.permutation_map[pair.first].get_str().c_str());
        fprintf(json, "\"%s\": %s", pair.first.get_str().c_str(), pair.second.get_str().c_str());
        ++i;
        if (i < result.map.size()) {
            fprintf(json, ",\n");
        } else {
            fprintf(json, "\n");
        }
    }
    fprintf(json, "}");

    fprintf(json, ",\n\"cycles\": [\n");
    for (unsigned int i = 0; i < result.cycles.size(); ++i) {
        const auto &cycle = result.cycles[i];
        fprintf(json, "[");
        for (unsigned int j = 0; j < cycle.size(); ++j) {
            fprintf(json, "%s", cycle[j].get_str().c_str());
            if (j != cycle.size() - 1) {
                fprintf(json, ", ");
            }
        }
        fprintf(json, "]");
        if (i != result.cycles.size() - 1) {
            fprintf(json, ",\n");
        }
    }
    fprintf(json, "]");

    fprintf(json, "}\n");
    fclose(json);
    return 0;
}
