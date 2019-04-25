#include <nlohmann/json.hpp>
#include <NTL/GF2X.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>

using namespace std;
using namespace NTL;

using json = nlohmann::json;

string to_string(vec_GF2 input) {
    stringstream buffer;
    buffer << input;
    return buffer.str();
}

string to_string(GF2X input) {
    stringstream buffer;
    buffer << input;
    return buffer.str();
}

/**
 * Calculates minimal polynomial and linear complexity
 * of linear reccurent sequence in GF2
 */
pair<int, GF2X> berlekamp_massey(const vec_GF2 &sequence) {
    long N = sequence.length();

    /** 
     * ? Step 2:
     * ? Initialise two arrays b and c of Length N to be zeros,
     * ? except b_0 and c_0, which should be 1
     */
    vec_GF2 b = vec_GF2(INIT_SIZE, N);
    vec_GF2 c = vec_GF2(INIT_SIZE, N);
    b.at(0) = GF2(1);
    c.at(0) = GF2(1);

    // ? Step 3
    int L = 0;
    int m = -1;

    // ? Step 4
    for(int n = 0; n < N; n++) {
        // ? Step 4.1 
        GF2 d = GF2(0);
        for(int i = 0; i <= L; i++) d += c.at(i) * sequence.at(n - i);

        // ! Step 4.2 skipped, because we are in GF2

        // ? Step 4.3
        if(d == GF2(1)) {
            // ? Step 4.3.1
            vec_GF2 t = vec_GF2(c);

            // ? Step 4.3.2
            int n_m = n - m;
            for(int j = 0; j < N - n_m; j++) c.at(n_m + j) += b.at(j);

            // ? Step 4.3.3
            if(L <= n / 2) {
                L = n + 1 - L;
                m = n;
                b = vec_GF2(t);
            }
        }
    }

    // ? Conversion of GF2 vector to polynomial
    GF2X min_polynomial = reverse(conv<GF2X>(c));
    if(L > deg(min_polynomial)) min_polynomial = LeftShift(min_polynomial, L - deg(min_polynomial));

    return pair<int, GF2X>(L, min_polynomial);
}

int main(int argc,  char** argv) {
    if(argc != 2) return 1;

    json output;
    vec_GF2 lin_seq;
    vector<long> vec_lin_comp;
    string input = argv[1];


    for(char c : input) {
        long coeff = c - '0';
        lin_seq.append(GF2(coeff));
    }

    auto result = berlekamp_massey(lin_seq);
    output["sequences"].push_back({
        {"sequence", to_string(lin_seq)},
        {"min_polynomial", to_string(result.second)},
        {"linear_complexity", result.first}
    });

    for(int i = 0; i < lin_seq.length(); ++i) {
        vec_GF2 temp_seq = lin_seq;
        temp_seq.at(i) += GF2(1);

        result = berlekamp_massey(temp_seq);
        vec_lin_comp.push_back(result.first);

        output["sequences"].push_back({
            {"sequence", to_string(temp_seq)},
            {"min_polynomial", to_string(result.second)},
            {"linear_complexity", result.first}
        });
    }
    output["spherical_complexity"] = *min_element(vec_lin_comp.begin(), vec_lin_comp.end());
    cout << output.dump(4) << endl;
    
    return 0;
}