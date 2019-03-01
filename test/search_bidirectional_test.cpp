#include "common.hpp"
#include <sdsl/suffix_arrays.hpp>
#include "gtest/gtest.h"
#include <iostream>
#include <vector>
#include <string>

namespace
{

using namespace sdsl;
using namespace std;

typedef sdsl::int_vector<>::size_type size_type;

tMSS   test_case_file_map;
tMSS   test_case_file_map_rev;
string test_file;
string test_file_rev;
string temp_file;
string temp_file_rev;
string temp_dir;

template<class T>
class search_bidirectional_test : public ::testing::Test { };

using testing::Types;

typedef Types<
    csa_wt<wt_blcd<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet>,
    csa_wt<wt_blcd<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> >,
    csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet>,
    csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> >,
    csa_wt<wt_hutu<bit_vector_il<> >, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet>
> Implementations;

TYPED_TEST_CASE(search_bidirectional_test, Implementations);

//! Compare bidirectional search and backward search
TYPED_TEST(search_bidirectional_test, bidirectional_search)
{
    bool debug = false;

    TypeParam csa1;
    TypeParam csa1_rev;

    temp_file = sdsl::tmp_file(temp_dir+"/"+util::basename(test_file),util::basename(test_file));
    cache_config config(false, temp_dir, util::basename(temp_file));
    temp_file_rev = sdsl::tmp_file(temp_dir+"/"+util::basename(test_file_rev),util::basename(test_file_rev));
    cache_config config_rev(false, temp_dir, util::basename(temp_file_rev));

    construct(csa1, test_file, config, 1);
    construct(csa1_rev, test_file_rev, config_rev, 1);

    test_case_file_map = config.file_map;
    test_case_file_map_rev = config_rev.file_map;

    std::mt19937_64 rng(13);
    std::uniform_int_distribution<uint64_t> distribution(0, csa1.size()-1);

    for (size_type h = 0; h<1000; ++h) {
        //search for an existing pattern forward and backward using bidirectional_search:
        size_type x = 4; // number of characters that are added to the pattern in each step
        size_type steps = 10; // maximal number of patternparts that are searched for
        size_type start = distribution(rng); //inclusive
        size_type end = start;  //exclusive
        bool forward = false;

        size_type l_rev = 0;
        size_type r_rev = csa1_rev.size()-1;
        size_type l = 0;
        size_type r = csa1.size()-1;
        size_type occ = csa1.size(); // size of the interval equals the number of occurrences of the pattern in the text (init empty pattern)
        size_type i,pos;

        // alternating forward and backward search using bidirectional_search: alternating every x characters
        for (size_type rep = 0; rep<steps and start>=x and end<=csa1.size()-x ; ++rep) {
            string newpat = "";
            if (forward) {
                //forward
                i = 0;
                pos = end;
                for (size_type j=0; j<x; ++j) {
                    newpat.push_back(csa1.text[pos+j]);
                }
                occ = bidirectional_search_forward(csa1, csa1_rev, l, r, l_rev, r_rev, newpat.begin(), newpat.end(), l, r, l_rev, r_rev);
                i = newpat.size();
                end += i;
            } else {
                //backward
                i = 0;
                pos = start-1;
                for (size_type j=0; j<x; ++j) {
                    newpat.push_back(csa1.text[pos-x+1+j]);
                }
                occ = bidirectional_search_backward(csa1, csa1_rev, l, r, l_rev, r_rev, newpat.begin(), newpat.end(), l, r, l_rev, r_rev);
                i = newpat.size();
                start -= i;
            }

            //output
            if (debug) {
                cout << "pattern (at text[" << start << ".." << end-1 << "]):" << endl;
                for (size_type j=start; j<end; ++j) {
                    cout << csa1.text[j];
                }
                cout << endl;
                if (occ) {
                    cout << "interval of pattern in csa1 is [" << l   << ".." << r   << "]" << endl;
                    cout << "interval of reverse pattern in csa1_rev is [" << l_rev << ".." << r_rev << "]" << endl;
                }
                cout << endl;
            }
            ASSERT_LT((size_type)0, occ) << "Pattern not found in input."; // make sure pattern was found in input (it has to be because we took part of the input as pattern)

            {
                //check using backward_search
                string pat = "";
                for (size_type j=0; j<end-start; ++j) {
                    pat.push_back(csa1.text[start+j]);
                }
                size_type b_l,b_r;
                size_type b_occ = backward_search(csa1, 0, csa1.size()-1, pat.begin(), pat.end(), b_l, b_r);
                ASSERT_EQ(b_occ, occ) << "Bidirectional_search and backward_search found different number of occurrences of the pattern.";
                ASSERT_EQ(b_l, l) << "Bidirectional_search and backward_search found different left border";
                ASSERT_EQ(b_r, r) << "Bidirectional_search and backward_search found different right border";
            }

            //change direction
            forward = !forward;
        }
    }
}

TYPED_TEST(search_bidirectional_test, delete_)
{
    sdsl::remove(temp_file);
    sdsl::remove(temp_file_rev);
    util::delete_all_files(test_case_file_map);
    util::delete_all_files(test_case_file_map_rev);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (init_2_arg_test(argc, argv, "SEARCH_BIDIRECTIONAL", test_file, temp_dir, temp_file) != 0 ) {
        return 1;
    }
    test_file_rev = test_file + "_rev";

    {
        //reverse input
        int_vector<8> text;
        load_vector_from_file(text, test_file, 1);
        size_type n = text.size();
        int_vector<8> text_rev(n);
        for (size_type i=0; i<n; i++) {
            text_rev[n-1-i] = text[i];
        }
        char* text2 = (char*)text_rev.data();
        ofstream of(test_file_rev, ofstream::binary);
        of.write(text2, n);
        of.close();
    }
    int result = RUN_ALL_TESTS();
    sdsl::remove(test_file_rev);
    return result;
}
