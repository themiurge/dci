#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../core/dci.hpp"
#include <fstream>

TEST_CASE( "Empty system creation", "[system]") {
    dci::system d;
    REQUIRE( d.num_samples() == 0 );
    REQUIRE( d.num_agents() == 0 );
}

TEST_CASE( "Default system creation", "[system]") {
    dci::system d(10, 100);
    REQUIRE( d.num_samples() == 100 );
    REQUIRE( d.num_agents() == 10 );
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            REQUIRE( d[s][a] == 0U );
}

TEST_CASE( "System fill", "[system]") {
    dci::system d(10, 100);
    REQUIRE( d.num_samples() == 100 );
    REQUIRE( d.num_agents() == 10 );
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = (a+s)%2;
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            REQUIRE( d[s][a] == (a+s)%2 );

}

TEST_CASE( "Default system creation - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    REQUIRE( d.num_samples() == 100 );
    REQUIRE( d.num_agents() == 10 );
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            REQUIRE( d[s][a] == 0U );
}

TEST_CASE( "System fill - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = (a+s)%16;
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            REQUIRE( d[s][a] == (a+s)%16 );
}

TEST_CASE( "System equality", "[system]") {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = 1U;
    dci::system e(10, 100);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            e[s][a] = 1U;
    REQUIRE( d == d );
    REQUIRE( !(d != d) );
    REQUIRE( d == e );
    REQUIRE( !(d != e) );
    REQUIRE( e == d );
    REQUIRE( !(e != d) );
    e[0][0] = 0U;
    REQUIRE( d != e );
    REQUIRE( e != d );
}

TEST_CASE( "System equality - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = (a+s)%16;
    dci::system e(10, 100, 4);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            e[s][a] = (a+s)%16;
    REQUIRE( d == d );
    REQUIRE( !(d != d) );
    REQUIRE( d == e );
    REQUIRE( !(d != e) );
    REQUIRE( e == d );
    REQUIRE( !(e != d) );
    e[0][0] = 1U;
    REQUIRE( d != e );
    REQUIRE( e != d );
}

TEST_CASE( "System save/load", "[system]" ) {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = (a+s)%2;
    ofstream fout("test.txt");
    fout << d;
    fout.close();
    dci::system e;
    ifstream fin("test.txt");
    fin >> e;
    fin.close();
    REQUIRE( d == e );
}

TEST_CASE( "System save/load - multilevel", "[system]" ) {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.num_samples(); ++s)
        for (size_t a = 0; a != d.num_agents(); ++a)
            d[s][a] = (a+s)%16;
    ofstream fout("test.txt");
    fout << d;
    fout.close();
    dci::system e;
    ifstream fin("test.txt");
    fin >> e;
    fin.close();
    REQUIRE( d == e );
}
