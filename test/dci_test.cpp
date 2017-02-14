#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../core/dci.hpp"
#include <fstream>

TEST_CASE( "Enter" ) {
    cout << "Entering test suite\n";
}

TEST_CASE( "Empty system creation", "[system]") {
    dci::system d;
    REQUIRE( d.samples().size() == 0 );
    REQUIRE( d.agents().size() == 0 );
}

TEST_CASE( "Default system creation", "[system]") {
    dci::system d(10, 100);
    REQUIRE( d.samples().size() == 100 );
    REQUIRE( d.agents().size() == 10 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            REQUIRE( d[s][a] == 0 );
    for (const auto& a : d.agents())
        REQUIRE ( a->size() == 1 );
}

TEST_CASE( "System fill", "[system]") {
    dci::system d(10, 100);
    REQUIRE( d.samples().size() == 100 );
    REQUIRE( d.agents().size() == 10 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (a+s)%2;
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            REQUIRE( d[s][a] == (a+s)%2 );

}

TEST_CASE( "Default system creation - multilevel", "[system]") {
    constexpr size_t agent_size_in_bits = 200;
    constexpr size_t agent_size_in_regs = dci::register_utils::number_of_registers(agent_size_in_bits);
    constexpr size_t agent_size_in_bytes = agent_size_in_regs * dci::bytes_per_reg;
    dci::reg_type all_zeros[agent_size_in_regs];
    memset(all_zeros, 0, agent_size_in_bytes);

    dci::system d(10, 100, agent_size_in_bits);
    REQUIRE( d.samples().size() == 100 );
    REQUIRE( d.agents().size() == 10 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            REQUIRE( d[s][a] == all_zeros );
    for (const auto& a : d.agents())
        REQUIRE ( a->size() == agent_size_in_bits );
}

TEST_CASE( "System fill - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (a+s)%16;
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            REQUIRE( d[s][a] == (a+s)%16 );
}

TEST_CASE( "System equality", "[system]") {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = 1;
    dci::system e(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            e[s][a] = 1;
    REQUIRE( d == d );
    REQUIRE( !(d != d) );
    REQUIRE( d == e );
    REQUIRE( !(d != e) );
    REQUIRE( e == d );
    REQUIRE( !(e != d) );
    e[0][0] = 0;
    REQUIRE( d != e );
    REQUIRE( e != d );
}

TEST_CASE( "System equality - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (a+s)%16;
    dci::system e(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            e[s][a] = (a+s)%16;
    REQUIRE( d == d );
    REQUIRE( !(d != d) );
    REQUIRE( d == e );
    REQUIRE( !(d != e) );
    REQUIRE( e == d );
    REQUIRE( !(e != d) );
    e[0][0] = 1;
    REQUIRE( d != e );
    REQUIRE( e != d );
}

TEST_CASE( "System save/load", "[system]" ) {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (a+s)%2;
    d.save_to_file("test.txt");
    dci::system e = dci::system::load_from_file("test.txt");
    REQUIRE( d == e );
}

TEST_CASE( "System save/load - multilevel", "[system]" ) {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (a+s)%16;
    d.save_to_file("test.txt");
    dci::system e = dci::system::load_from_file("test.txt");
    REQUIRE( d == e );
}

TEST_CASE( "Copy constructor and assignment", "[system]") {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = 1;
    dci::system e(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            e[s][a] = 0;
    dci::system f = d;
    REQUIRE( d == f );
    f = e;
    REQUIRE( f == e );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            f[s][a] = 1;
    REQUIRE( f == d );
    REQUIRE( f != e );
}

TEST_CASE( "Copy constructor and assignment - multilevel", "[system]") {
    dci::system d(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = 15;
    dci::system e(10, 100, 4);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            e[s][a] = 0;
    dci::system f(d);
    REQUIRE( d == f );
    f = e;
    REQUIRE( f == e );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            f[s][a] = 15;
    REQUIRE( f == d );
    REQUIRE( f != e );
}

TEST_CASE( "Move constructor and assignment", "[system]") {
    dci::system empty;
    dci::system d;
    REQUIRE( d == empty );
    d = dci::system(10, 100);
    REQUIRE( d.agents().size() == 10 );
    REQUIRE( d.samples().size() == 100 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = 1;
    
    dci::system e = d;
    
    dci::system f = move(d);
    for (size_t s = 0; s != e.samples().size(); ++s)
        for (size_t a = 0; a != e.agents().size(); ++a)
            REQUIRE( e[s][a] == 1 );    
    REQUIRE( e == f );
    REQUIRE( d == empty );
}

TEST_CASE( "Move constructor and assignment - multilevel", "[system]") {
    dci::system empty;
    dci::system d;
    REQUIRE( d == empty );
    d = dci::system(10, 100, 4);
    REQUIRE( d.agents().size() == 10 );
    REQUIRE( d.samples().size() == 100 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = 15;
    
    dci::system e = d;
    
    dci::system f = move(d);
    for (size_t s = 0; s != e.samples().size(); ++s)
        for (size_t a = 0; a != e.agents().size(); ++a)
            REQUIRE( e[s][a] == 15 );    
    REQUIRE( e == f );
    REQUIRE( d == empty );

}

TEST_CASE( "Agent statistics - multilevel", "[system]") {
    constexpr size_t n_agents = 200;
    constexpr size_t n_samples = 100;
    constexpr size_t bits_per_agent = 8;
    dci::system d(n_agents, n_samples, bits_per_agent);
    REQUIRE( d.samples().size() == n_samples );
    REQUIRE( d.agents().size() == n_agents );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = s;
    d.compute_agent_statistics();
    for (const auto& a : d.agents())
    {

        REQUIRE( a->pdf().size() == n_samples );
        for (const auto& item : a->pdf())
            REQUIRE( item.second == 1 );

    }

}

TEST_CASE( "Agent statistics", "[system]") {
    dci::system d(10, 100);
    REQUIRE( d.samples().size() == 100 );
    REQUIRE( d.agents().size() == 10 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = a % 2;
    d.compute_agent_statistics();
    for (const auto& a : d.agents())
    {
        REQUIRE( a->pdf().size() == 1 );
        REQUIRE( a->pdf().begin()->second == 100 );
    }

}

TEST_CASE( "Exit" ) {
    cout << "Leaving test suite...\n";
    dci::print_allocation_stats(cout);
}
