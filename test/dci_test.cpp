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
        REQUIRE( a->entropy() == Approx(0.0) );
    }

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
        REQUIRE( a->entropy() == Approx(std::log2((dci::fp_type)n_samples)) );
    }

}

TEST_CASE( "Homogeneous system", "[system]" ) {
    dci::system d(10, 100);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = a % 2;
    d.compute_agent_statistics();
    
    dci::system hs1 = d.generate_homogeneous_system(123456); 
    dci::system hs2 = d.generate_homogeneous_system(654321); 
    REQUIRE( d == hs1 );
    REQUIRE( d == hs2 );

}

TEST_CASE( "Homogeneous system - multilevel", "[system]" ) {
    constexpr size_t n_agents = 200;
    constexpr size_t n_samples = 100;
    constexpr size_t bits_per_agent = 8;
    dci::system d(n_agents, n_samples, bits_per_agent);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = (s < n_samples / 2) ? 1 : 0;
    d.compute_agent_statistics();

    dci::system hs1 = d.generate_homogeneous_system(123456); 
    dci::system hs2 = d.generate_homogeneous_system(453664536); 
    REQUIRE( hs1.samples().size() == n_samples );
    REQUIRE( hs1.agents().size() == n_agents );
    REQUIRE( hs2.samples().size() == n_samples );
    REQUIRE( hs2.agents().size() == n_agents );

    /*
    hs1.compute_agent_statistics();
    hs2.compute_agent_statistics();
    cout << d.agents()[0]->pdf() << endl;
    cout << hs1.agents()[0]->pdf() << endl;
    cout << hs2.agents()[0]->pdf() << endl;
    */

    REQUIRE( hs1 != hs2 );

}

TEST_CASE( "Real life system - arabidopsis", "[system]" ) {
    dci::system d = dci::system::load_from_file("data/arabidopsis.txt");
    REQUIRE( d.samples().size() == 5000 );
    REQUIRE( d.agents().size() == 15 );

    for (const auto& a : d.agents())
        REQUIRE( a->size() == 1 );
}

TEST_CASE( "Clusters", "[cluster]" ) {
    dci::system d(10, 100);
    dci::cluster c1(&d, { 1, 3, 5 });
    dci::cluster c2(&d, "0101010000");
    dci::reg_type t = 42;
    dci::cluster c3(&d, &t);
    dci::cluster c4(&d, { 0, 2, 4, 6, 7, 8, 9 } );
    dci::cluster c5 = !c2;
    REQUIRE( c1 == c2 );
    REQUIRE( c2 == c3 );
    REQUIRE( c4 == c5 );
    REQUIRE( c5 == !c3 );
    REQUIRE( c4 != c1 );
    REQUIRE( c5 != c2 );
}

TEST_CASE( "Clusters - multilevel", "[cluster]" ) {
    dci::system d(10, 100, 47);
    dci::cluster c1(&d, { 1, 3, 5 });
    dci::cluster c2(&d, "0101010000");
    dci::reg_type t = 42;
    dci::cluster c3(&d, &t);
    dci::cluster c4(&d, { 0, 2, 4, 6, 7, 8, 9 } );
    dci::cluster c5 = !c2;
    REQUIRE( c1 == c2 );
    REQUIRE( c2 == c3 );
    REQUIRE( c4 == c5 );
    REQUIRE( c5 == !c3 );
    REQUIRE( c4 != c1 );
    REQUIRE( c5 != c2 );
}

TEST_CASE( "Cluster histograms", "[cluster]" ) {
    dci::system d(10, 100);
    REQUIRE( d.samples().size() == 100 );
    REQUIRE( d.agents().size() == 10 );
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = a % 2;

    dci::cluster c1(&d, { 1, 3, 5 });
    dci::cluster c2 = !c1;

    REQUIRE( c1.histogram().size() == 1 );
    REQUIRE( c1.histogram().begin()->second == 100 );
    REQUIRE( c1.entropy() == Approx(0.0) );

    REQUIRE( c2.histogram().size() == 1 );
    REQUIRE( c2.histogram().begin()->second == 100 );
    REQUIRE( c2.entropy() == Approx(0.0) );

}

TEST_CASE( "Cluster histograms - multilevel", "[system]") {
    constexpr size_t n_agents = 200;
    constexpr size_t n_samples = 100;
    constexpr size_t bits_per_agent = 8;
    dci::system d(n_agents, n_samples, bits_per_agent);
    for (size_t s = 0; s != d.samples().size(); ++s)
        for (size_t a = 0; a != d.agents().size(); ++a)
            d[s][a] = s;

    dci::cluster c1(&d, { 1, 3, 5 });
    dci::cluster c2 = !c1;

    REQUIRE( c1.histogram().size() == n_samples );
    for (const auto& item : c1.histogram())
        REQUIRE( item.second == 1 );
    REQUIRE( c1.entropy() == Approx(std::log2((dci::fp_type)n_samples)) );

    REQUIRE( c2.histogram().size() == n_samples );
    for (const auto& item : c2.histogram())
        REQUIRE( item.second == 1 );
    REQUIRE( c2.entropy() == Approx(std::log2((dci::fp_type)n_samples)) );

}

TEST_CASE( "Exit" ) {
    cout << "Leaving test suite...\n";
    dci::print_allocation_stats(cout);
}
