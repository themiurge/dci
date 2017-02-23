#include "dci.hpp"

int main(int argc, char** argv)
{

    dci::system d = dci::system::load_from_file("data/boolean_10variables.txt");
    d.compute_agent_statistics();
    
    dci::system hs = d.generate_homogeneous_system(123456);
    hs.compute_agent_statistics();
    hs.compute_cluster_index_statistics();
    
    d.load_cluster_index_statistics_from_homogeneous_system(hs);
    
    priority_queue<dci::cluster> results = d.compute_system_statistics(30);
    while (!results.empty())
    {
        cout << results.top() << " --> " << results.top().statistical_index() << endl;
        results.pop();
    }

    return 0;

}

