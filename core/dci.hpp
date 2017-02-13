#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <utility>

using namespace std;

namespace dci
{
    
    // Base register type definition. Each system sample is stored in a number of contiguous reg_type
    // memory cells. Basic constexpr utilities for register handling are also defined.
    typedef unsigned int reg_type;
    constexpr size_t bytes_per_reg = sizeof(reg_type);
    constexpr size_t bits_per_reg = bytes_per_reg * 8;
    constexpr reg_type selection_mask = 1;
    constexpr reg_type zero_mask = 0;
    constexpr reg_type one_mask = ~zero_mask;

    // Register utilities
    namespace register_utils 
    {
        inline size_t number_of_registers(const size_t& number_of_bits) { return (number_of_bits - 1) / bits_per_reg + 1; }
        inline void set_bit(reg_type* reg, const size_t& pos, const bool& value) 
        { 
            if (value) *(reg + pos / bits_per_reg) |= selection_mask << (pos % bits_per_reg); 
            else *(reg + pos / bits_per_reg) &= ~(selection_mask << (pos % bits_per_reg));
        }
        inline bool get_bit(const reg_type* reg, const size_t& pos)
        { 
            return *(reg + pos / bits_per_reg) & (selection_mask << (pos % bits_per_reg));
        }
        inline ostream& write_bits(ostream& out, const reg_type* reg, const size_t& n)
        {
            for (size_t i = 0; i != n; ++i)
                out << get_bit(reg, i);
            return out;
        }
        inline void set_bits_from_string(reg_type* reg, const string& s)
        {
            for (size_t pos = 0; pos != s.size(); ++pos)
                set_bit(reg, pos, s[pos] - '0');
        }
    }

    // Struct containing system properties
    struct system_properties {
        size_t N; // number of agents
        size_t M; // number of samples
        size_t S; // number of registers per sample
        size_t NB; // number of bits per sample
        size_t sample_bytes;
        size_t data_bytes;
        size_t agent_pool_bytes;
        inline void init(const size_t& N, const size_t& M, const size_t& NB)
        {
            this->N = N; this->M = M; this->NB = NB;
            S = register_utils::number_of_registers(NB);
            sample_bytes = S * bytes_per_reg;
            data_bytes = sample_bytes * M;
            agent_pool_bytes = sample_bytes * N;
        }
    };

    // This class represents a dynamic system, defined by N agents and M samples.
    class system;

    class agent_value;

    // This class represents a system agent, defined by its bitmask.
    class agent {
    private:
        system* parent;
        reg_type* bitmask;
        string name;
        agent() : bitmask(nullptr), name() { }
        agent(reg_type* agent_bitmask, const string& agent_name, system* parent_system) : bitmask(agent_bitmask), name(agent_name), parent(parent_system) { }
    public:

        string get_name() { return name; }
        const string& get_name() const { return name; }

        // friend classes
        friend class system;
        friend class sample; 
        friend class agent_value;   

        // friend operators
        friend ostream& operator<<(ostream& out, const agent& a);
        friend istream& operator>>(istream& in, system& sys);    

    };

    // This class representa a single system sample. No data is actually stored here
    // except for the pointer to the first reg_type cell of the sequence containing
    // the sample. Data is stored in the parent system object.
    class sample {
    private:
        system* parent;
        reg_type* data;
        sample() : data(nullptr), parent(nullptr) { }
        sample(reg_type* sample_data, system* parent_system) : data(sample_data), parent(parent_system) { }
    public:

        // access agent value in sample
        agent_value operator[](const size_t& a);

        // friend classes
        friend class system;
        friend class agent;
        friend class agent_value;

        // friend operators
        friend ostream& operator<<(ostream& out, const sample& s);
        friend istream& operator>>(istream& in, system& sys);
        
    };

    // This class holds a (usually) short-lived pointer to an agent's value in a sample.
    class agent_value {
    private:
        sample* s;
        agent* a;
        agent_value(sample* s_in, agent* a_in) : s(s_in), a(a_in) { }
        bool compare_bits(const reg_type* t) const;
        void set_from_bits(const reg_type* t);
    public:

        // NOTE: these operators assume that agent size in bits <= bits_per_reg
        agent_value operator=(const reg_type& t);
        bool operator==(const reg_type& t) const;

        // NOTE: these operators assume that t points to a buffer of at least S reg_type
        agent_value operator=(const reg_type* t);
        bool operator==(const reg_type* t) const;

        // friend classes
        friend class agent;
        friend class sample;
        friend class system;

    };

    class system {
    private:
        vector<sample*> samples;
        vector<agent*> agents;
        reg_type* data;
        reg_type* agent_pool;
        system_properties props;
        inline void reset_fields()
        {
            samples.resize(0);
            agents.resize(0);
            data = nullptr;
            agent_pool = nullptr;
            props.init(0, 0, 0);            
        }
        void free_data(const bool& reset = false)
        {
            for (auto& s : samples) if (s) delete s; for (auto& a : agents) if (a) delete a; if (data) delete[] data; if (agent_pool) delete[] agent_pool;
            if (reset) reset_fields();
        }
        void allocate_data(const size_t& N, const size_t& M, const size_t& NB) 
        {
            props.init(N, M, NB);
            samples.resize(props.M, nullptr);
            agents.resize(props.N, nullptr);
            
            // empty system
            if (!props.N)
            {
                data = nullptr;
                agent_pool = nullptr;
                return;
            }

            data = new reg_type[props.data_bytes];
            agent_pool = new reg_type[props.agent_pool_bytes];

            // reset memory
            memset(data, 0, props.data_bytes);
            memset(agent_pool, 0, props.agent_pool_bytes);

            // set internal pointers for agents and samples
            for (size_t a = 0; a != props.N; ++a)
                agents[a] = new agent(agent_pool + a * props.S, string(1, '0' + a), this);
            for (size_t s = 0; s != props.M; ++s)
                samples[s] = new sample(data + s * props.S, this);
        }
        void copy_from(const system& s)
        {
            allocate_data(s.props.N, s.props.M, s.props.NB);

            // copy raw data and agent pool
            if (data) memcpy(data, s.data, props.data_bytes);
            if (agent_pool) memcpy(agent_pool, s.agent_pool, props.agent_pool_bytes);
        }
        void move_from(system&& s)
        {
            props = s.props;
            data = s.data;
            agent_pool = s.agent_pool;
            samples = move(s.samples);
            agents = move(s.agents);
            for (auto& a : agents) a->parent = this;
            for (auto& a : samples) a->parent = this;

            s.reset_fields();
        }
    public:

        system() : samples(0), agents(0), data(nullptr), agent_pool(nullptr) { props.init(0, 0, 0); }
        system(const size_t& N, const size_t& M, const size_t& bits_per_agent = 1)
        {  
            allocate_data(N, M, N * bits_per_agent);

            // build agent pool
            for (size_t a = 0; a != props.N; ++a)
                // set bits in agent bitmask
                for (size_t b = 0; b != bits_per_agent; ++b)
                    register_utils::set_bit(agent_pool + a * props.S, b + a * bits_per_agent, 1);

        }

        // copy constructor/assignment
        system(const system& s) 
        {
            copy_from(s);
        }
        const system& operator=(const system& s)
        {
            if (&s != this)
            {
                free_data();
                copy_from(s);
            }
            return *this;
        }

        // move constructor/assignment
        system(system&& s) 
        {
            move_from(move(s));
        }
        const system& operator=(system&& s)
        {
            if (&s != this)
            {
                free_data();
                move_from(move(s));
            }
            return *this;
        }

        // number of samples, agents
        inline size_t num_samples() const { return samples.size(); }
        inline size_t num_agents() const { return agents.size(); }

        // access sample by index
        inline sample& operator[](const size_t& index) { if (index >= props.M) throw out_of_range("index"); return *samples[index]; }
        inline const sample& operator[](const size_t& index) const { if (index >= props.M) throw out_of_range("index"); return *samples[index]; }

        // friend classes
        friend class sample;
        friend class agent;
        friend class agent_value;

        // friend operators
        friend ostream& operator<<(ostream& out, const system& sys) 
        { 
            // check if header needs to be printed
            if (sys.props.NB > sys.props.N)
            {
                out << "%%\n";
                for (const auto& a : sys.agents)
                    out << (*a) << endl;
                out << "%%\n";
            }

            // print rows
            for (const auto& s : sys.samples) 
                out << (*s) << endl; 
            return out; 
        }
        friend ostream& operator<<(ostream& out, const sample& s);
        friend ostream& operator<<(ostream& out, const agent& a);
        friend bool operator==(const system& s1, const system& s2);
        friend bool operator!=(const system& s1, const system& s2);
        friend istream& operator>>(istream& in, system& s);

        // save functions
        void save_to_file(const string& file_path)
        {
            ofstream out(file_path);
            out << *this;
        }

        // homogeneous system generation
        system generate_homogeneous_system()
        {


            system s;
            s.allocate_data(props.N, props.M, props.NB);

            // just copy agent pool
            if (agent_pool) memcpy(s.agent_pool, agent_pool, props.agent_pool_bytes);

            return s;

        }

        virtual ~system() { free_data(); }
    };
    
    //
    // agent_value class implementation
    //
    inline void agent_value::set_from_bits(const reg_type* t)
    {
        size_t pt = 0;
        for (size_t p = 0; p != s->parent->props.NB; ++p)
            if (register_utils::get_bit(a->bitmask, p))
                register_utils::set_bit(s->data, p, register_utils::get_bit(t, pt++));
    }

    inline bool agent_value::compare_bits(const reg_type* t) const
    {
        size_t pt = 0;
        for (size_t p = 0; p != s->parent->props.NB; ++p)
            if (register_utils::get_bit(a->bitmask, p))
                if (register_utils::get_bit(t, pt++) != register_utils::get_bit(s->data, p)) return false;
        return true;
    }

    inline agent_value agent_value::operator=(const reg_type& t)
    {
        set_from_bits(&t);
        return *this;
    }

    inline bool agent_value::operator==(const reg_type& t) const
    {
        return compare_bits(&t);
    }

    inline agent_value agent_value::operator=(const reg_type* t)
    {
        set_from_bits(t);
        return *this;
    }

    inline bool agent_value::operator==(const reg_type* t) const
    {
        return compare_bits(t);
    }

    //
    // agent class implementation
    //
    ostream& operator<<(ostream& out, const agent& a) 
    { 
        return register_utils::write_bits(out, a.bitmask, (a.parent)->props.NB);
    }

    //
    // sample class implementation
    //
    agent_value sample::operator[](const size_t& a)
    {
        return agent_value(this, parent->agents[a]);
    }

    ostream& operator<<(ostream& out, const sample& s) 
    { 
        return register_utils::write_bits(out, s.data, (s.parent)->props.NB);
    }

    //
    // system class implementation
    //
    inline bool operator==(const system& s1, const system& s2)
    {
        return s1.props.N == s2.props.N && s1.props.M == s2.props.M && s1.props.NB == s1.props.NB && 
            !memcmp(s1.data, s2.data, s1.props.data_bytes) &&
            !memcmp(s1.agent_pool, s2.agent_pool, s1.props.agent_pool_bytes);
    }

    inline bool operator!=(const system& s1, const system& s2)
    {
        return s1.props.N != s2.props.N || s1.props.M != s2.props.M || s1.props.NB != s1.props.NB || 
            memcmp(s1.data, s2.data, s1.props.data_bytes) ||
            memcmp(s1.agent_pool, s2.agent_pool, s1.props.agent_pool_bytes);
    }

    istream& operator>>(istream& in, system& s)
    {

        size_t N = 0, M = 0, NB = 0;
        vector<string> lines;
        string line;

        // free all data in previous system object
        s.free_data(true);

        // store all lines of input file in memory
        while (in >> line)
            lines.emplace_back(line);

        // index of first data line
        size_t offset = 0;

        // check header
        if (lines[0][0] == '%' && lines[0][1] == '%')
            for (offset = 2; lines[offset - 1][0] != '%'; ++offset);

        NB = lines[offset].size();
        N = offset ? (offset - 2) : NB;
        M = lines.size() - offset;
        s.allocate_data(N, M, NB);

        // build agent pool
        for (size_t a = 0; a != s.props.N; ++a)
        {
            // create agent pointing to corresponding agent pool position
            reg_type* cur_agent = s.agent_pool + a * s.props.S;
            s.agents[a] = new agent(cur_agent, string(1, '0' + a), &s);

            // set bits in agent bitmask
            if (offset)
                register_utils::set_bits_from_string((s.agents[a])->bitmask, lines[a + 1]);
            else
                register_utils::set_bit(cur_agent, a, 1);
        }

        // build sample pool
        for (size_t i = 0; i != s.props.M; ++i)
        {
            s.samples[i] = new sample(s.data + i * s.props.S, &s);
            register_utils::set_bits_from_string((s.samples[i])->data, lines[i + offset]);
        }

        return in;
    }

    //
    // System factories
    //

    // load from file
    system load_from_file(const string& file_path)
    {
        ifstream in(file_path);
        system s;
        in >> s;
        return s;
    }

}
