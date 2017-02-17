#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <cmath>
#include <random>

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

    // Floating point operations type definition
    typedef float fp_type;

    // Register utilities
    namespace register_utils 
    {

        // defines for murmur3 hash
        #define m3_c1 0xcc9e2d51
        #define m3_c2 0x1b873593
        #define m3_r1 15
        #define m3_r2 13
        #define m3_m 5
        #define m3_n 0xe6546b64
        #define m3_seed 0x422a0b8f
        #define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effort

        constexpr inline size_t number_of_registers(const size_t& number_of_bits) { return (number_of_bits - 1) / bits_per_reg + 1; }
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
        template<bool b> inline size_t count(const reg_type& reg)        
        {
            size_t n = 0;
            for(reg_type mask = selection_mask; mask; mask <<= 1)
                n += b ? ((reg & mask) != 0) : ((reg & mask) == 0);
            return n;
        }
        template<bool b> inline size_t count(const reg_type* reg, const size_t& s)
        {
            size_t n = 0;
            for (size_t i = 0; i != s; ++i)
                n += count<b>(*(reg++));
            return n;
        }

        //
        // template, template specializations and general function for comparison
        //
        inline bool equal(const reg_type* a, const reg_type* b, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) if (a[i] != b[i]) return false; return true;
        }
        template<size_t S> bool equal(const reg_type* a, const reg_type* b)
        {
            return equal(a, b, S);
        }
        template<> inline bool equal<1>(const reg_type* a, const reg_type* b) { return a[0] == b[0]; }
        template<> inline bool equal<2>(const reg_type* a, const reg_type* b) { return a[0] == b[0] && a[1] == b[1]; }
        template<> inline bool equal<3>(const reg_type* a, const reg_type* b) 
        { 
            /*write_bits(cout, a, 3 * bits_per_reg) << endl; write_bits(cout, b, 3 * bits_per_reg) << endl << endl;*/ 
            //write_bits(cout, a, 3 * bits_per_reg) << " equal" << endl;
            return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; 
        }
        template<> inline bool equal<4>(const reg_type* a, const reg_type* b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]; }

        
        //
        // template, template specializations and general function for assignment
        //
        inline void assign(reg_type* dest, const reg_type* src, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] = src[i];
        }
        template<size_t S> inline void assign(reg_type* dest, const reg_type* src)
        {
            assign(dest, src, S);
        }
        template<> inline void assign<1>(reg_type* dest, const reg_type* src) { dest[0] = src[0]; }
        template<> inline void assign<2>(reg_type* dest, const reg_type* src) { dest[0] = src[0]; dest[1] = src[1]; }
        template<> inline void assign<3>(reg_type* dest, const reg_type* src) { dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2]; }
        template<> inline void assign<4>(reg_type* dest, const reg_type* src) { dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2]; dest[3] = src[3]; }

        //
        // template, template specializations and general function for combination (bitwise or)
        //
        inline void combine(reg_type* dest, const reg_type* src, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] |= src[i];
        }
        template<size_t S> inline void combine(reg_type* dest, const reg_type* src)
        {
            combine(dest, src, S);
        }
        template<> inline void combine<1>(reg_type* dest, const reg_type* src) { dest[0] |= src[0]; }
        template<> inline void combine<2>(reg_type* dest, const reg_type* src) { dest[0] |= src[0]; dest[1] |= src[1]; }
        template<> inline void combine<3>(reg_type* dest, const reg_type* src) { dest[0] |= src[0]; dest[1] |= src[1]; dest[2] |= src[2]; }
        template<> inline void combine<4>(reg_type* dest, const reg_type* src) { dest[0] |= src[0]; dest[1] |= src[1]; dest[2] |= src[2]; dest[3] |= src[3]; }

        //
        // template, template specializations and general function for mask-limited assignment
        //
        inline void assign_in_mask(reg_type* dest, const reg_type* src, const reg_type* mask, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] = (dest[i] & ~mask[i]) | src[i];
        }
        template<size_t S> inline void assign_in_mask(reg_type* dest, const reg_type* src, const reg_type* mask)
        {
            assign_in_mask(dest, src, mask, S);
        }
        template<> inline void assign_in_mask<1>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = (dest[0] & ~mask[0]) | src[0]; }
        template<> inline void assign_in_mask<2>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = (dest[0] & ~mask[0]) | src[0]; dest[1] = (dest[1] & ~mask[1]) | src[1]; }
        template<> inline void assign_in_mask<3>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = (dest[0] & ~mask[0]) | src[0]; dest[1] = (dest[1] & ~mask[1]) | src[1]; dest[2] = (dest[2] & ~mask[2]) | src[2]; }
        template<> inline void assign_in_mask<4>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = (dest[0] & ~mask[0]) | src[0]; dest[1] = (dest[1] & ~mask[1]) | src[1]; dest[2] = (dest[2] & ~mask[2]) | src[2]; dest[3] = (dest[3] & ~mask[3]) | src[3]; }

        //
        // template, template specializations and general function for mask assignment
        //
        inline void assign_from_mask(reg_type* dest, const reg_type* src, const reg_type* mask, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] = src[i] & mask[i];
        }
        template<size_t S> inline void assign_from_mask(reg_type* dest, const reg_type* src, const reg_type* mask)
        {
            assign_from_mask(dest, src, mask, S);
        }
        template<> inline void assign_from_mask<1>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = src[0] & mask[0]; }
        template<> inline void assign_from_mask<2>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = src[0] & mask[0]; dest[1] = src[1] & mask[1]; }
        template<> inline void assign_from_mask<3>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = src[0] & mask[0]; dest[1] = src[1] & mask[1]; dest[2] = src[2] & mask[2]; }
        template<> inline void assign_from_mask<4>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] = src[0] & mask[0]; dest[1] = src[1] & mask[1]; dest[2] = src[2] & mask[2]; dest[3] = src[3] & mask[3]; }

        //
        // template, template specializations and general function for mask combination
        //
        inline void combine_from_mask(reg_type* dest, const reg_type* src, const reg_type* mask, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] |= src[i] & mask[i];
        }
        template<size_t S> inline void combine_from_mask(reg_type* dest, const reg_type* src, const reg_type* mask)
        {
            combine_from_mask(dest, src, mask, S);
        }
        template<> inline void combine_from_mask<1>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] |= src[0] & mask[0]; }
        template<> inline void combine_from_mask<2>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] |= src[0] & mask[0]; dest[1] |= src[1] & mask[1]; }
        template<> inline void combine_from_mask<3>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] |= src[0] & mask[0]; dest[1] |= src[1] & mask[1]; dest[2] |= src[2] & mask[2]; }
        template<> inline void combine_from_mask<4>(reg_type* dest, const reg_type* src, const reg_type* mask) { dest[0] |= src[0] & mask[0]; dest[1] |= src[1] & mask[1]; dest[2] |= src[2] & mask[2]; dest[3] |= src[3] & mask[3]; }

        //
        // hashing functions
        //
        inline unsigned int murmur3_hash(const reg_type* reg, const size_t& S)
        {
            char* key = (char*)reg;
            unsigned int hash = m3_seed;
            unsigned int len = S * bytes_per_reg;

            const int nblocks = len / 4;
            const unsigned int *blocks = (const unsigned int *) key;
            int i;
            unsigned int k;
            for (i = 0; i < nblocks; i++) {
                k = blocks[i];
                k *= m3_c1;
                k = ROT32(k, m3_r1);
                k *= m3_c2;

                hash ^= k;
                hash = ROT32(hash, m3_r2) * m3_m + m3_n;
            }

            const unsigned char *tail = (const unsigned char *) (key + nblocks * 4);
            unsigned int k1 = 0;

            switch (len & 3) {
            case 3:
                k1 ^= tail[2] << 16;
            case 2:
                k1 ^= tail[1] << 8;
            case 1:
                k1 ^= tail[0];

                k1 *= m3_c1;
                k1 = ROT32(k1, m3_r1);
                k1 *= m3_c2;
                hash ^= k1;
            }

            hash ^= len;
            hash ^= (hash >> 16);
            hash *= 0x85ebca6b;
            hash ^= (hash >> 13);
            hash *= 0xc2b2ae35;
            hash ^= (hash >> 16);

            return hash;
        }
        inline reg_type basic_hash(const reg_type* reg, const size_t& S) { reg_type h = reg[0]; for (size_t i = 1; i < S; ++i) h += reg[i]; return h; }

    }

    // Class for murmur3 hashing
    class murmur3_hasher {
    private:
        size_t _s;
        unsigned int (*_f)(const reg_type*, const size_t&);
    public:
        murmur3_hasher(const size_t& s = 0) : _s(s) { _f = register_utils::murmur3_hash; }
        inline unsigned int operator()(const reg_type* reg) const { 
            return _f(reg, _s); 
            //unsigned int hash = _f(reg, _s); 
            //register_utils::write_bits(cout, reg, _s * bits_per_reg) << " -> " << hash << endl;
            //return hash;
        }
    };

    // Class for basic hashing
    class basic_hasher {
    private:
        size_t _s;
        reg_type (*_f)(const reg_type*, const size_t&);
    public:
        basic_hasher(const size_t& s = 0) : _s(s) { _f = register_utils::basic_hash; }
        inline reg_type operator()(const reg_type* reg) const { return _f(reg, _s); }
    };

    // Class for register comparison
    class register_compare {
    private:
        size_t _s;
        bool (*_f)(const reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static bool _wrapper(const reg_type* a, const reg_type* b, const size_t& s) { return register_utils::equal<S>(a, b); }
    public:
        register_compare(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::equal; break;
            }
        }
        inline bool operator()(const reg_type* a, const reg_type* b) const { return _f(a, b, _s); }
    };

    // Class for register assignment
    class register_assign {
    private:
        size_t _s;
        void (*_f)(reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* a, const reg_type* b, const size_t& s) { register_utils::assign<S>(a, b); }
    public:
        register_assign(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::assign; break;
            }
        }
        inline void operator()(reg_type* a, const reg_type* b) const { _f(a, b, _s); }
    };

    // Class for register combination
    class register_combine {
    private:
        size_t _s;
        void (*_f)(reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* a, const reg_type* b, const size_t& s) { register_utils::combine<S>(a, b); }
    public:
        register_combine(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::combine; break;
            }
        }
        inline void operator()(reg_type* a, const reg_type* b) const { _f(a, b, _s); }
    };

    // Class for register mask assignment
    class register_assign_from_mask {
    private:
        size_t _s;
        void (*_f)(reg_type*, const reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* a, const reg_type* b, const reg_type* mask, const size_t& s) { register_utils::assign_from_mask<S>(a, b, mask); }
    public:
        register_assign_from_mask(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::assign_from_mask; break;
            }
        }
        inline void operator()(reg_type* a, const reg_type* b, const reg_type* mask) const { _f(a, b, mask, _s); }
    };

    // Class for register mask-delimited assignment
    class register_assign_in_mask {
    private:
        size_t _s;
        void (*_f)(reg_type*, const reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* a, const reg_type* b, const reg_type* mask, const size_t& s) { register_utils::assign_in_mask<S>(a, b, mask); }
    public:
        register_assign_in_mask(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::assign_in_mask; break;
            }
        }
        inline void operator()(reg_type* a, const reg_type* b, const reg_type* mask) const { _f(a, b, mask, _s); }
    };

    // Class for register mask combination
    class register_combine_from_mask {
    private:
        size_t _s;
        void (*_f)(reg_type*, const reg_type*, const reg_type*, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* a, const reg_type* b, const reg_type* mask, const size_t& s) { register_utils::combine_from_mask<S>(a, b, mask); }
    public:
        register_combine_from_mask(const size_t& s = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::combine_from_mask; break;
            }
        }
        inline void operator()(reg_type* a, const reg_type* b, const reg_type* mask) const { _f(a, b, mask, _s); }
    };

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

    // Struct containing a full cluster/agent register, used for histograms
    class reg 
    {
    public:
        static size_t n_new;
        static size_t n_del;
    private:
        reg_type* _data;
        size_t _s;
        void reset_fields()
        {
            _s = 0;
            _data = nullptr;
        }
        void copy_from(const reg_type* data, const size_t& s)
        {
            _s = s;
            if (_s)
            {
                _data = new reg_type[_s]; reg::n_new++;
                memcpy(_data, data, _s * bytes_per_reg);
            }
            else
                _data = nullptr;
        }
        void free_data(const bool& reset = false)
        {
            if (_s) { delete[] _data; reg::n_del++; }
            if (reset) reset_fields();
        }
        void move_from(reg&& r)
        {
            _s = r._s;
            _data = r._data;

            r.reset_fields();
        }
    public:
        reg() : _data(nullptr), _s(0) { }
        reg(const reg_type* data, const size_t& s) { copy_from(data, s); }

        reg(const reg& r) { copy_from(r._data, r._s); }
        inline reg& operator=(const reg& r) { if (this != &r) { free_data(); copy_from(r._data, r._s); } return *this; }
        reg(reg&& r) { move_from(move(r)); }
        inline reg& operator=(reg&& r) { if (this != &r) { free_data(); move_from(move(r)); } return *this; }

        inline operator const reg_type*() const { return _data; }
        inline operator reg_type*() { return _data; }

        friend ostream& operator<<(ostream&, const reg&);

        virtual ~reg() { free_data(); }
    };

    // This class represents a dynamic system, defined by N agents and M samples.
    class system;

    // This class holds a short-lived agent value inside a sample
    class agent_value;

    // This class represents a group of agents
    class cluster;

    // This is just a typedef for a categorical distribution
    template<class hash_func> using hashmap = unordered_map<reg, size_t, hash_func, register_compare>;
    using categorical_distribution = hashmap<murmur3_hasher>;

    // Template function used to compute entropy over a histogram
    template<class hash_func> 
    fp_type entropy(const hashmap<hash_func>& histogram, const size_t& tot)
    {
        if (!tot) return 0.0;
        fp_type ent = 0.0;
        for (const auto& elem : histogram)
            if (elem.second) 
            {
                fp_type freq = (fp_type)(elem.second) / (fp_type)tot;
                ent -= freq * log2(freq);
            }
        return ent;
    }

    // Template function used to compute a histogram
    template<class hash_func> hashmap<hash_func> compute_histogram(const system&, const reg_type*, const hash_func&);

    // This class represents a system agent, defined by its bitmask.
    class agent {
    public:
        static size_t n_new;
        static size_t n_del;
    private:
        system* _parent;
        reg_type* _bitmask;
        string _name;
        size_t _id;
        categorical_distribution _pdf;
        agent() : _bitmask(nullptr), _name(), _pdf(), _id(0) { }
        agent(reg_type* agent_bitmask, const string& agent_name, const size_t& id, system* parent_system) : 
            _bitmask(agent_bitmask), _name(agent_name), _parent(parent_system), _pdf(), _id(id)
        { agent::n_new++; }
        void compute_pdf();
    public:

        inline const string& name() const { return _name; }
        inline const size_t& id() const { return _id; }
        inline const categorical_distribution& pdf() const { return _pdf; }
        fp_type entropy() const;
        size_t size() const;

        // friend classes
        friend class system;
        friend class sample; 
        friend class agent_value;   
        friend class cluster;   

        // friend operators
        friend ostream& operator<<(ostream& out, const agent& a);
        friend istream& operator>>(istream& in, system& sys); 

        virtual ~agent() { if (_parent) agent::n_del++; }   

    };

    // This class represents a single system sample. No data is actually stored here
    // except for the pointer to the first reg_type cell of the sequence containing
    // the sample. Data is stored in the parent system object.
    class sample {
    public:
        static size_t n_new;
        static size_t n_del;
    private:
        system* _parent;
        reg_type* _data;
        sample() : _data(nullptr), _parent(nullptr) { }
        sample(reg_type* sample_data, system* parent_system) : _data(sample_data), _parent(parent_system) { sample::n_new++; }
    public:

        // access agent value in sample
        agent_value operator[](const size_t& a);

        // friend classes
        friend class system;
        friend class agent;
        friend class agent_value;
        friend class cluster;   

        // friend template functions
        template<class hash_func> friend hashmap<hash_func> compute_histogram(const system&, const reg_type*, const hash_func&);

        // friend operators
        friend ostream& operator<<(ostream& out, const sample& s);
        friend istream& operator>>(istream& in, system& sys);
        
        virtual ~sample() { if (_parent) sample::n_del++; }   

    };

    // This class holds a (usually) short-lived pointer to an agent's value in a sample.
    class agent_value {
    private:
        sample* _s;
        agent* _a;
        agent_value(sample* s, agent* a) : _s(s), _a(a) { }
        bool compare_bits(const reg_type* t) const;
        void set_from_bits(const reg_type* t);
        void set_from_sample_value(const reg_type* t);
    public:

        // NOTE: these operators assume that agent size in bits <= bits_per_reg
        template<typename T> agent_value operator=(T t);
        template<typename T> bool operator==(T t) const;

        // friend classes
        friend class agent;
        friend class sample;
        friend class system;
        friend class cluster;   

    };

    class cluster {
    private:

    };

    class system {
    public:
        static size_t n_new;
        static size_t n_del;
    private:
        vector<sample*> _samples;
        vector<agent*> _agents;
        reg_type* _data;
        reg_type* _agent_pool;
        system_properties _props;
        register_compare _equal;
        register_assign _assign;
        register_assign_from_mask _assign_from_mask;
        register_assign_in_mask _assign_in_mask;
        register_combine _combine;
        register_combine_from_mask _combine_from_mask;
        inline void reset_fields()
        {
            _samples.resize(0);
            _agents.resize(0);
            _data = nullptr;
            _agent_pool = nullptr;
            _props.init(0, 0, 0);            
        }
        void free_data(const bool& reset = false)
        {
            for (auto& s : _samples) if (s) delete s; for (auto& a : _agents) if (a) delete a; 
            if (_data) { delete[] _data; system::n_del++; }
            if (_agent_pool) delete[] _agent_pool;
            if (reset) reset_fields();
        }
        void allocate_data(const size_t& N, const size_t& M, const size_t& NB) 
        {
            _props.init(N, M, NB);
            _samples.resize(_props.M, nullptr);
            _agents.resize(_props.N, nullptr);
            _equal = register_compare(_props.S);
            _assign = register_assign(_props.S);
            _assign_from_mask = register_assign_from_mask(_props.S);
            _assign_in_mask = register_assign_in_mask(_props.S);
            _combine = register_combine(_props.S);
            _combine_from_mask = register_combine_from_mask(_props.S);
            
            // empty system
            if (!_props.N)
            {
                _data = nullptr;
                _agent_pool = nullptr;
                return;
            }

            _data = new reg_type[_props.data_bytes]; system::n_new++;
            _agent_pool = new reg_type[_props.agent_pool_bytes];

            // reset memory
            memset(_data, 0, _props.data_bytes);
            memset(_agent_pool, 0, _props.agent_pool_bytes);

            // set internal pointers for agents and samples
            for (size_t a = 0; a != _props.N; ++a)
                _agents[a] = new agent(_agent_pool + a * _props.S, string(1, '0' + a), a, this);
            for (size_t s = 0; s != _props.M; ++s)
                _samples[s] = new sample(_data + s * _props.S, this);
        }
        void copy_from(const system& s)
        {
            allocate_data(s._props.N, s._props.M, s._props.NB);

            // copy raw data and agent pool
            if (_data) memcpy(_data, s._data, _props.data_bytes);
            if (_agent_pool) memcpy(_agent_pool, s._agent_pool, _props.agent_pool_bytes);
        }
        void move_from(system&& s)
        {
            _props = s._props;
            _data = s._data;
            _agent_pool = s._agent_pool;
            _samples = move(s._samples);
            _agents = move(s._agents);
            for (auto& a : _agents) a->_parent = this;
            for (auto& a : _samples) a->_parent = this;

            s.reset_fields();
        }
        size_t get_max_values(const size_t& nbits) const
        {
            if (nbits > sizeof(size_t) * 8 || (1 << nbits) > _props.M) return _props.M;
            return 1 << nbits;
        }
    public:

        system() : _samples(0), _agents(0), _data(nullptr), _agent_pool(nullptr) { _props.init(0, 0, 0); }
        system(const size_t& N, const size_t& M, const size_t& bits_per_agent = 1)
        {  
            allocate_data(N, M, N * bits_per_agent);

            // build agent pool
            for (size_t a = 0; a != _props.N; ++a)
                // set bits in agent bitmask
                for (size_t b = 0; b != bits_per_agent; ++b)
                    register_utils::set_bit(_agent_pool + a * _props.S, b + a * bits_per_agent, 1);

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

        // access agents and samples
        const vector<sample*>& samples() const { return _samples; }
        const vector<agent*>& agents() const { return _agents; }

        // access sample by index
        inline sample& operator[](const size_t& index) { if (index >= _props.M) throw out_of_range("index"); return *_samples[index]; }
        inline const sample& operator[](const size_t& index) const { if (index >= _props.M) throw out_of_range("index"); return *_samples[index]; }

        // friend classes
        friend class sample;
        friend class agent;
        friend class agent_value;
        friend class cluster;   

        // friend template functions
        template<class hash_func> friend hashmap<hash_func> compute_histogram(const system&, const reg_type*, const hash_func&);

        // friend operators
        friend ostream& operator<<(ostream& out, const system& sys) 
        { 
            // check if header needs to be printed
            if (sys._props.NB > sys._props.N)
            {
                out << "%%\n";
                for (const auto& a : sys._agents)
                    out << (*a) << endl;
                out << "%%\n";
            }

            // print rows
            for (const auto& s : sys._samples) 
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

        //
        // clusters
        //


        // homogeneous system generation
        system generate_homogeneous_system(const unsigned int& random_seed)
        {


            system hs;
            hs.allocate_data(_props.N, _props.M, _props.NB);

            // just copy agent pool
            if (_agent_pool) memcpy(hs._agent_pool, _agent_pool, _props.agent_pool_bytes);

            // setup random generator
            mt19937 rng(random_seed); // set random seed
            uniform_int_distribution<size_t> gen(0, _props.M - 1); // setup random number generator  

            // set sample value for each agent according to original distribution
            for (auto& s : hs._samples)
                for (size_t a = 0; a != _agents.size(); ++a)
                {

                    // get a random number between 0 and M-1
                    size_t val = gen(rng);

                    // retrieve corresponding value in distribution
                    for (const auto& elem : _agents[a]->pdf())
                    {
                        if (val < elem.second)
                        {
                            _combine_from_mask(s->_data, elem.first, _agents[a]->_bitmask);
                            break;
                        }
                        val -= elem.second;
                    }

                }


            return hs;

        }

        // stats
        void compute_agent_statistics()
        {
            for (auto& a : _agents) a->compute_pdf();
        }

        //
        // System factories
        //

        // load from file
        static system load_from_file(const string& file_path)
        {
            ifstream in(file_path);
            system s;
            in >> s;
            return s;
        }

        virtual ~system() { free_data(); }
    };
    
    //
    // agent_value class implementation
    //
    inline void agent_value::set_from_sample_value(const reg_type* t)
    {
        _s->_parent->_assign_in_mask(_s->_data, t, _a->_bitmask);
    }

    inline void agent_value::set_from_bits(const reg_type* t)
    {
        size_t pt = 0;
        for (size_t p = 0; p != _s->_parent->_props.NB; ++p)
            if (register_utils::get_bit(_a->_bitmask, p))
                register_utils::set_bit(_s->_data, p, register_utils::get_bit(t, pt++));
    }

    inline bool agent_value::compare_bits(const reg_type* t) const
    {
        size_t pt = 0;
        for (size_t p = 0; p != _s->_parent->_props.NB; ++p)
            if (register_utils::get_bit(_a->_bitmask, p))
                if (register_utils::get_bit(t, pt++) != register_utils::get_bit(_s->_data, p)) return false;
        return true;
    }

    template<typename T> inline agent_value agent_value::operator=(T t)
    {
        reg_type temp = static_cast<reg_type>(t);
        set_from_bits(&temp);
        return *this;
    }

    template<> inline agent_value agent_value::operator=(reg_type t)
    {
        set_from_bits(&t);
        return *this;
    }

    template<> inline agent_value agent_value::operator=(reg_type* t)
    {
        set_from_bits(t);
        return *this;
    }

    template<typename T> inline bool agent_value::operator==(T t) const
    {
        reg_type temp = static_cast<reg_type>(t);
        return compare_bits(&temp);
    }

    template<> inline bool agent_value::operator==(reg_type t) const
    {
        return compare_bits(&t);
    }

    template<> inline bool agent_value::operator==(reg_type* t) const
    {
        return compare_bits(t);
    }

    //
    // agent class implementation
    //
    inline size_t agent::size() const { return register_utils::count<1>(_bitmask, _parent->_props.S); }

    ostream& operator<<(ostream& out, const agent& a) 
    { 
        return register_utils::write_bits(out, a._bitmask, (a._parent)->_props.NB);
    }

    void agent::compute_pdf()
    {
        _pdf = compute_histogram(*_parent, _bitmask, murmur3_hasher(_parent->_props.S));
    }

    inline fp_type agent::entropy() const { return dci::entropy(_pdf, _parent->_props.M); }

    //
    // sample class implementation
    //
    agent_value sample::operator[](const size_t& a)
    {
        return agent_value(this, _parent->_agents[a]);
    }

    ostream& operator<<(ostream& out, const sample& s) 
    { 
        return register_utils::write_bits(out, s._data, (s._parent)->_props.NB);
    }

    //
    // cluster class implementation
    //


    //
    // system class implementation
    //
    inline bool operator==(const system& s1, const system& s2)
    {
        return s1._props.N == s2._props.N && s1._props.M == s2._props.M && s1._props.NB == s1._props.NB && 
            !memcmp(s1._data, s2._data, s1._props.data_bytes) &&
            !memcmp(s1._agent_pool, s2._agent_pool, s1._props.agent_pool_bytes);
    }

    inline bool operator!=(const system& s1, const system& s2)
    {
        return s1._props.N != s2._props.N || s1._props.M != s2._props.M || s1._props.NB != s1._props.NB || 
            memcmp(s1._data, s2._data, s1._props.data_bytes) ||
            memcmp(s1._agent_pool, s2._agent_pool, s1._props.agent_pool_bytes);
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
        for (size_t a = 0; a != s._props.N; ++a)
        {
            // create agent pointing to corresponding agent pool position
            reg_type* cur_agent = s._agent_pool + a * s._props.S;
            
            // set bits in agent bitmask
            if (offset)
                register_utils::set_bits_from_string((s._agents[a])->_bitmask, lines[a + 1]);
            else
                register_utils::set_bit(cur_agent, a, 1);
        }

        // build sample pool
        for (size_t i = 0; i != s._props.M; ++i)
            register_utils::set_bits_from_string((s._samples[i])->_data, lines[i + offset]);

        return in;
    }

    // memory leak detector helpers
    size_t reg::n_new = 0;
    size_t reg::n_del = 0;
    size_t system::n_new = 0;
    size_t system::n_del = 0;
    size_t agent::n_new = 0;
    size_t agent::n_del = 0;
    size_t sample::n_new = 0;
    size_t sample::n_del = 0;

    ostream& print_allocation_stats(ostream& out)
    {
        out << "reg class performed " << dci::reg::n_new << " new[] and " << dci::reg::n_del << " delete[] statements\n";
        out << "system class performed " << dci::system::n_new << " new[] and " << dci::system::n_del << " delete[] statements\n";
        out << "agent class performed " << dci::agent::n_new << " new and " << dci::agent::n_del << " delete statements\n";
        out << "sample class performed " << dci::sample::n_new << " new and " << dci::sample::n_del << " delete statements\n"; 
        return out;
    }

    // miscellaneous operators
    ostream& operator<<(ostream& out, const reg& r)
    {
        return register_utils::write_bits(out, r._data, r._s * bits_per_reg);
    }

    ostream& operator<<(ostream& out, const categorical_distribution& pdf)
    {
        for(const auto& item : pdf) out << item.first << " -> " << item.second << endl;
            return out;
    }

    // utility functions
    template<class hash_func> hashmap<hash_func> compute_histogram(const system& parent, const reg_type* mask, const hash_func& h)
    {

        // retrieve number of '1's
        size_t dim = register_utils::count<1>(mask, parent._props.S);

        // create histogram using given hash functor
        hashmap<hash_func> histogram = hashmap<hash_func>(parent.get_max_values(dim), h, register_compare(parent._props.S));

        // buffer for temporary values
        reg_type* buf = new reg_type[parent._props.S];

        // cycle all samples, apply mask and update frequency of occurrence
        for (const auto& s : parent._samples)
        {
            parent._assign_from_mask(buf, s->_data, mask);
            reg temp(buf, parent._props.S);
            histogram[temp]++;
        }
        delete[] buf;

        return histogram;
    }

}
