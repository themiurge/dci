#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <cmath>
#include <random>
#include <initializer_list>
#include <algorithm>
#include <limits>
#include <queue>

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
    constexpr fp_type eps = 1e-6;
    constexpr fp_type nan = numeric_limits<fp_type>::quiet_NaN();

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
        inline reg_type get_last_reg_mask(const size_t& n) { return ((reg_type)1 << (n % bits_per_reg)) - 1; }

        //
        // template, template specializations and general function for negation
        //
        inline void negate(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask, const size_t& S)
        {
            for (size_t i = 0; i != S; ++i) dest[i] = ~src[i]; dest[S-1] &= last_reg_mask;
        }
        template<size_t S> void negate(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask)
        {
            return negate(dest, src, last_reg_mask, S);
        }
        template<> inline void negate<1>(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask) { dest[0] = ~src[0] & last_reg_mask; }
        template<> inline void negate<2>(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask) { dest[0] = ~src[0]; dest[1] = ~src[1] & last_reg_mask; }
        template<> inline void negate<3>(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask) { dest[0] = ~src[0]; dest[1] = ~src[1]; dest[2] = ~src[2] & last_reg_mask; }
        template<> inline void negate<4>(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask) { dest[0] = ~src[0]; dest[1] = ~src[1]; dest[2] = ~src[2]; dest[3] = ~src[3] & last_reg_mask; }

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
        template<> inline bool equal<3>(const reg_type* a, const reg_type* b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; }
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
    
    //
    // This class is used to compute statistics
    // over a running sequence of fp_type values
    //
    class running_stat 
    {
    public:
        running_stat() : m_n(0) {}

        void clear()
        {
            m_n = 0;
        }

        void push(fp_type x)
        {
            m_n++;

            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (m_n == 1)
            {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
            }
            else
            {
                m_newM = m_oldM + (x - m_oldM)/m_n;
                m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

                // set up for next iteration
                m_oldM = m_newM; 
                m_oldS = m_newS;
            }
        }

        int num_data_values() const
        {
            return m_n;
        }

        fp_type mean() const
        {
            return (m_n > 0) ? m_newM : 0.0;
        }

        fp_type variance() const
        {
            return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
        }

        fp_type standard_deviation() const
        {
            return sqrt( variance() );
        }

    private:

        size_t m_n;
        fp_type m_oldM, m_newM, m_oldS, m_newS;
    };

    // Class for register negation
    class register_negate {
    private:
        size_t _s;
        reg_type _last_reg_mask;
        void (*_f)(reg_type*, const reg_type*, const reg_type&, const size_t&);
        template<size_t S> inline static void _wrapper(reg_type* dest, const reg_type* src, const reg_type& last_reg_mask, const size_t& s) { register_utils::negate<S>(dest, src, last_reg_mask); }
    public:
        register_negate(const size_t& s = 0, const size_t& n = 0) : _s(s)
        { 
            // use template specializations for small values of S
            switch (s)
            {
                case 1: _f = _wrapper<1>; break;
                case 2: _f = _wrapper<2>; break;
                case 3: _f = _wrapper<3>; break;
                case 4: _f = _wrapper<4>; break;
                default: _f = register_utils::negate; break;
            }

            // build last reg mask
            _last_reg_mask = register_utils::get_last_reg_mask(n);
        }
        inline void operator()(reg_type* dest, const reg_type* src) const { _f(dest, src, _last_reg_mask, _s); }
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
        size_t SA; // number of registers per cluster (compact format, i.e. one bit per agent)
        size_t NB; // number of bits per sample
        size_t sample_bytes;
        size_t data_bytes;
        size_t agent_pool_bytes;
        size_t compact_cluster_bytes;
        inline void init(const size_t& N, const size_t& M, const size_t& NB)
        {
            this->N = N; this->M = M; this->NB = NB;
            S = register_utils::number_of_registers(NB);
            SA = register_utils::number_of_registers(N);
            compact_cluster_bytes = SA * bytes_per_reg;
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

    // This class represents a cluster generator (by lexicographic permutation)
    class cluster_generator;

    // This is just a typedef for a categorical distribution
    template<class hash_func> using hashmap = unordered_map<reg, size_t, hash_func, register_compare>;
    using categorical_distribution = hashmap<murmur3_hasher>;
    using histogram_type = hashmap<basic_hasher>;

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
        fp_type _entropy;
        bool _got_entropy;

        categorical_distribution _pdf;
        agent() : _bitmask(nullptr), _name(), _pdf(), _id(0), _got_entropy(false) { }
        agent(reg_type* agent_bitmask, const string& agent_name, const size_t& id, system* parent_system) : 
            _bitmask(agent_bitmask), _name(agent_name), _parent(parent_system), _pdf(), _id(id), _got_entropy(false)
        { agent::n_new++; }
        void compute_pdf();
    public:

        inline const string& name() const { return _name; }
        inline const size_t& id() const { return _id; }
        inline const categorical_distribution& pdf() const { return _pdf; }
        fp_type entropy();
        size_t size() const;

        // friend classes
        friend class system;
        friend class sample; 
        friend class agent_value;   
        friend class cluster;   

        // friend operators
        friend ostream& operator<<(ostream& out, const agent& a);
        friend istream& operator>>(istream& in, system& sys); 
        friend ostream& operator<<(ostream& out, const cluster& c);

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

        template<typename T> agent_value operator=(T t);
        template<typename T> bool operator==(T t) const;

        // friend classes
        friend class agent;
        friend class sample;
        friend class system;
        friend class cluster;   

    };

    class cluster {
    public:
        static size_t n_new;
        static size_t n_del;
    private:
        system* _parent;
        reg_type* _bitmask;
        reg_type* _compact;
        mutable fp_type _entropy, _comp_entropy, _cluster_index, _statistical_index;
        mutable bool _got_entropy, _got_comp_entropy, _got_cluster_index, _got_statistical_index;
        size_t _size;
        mutable histogram_type _histogram;
        
        void build_bitmask();
        void allocate_data();
        void reset_data();
        inline void reset_fields() { _bitmask = nullptr; _compact = nullptr; _histogram.clear(); _size = 0; _got_entropy = false; }
        inline void free_data() { if (_bitmask) delete[] _bitmask; if (_compact) delete[] _compact; _histogram.clear(); cluster::n_del++; }
        void copy_from(const cluster&);
        void move_from(cluster&&);
        void build_from_compact_bitmask(const reg_type*);
        void build_from_compact_string(const string&);
        template<typename T> void build_from_initializer_list(initializer_list<T>);
    public:

        // default constructor, just never use it
        cluster();

        // system only - empty cluster
        cluster(system*);

        // system + compact bitmask
        cluster(system*, const reg_type*);

        // system + string bitmask
        cluster(system*, const string&);

        // system + initializer list
        template<typename T> cluster(system*, initializer_list<T>);

        // copy constructor
        cluster(const cluster&);

        // move constructor
        cluster(cluster&&);

        // compact bitmask - assumes same parent system
        cluster& operator=(const reg_type*);

        // string bitmask - assumes same parent system
        cluster& operator=(const string&);

        // initializer list - assumes same parent system
        template<typename T> cluster& operator=(initializer_list<T>);

        // copy assignment
        cluster& operator=(const cluster&);

        // move assignment
        cluster& operator=(cluster&&);
        
        // friend classes
        friend class cluster_generator;
        friend class system;

        // friend operators
        friend bool operator==(const cluster&, const cluster&);
        friend bool operator!=(const cluster&, const cluster&);
        friend cluster operator!(const cluster&);
        friend ostream& operator<<(ostream& out, const cluster& c);

        // member access
        inline size_t size() { return _size; }
        const histogram_type& histogram() const;
        fp_type entropy() const;
        fp_type comp_entropy() const;
        fp_type cluster_index() const;
        fp_type statistical_index() const;
        void set_comp_entropy(const fp_type& ent_val) { _comp_entropy = ent_val; _got_comp_entropy = true; }
        
        // comparison operator for priority queues
        inline bool operator<(const cluster& c) const { return statistical_index() > c.statistical_index(); }

        // destructor
        virtual ~cluster() { free_data(); }

    };
    
    class cluster_generator {
    private:
        system* _parent;
        size_t _size;
        string _current;
        bool _has_next;
        cluster_generator() { }
    public:
        cluster_generator(system* parent, size_t size);
        bool has_next() const { return _has_next; }
        cluster next();
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
        
        vector<fp_type> _dci_mean;
        vector<fp_type> _dci_sigma;

        system_properties _props;
        fp_type _system_entropy;
        size_t _system_diff_vals;
        
        bool _got_cluster_index_statistics;

        // register functors - standard (_props.S registers)
        register_compare _equal;
        register_negate _negate;
        register_assign _assign;
        register_assign_from_mask _assign_from_mask;
        register_assign_in_mask _assign_in_mask;
        register_combine _combine;
        register_combine_from_mask _combine_from_mask;

        // register functors - compact (_props.SA registers)
        register_compare _equal_compact;
        register_negate _negate_compact;
        register_assign _assign_compact;
        register_assign_from_mask _assign_from_mask_compact;
        register_assign_in_mask _assign_in_mask_compact;
        register_combine _combine_compact;
        register_combine_from_mask _combine_from_mask_compact;

        // hash functors
        basic_hasher _cluster_hasher;
        murmur3_hasher _agent_hasher;

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
            
            _system_diff_vals = _props.M;
            
            _dci_mean.resize(_props.N - 2, 0.0);
            _dci_sigma.resize(_props.N - 2, 0.0);

            // initialize standard register functors
            _equal = register_compare(_props.S);
            _negate = register_negate(_props.S, _props.NB);
            _assign = register_assign(_props.S);
            _assign_from_mask = register_assign_from_mask(_props.S);
            _assign_in_mask = register_assign_in_mask(_props.S);
            _combine = register_combine(_props.S);
            _combine_from_mask = register_combine_from_mask(_props.S);

            // initialize compact register functors
            _equal_compact = register_compare(_props.SA);
            _negate_compact = register_negate(_props.SA, _props.N);
            _assign_compact = register_assign(_props.SA);
            _assign_from_mask_compact = register_assign_from_mask(_props.SA);
            _assign_in_mask_compact = register_assign_in_mask(_props.SA);
            _combine_compact = register_combine(_props.SA);
            _combine_from_mask_compact = register_combine_from_mask(_props.SA);
            
            // initialize flags
            _got_cluster_index_statistics = false;

            // hash functors
            _cluster_hasher = basic_hasher(_props.S);
            _agent_hasher = murmur3_hasher(_props.S);
            
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
                _agents[a] = new agent(_agent_pool + a * _props.S, to_string(a), a, this);
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
            if (nbits > sizeof(size_t) * 8 || (1 << nbits) > _system_diff_vals) return _system_diff_vals;
            return 1 << nbits;
        }
        inline void check_insert_cluster_in_queue(const cluster& c, priority_queue<cluster>& q, const size_t& num_results)
        {
            if (isinf(c.statistical_index())) return;
            if (q.size() < num_results)
            {
                q.push(c);
                return;
            }
            if (q.top().statistical_index() < c.statistical_index())
            {
                q.pop();
                q.push(c);
            }
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
        friend class cluster_generator;

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
        friend bool operator==(const cluster& c1, const cluster& c2);
        friend bool operator!=(const cluster& c1, const cluster& c2);
        friend cluster operator!(const cluster&);
        friend ostream& operator<<(ostream& out, const cluster& c);

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
        inline fp_type get_mean_dci(const size_t& cluster_size) { return _dci_mean[cluster_size - 2]; }
        inline fp_type get_sigma_dci(const size_t& cluster_size) { return _dci_sigma[cluster_size - 2]; }

        void compute_agent_statistics()
        {
            for (auto& a : _agents) a->compute_pdf();
            cluster sys(this, string(_props.N, '1'));
            _system_entropy = sys.entropy();
            _system_diff_vals = sys.histogram().size();
        }
        
        void compute_cluster_index_statistics()
        {
        
            if (_got_cluster_index_statistics) return;
        
            // create rolling stat objects
            vector<running_stat> stats(_props.N - 2);
            
            // cycle all cluster sizes
            for (size_t s = 1; s != _props.N / 2 + 1; ++s)
            {
                cluster_generator gen(this, s);
                while (gen.has_next())
                {
                    cluster cur = gen.next();
                    cluster comp = !cur;
                    cur.set_comp_entropy(comp.entropy());
                    comp.set_comp_entropy(cur.entropy());
                    if (s > 1 && !isnan(cur.cluster_index())) stats[s - 2].push(cur.cluster_index());
                    if ((_props.N % 2 || s != _props.N / 2) && !isnan(comp.cluster_index())) stats[_props.N - s - 2].push(comp.cluster_index());
                }
                
            }
            
            // copy stats to local arrays
            for (size_t i = 0; i != _props.N - 2; ++i)
            {
                _dci_mean[i] = stats[i].mean();
                _dci_sigma[i] = stats[i].standard_deviation();
            }
            
            _got_cluster_index_statistics = true;
            
        }
        
        void load_cluster_index_statistics_from_homogeneous_system(system& hs)
        {
        
            // compute hs stats if necessary
            hs.compute_cluster_index_statistics();
            
            // copy stats
            _dci_mean = hs._dci_mean;
            _dci_sigma = hs._dci_sigma;
            
        }
        
        ostream& write_cluster_index_statistics(ostream& out)
        {
            out << 0.0 << endl << 0.0 << endl;
            for (size_t i = 0; i != _props.N - 2; ++i)
                out << _dci_mean[i] << endl << _dci_sigma[i] << endl;
            out << 0.0 << endl << 0.0 << endl;
            return out;
        }
        
        //
        // Output computation
        //
        priority_queue<cluster> compute_system_statistics(const size_t& num_results)
        {
        
            priority_queue<cluster> results;
            
            // cycle all cluster sizes
            for (size_t s = 1; s != _props.N / 2 + 1; ++s)
            {
                cluster_generator gen(this, s);
                while (gen.has_next())
                {
                    cluster cur = gen.next();
                    cluster comp = !cur;
                    cur.set_comp_entropy(comp.entropy());
                    comp.set_comp_entropy(cur.entropy());
                    if (s > 1 && !isnan(cur.statistical_index())) check_insert_cluster_in_queue(cur, results, num_results);
                    if ((_props.N % 2 || s != _props.N / 2) && !isnan(comp.statistical_index())) check_insert_cluster_in_queue(comp, results, num_results);
                }
                
            }
            
            return results;
            
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
        _pdf = compute_histogram(*_parent, _bitmask, _parent->_agent_hasher);
    }

    inline fp_type agent::entropy() { if (!_got_entropy) { _entropy = dci::entropy(_pdf, _parent->_props.M); _got_entropy = true; } return _entropy; }

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
    inline void cluster::reset_data()
    {
        if (!_bitmask) return;
        memset(_bitmask, 0, _parent->_props.sample_bytes);
        memset(_compact, 0, _parent->_props.compact_cluster_bytes);
    }

    inline void cluster::allocate_data()
    {
        if (_parent && _parent->_props.S)
        {
            _bitmask = new reg_type[_parent->_props.S];
            _compact = new reg_type[_parent->_props.SA];
            reset_data();
        }
        else
        {
            _bitmask = nullptr;
            _compact = nullptr;
        }
        _got_entropy = false;
        _got_comp_entropy = false;
        _got_cluster_index = false;
        _got_statistical_index = false;
        cluster::n_new++;
    }

    inline void cluster::build_bitmask()
    {
        _size = 0;
        for (size_t a = 0; a != _parent->_props.N; ++a)
            if (register_utils::get_bit(_compact, a))
            {
                ++_size;
                _parent->_combine(_bitmask, _parent->_agents[a]->_bitmask);
            }
    }

    inline void cluster::copy_from(const cluster& c)
    {
        _parent = c._parent; 
        allocate_data(); 
        _parent->_assign_compact(_compact, c._compact);
        _parent->_assign(_bitmask, c._bitmask);

        _got_entropy = c._got_entropy;
        _entropy = c._entropy;
        _got_comp_entropy = c._got_comp_entropy;
        _comp_entropy = c._comp_entropy;
        _got_cluster_index = c._got_cluster_index;
        _cluster_index = c._cluster_index;
        _got_statistical_index = c._got_statistical_index;
        _statistical_index = c._statistical_index;

        _histogram = c._histogram;
        _size = c._size;
    }

    inline void cluster::move_from(cluster&& c)
    {
        _parent = c._parent; 
        _compact = c._compact;
        _bitmask = c._bitmask;
        
        _got_entropy = c._got_entropy;
        _entropy = c._entropy;
        _got_comp_entropy = c._got_comp_entropy;
        _comp_entropy = c._comp_entropy;
        _got_cluster_index = c._got_cluster_index;
        _cluster_index = c._cluster_index;
        _got_statistical_index = c._got_statistical_index;
        _statistical_index = c._statistical_index;

        _size = c._size;
        _histogram = move(c._histogram);

        c.reset_fields();
    }

    inline void cluster::build_from_compact_bitmask(const reg_type* compact)
    {
        _parent->_assign_compact(_compact, compact); build_bitmask();
    }

    inline void cluster::build_from_compact_string(const string& compact)
    {
        register_utils::set_bits_from_string(_compact, compact); build_bitmask();
    }

    template<typename T> inline void cluster::build_from_initializer_list(initializer_list<T> l)
    {
        for (const auto& a : l) register_utils::set_bit(_compact, a, 1); build_bitmask();
    }

    inline const histogram_type& cluster::histogram() const
    {
        if (_histogram.size() == 0) _histogram = compute_histogram(*_parent, _bitmask, _parent->_cluster_hasher);
        return _histogram;
    }

    inline fp_type cluster::entropy() const { if (!_got_entropy) { _entropy = dci::entropy(histogram(), _parent->_props.M); _got_entropy = true; } return _entropy; }
    
    inline fp_type cluster::comp_entropy() const { if (!_got_comp_entropy) { _comp_entropy = (!(*this)).entropy(); _got_comp_entropy = true; } return _comp_entropy; }
    
    inline fp_type cluster::cluster_index() const 
    { 
        if (!_got_cluster_index) 
        { 
            fp_type integration = -entropy();
            for (size_t a = 0; a != _parent->_props.N; ++a)
                if (register_utils::get_bit(_compact, a))
                    integration += _parent->_agents[a]->entropy();
            fp_type mutual_information = entropy() + comp_entropy() - _parent->_system_entropy;
            if (fabs(mutual_information) < eps) 
                _cluster_index = numeric_limits<fp_type>::quiet_NaN();
            else
                _cluster_index = integration / mutual_information;
            //if (_size == 14) cout << "H = " << entropy() << ", Hc = " << comp_entropy() << ", Hs = " << _parent->_system_entropy << ", I = " << integration << ", M = " << mutual_information << ", DCI = " << _cluster_index << endl;
            _got_cluster_index = true; 
        } 
        return _cluster_index; 
    }
    
    inline fp_type cluster::statistical_index() const { if (!_got_statistical_index) { _statistical_index = (cluster_index() - _parent->get_mean_dci(_size)) / _parent->get_sigma_dci(_size); _got_statistical_index = true; } return _statistical_index; }
    
    // default constructor, just never use it
    cluster::cluster() { _parent = nullptr; allocate_data(); }

    // system only - empty cluster
    cluster::cluster(system* s) { _parent = s; allocate_data(); }

    // system + compact bitmask
    cluster::cluster(system* s, const reg_type* compact) { _parent = s; allocate_data(); build_from_compact_bitmask(compact); }

    // system + string bitmask
    cluster::cluster(system* s, const string& compact) { _parent = s; allocate_data(); build_from_compact_string(compact); }

    // system + initializer list
    template<typename T> cluster::cluster(system* s, initializer_list<T> l) { _parent = s; allocate_data(); build_from_initializer_list(l); }

    // copy constructor
    cluster::cluster(const cluster& c) { copy_from(c); }

    // move constructor
    cluster::cluster(cluster&& c) { move_from(move(c)); }

    // compact bitmask - assumes same parent system
    inline cluster& cluster::operator=(const reg_type* compact) { reset_data(); build_from_compact_bitmask(compact); return *this; }

    // string bitmask - assumes same parent system
    inline cluster& cluster::operator=(const string& compact) { reset_data(); build_from_compact_string(compact); return *this; }

    // initializer list - assumes same parent system
    template<typename T> inline cluster& cluster::operator=(initializer_list<T> l) { reset_data(); build_from_initializer_list(l); return *this; }

    // copy assignment
    cluster& cluster::operator=(const cluster&c) { if (&c != this) { free_data(); copy_from(c); } return *this; }

    // move assignment
    cluster& cluster::operator=(cluster&& c) { if (&c != this) { free_data(); move_from(move(c)); } return *this; }

    // comparison
    inline bool operator==(const cluster& c1, const cluster& c2) 
    {
        return c1._parent == c2._parent && (!c1._parent ||
            c1._parent->_equal_compact(c1._compact, c2._compact)
            );
    }

    inline bool operator!=(const cluster& c1, const cluster& c2) 
    {
        return c1._parent != c2._parent || (c1._parent && !c1._parent->_equal_compact(c1._compact, c2._compact));
    }

    inline cluster operator!(const cluster& c)
    {
        cluster comp(c._parent);
        comp._parent->_negate_compact(comp._compact, c._compact);
        comp._parent->_negate(comp._bitmask, c._bitmask);
        comp._size = comp._parent->_props.N - c._size;
        return comp;
    }

    ostream& operator<<(ostream& out, const cluster& c)
    {
        for (size_t i = 0; i != c._parent->_props.N; ++i)
            if (register_utils::get_bit(c._compact, i))
                out << '[' << c._parent->_agents[i]->_name << ']';
        return out;
    }
    
    //
    // cluster_generator class implementation
    //
    cluster_generator::cluster_generator(system* parent, size_t size) : _parent(parent), _size(size), _has_next(true) 
    { 
        _current = string(_parent->_props.N, '0'); for (size_t i = 0; i != _size; ++i) _current[i] = '1'; 
    }
    
    inline cluster cluster_generator::next() 
    { 
        cluster c(_parent, _current); 
        _has_next = prev_permutation(_current.begin(), _current.end()); 
        return c; 
    }

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
    size_t cluster::n_new = 0;
    size_t cluster::n_del = 0;

    ostream& print_allocation_stats(ostream& out)
    {
        out << "reg class performed " << dci::reg::n_new << " new[] and " << dci::reg::n_del << " delete[] statements\n";
        out << "system class performed " << dci::system::n_new << " new[] and " << dci::system::n_del << " delete[] statements\n";
        out << "cluster class performed " << dci::cluster::n_new << " new[] and " << dci::cluster::n_del << " delete[] statements\n";
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
