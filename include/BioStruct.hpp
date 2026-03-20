/** 
 * @brief       Bioinformatics data structures. 
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     1.0.3
 * @date        2025-11-30
 * @modified    2026-02-28
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef BIO_STRUCT
#define BIO_STRUCT
/* no throw memory error new */
#define NEW new (std::nothrow)
/* memory error info */
#define MEM_ERR_INFO "\nError: memory error (cannot alloc memory)\n"
/* epsilon value for floating point comparisons. */
#define EPSILON 1E-12

#include <string>
#include <vector>

/* utility array types. */
typedef std::vector<int>          int_array;
typedef std::vector<float>        flt_array;
typedef std::vector<std::string>  str_array;
typedef std::vector<char *>       pch_array;

/* overprinting types */
enum op_type { DISJOINT, INTERSECT, INCLUDE };

/* bioinformatics data structures */
namespace bio {
    /*  scaffold record type */
    typedef struct record {
        std::string name;            // scaffold name
        std::string sequence;        // scaffold sequence
        char *complement = nullptr;  // complement sequence
        record() noexcept: name("anonymous") {}
        record(std::string &&name, std::string &&sequence) noexcept:
        name(std::move(name)), sequence(std::move(sequence)) {}
        ~record() noexcept { delete[] complement; }
    } record;
    /*  ORF type */
    typedef struct orf {
        char *     host;            // host scaffold name
        int        host_len;        // length of host
        int_array  starts;          // alter start positions
        int_array  types;           // alter start types 
        int_array  scores;          // scores for each start
        int        t_start;         // true start position
        int        end;             // true end position
        int        len;             // length
        char       strand;          // strand direction
        char *     pstr;            // pointer on genome
        char *     seq;             // nucleotide sequence
        float      gc_frac;         // GC fraction
        float      score=0;         // zcurve score
        bool       partial5=false;  // partial 5'-end
        bool       partial3=false;  // partial 3'-end
        orf() noexcept: host((char*)"anonymous") {}
        orf(char *host, int host_len, int_array &&starts, int_array &&types, int end, 
        int t_start, int len, char strand, char *pstr, float gc_frac) noexcept: seq(pstr),
        host(host), host_len(host_len),starts(std::move(starts)),types(std::move(types)), 
        len(len), t_start(t_start),end(end),strand(strand),pstr(pstr),gc_frac(gc_frac){}
    } orf;
    /*  scaffold record array type */
    typedef std::vector<record>  record_array;
    /*  ORF array type */
    typedef std::vector<orf>     orf_array;
}

#endif