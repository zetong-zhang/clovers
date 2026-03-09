#include "BioUtil.hpp"
/* complement base table */
static char COMP[] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N',
  /* A    B    C    D    E    F    G    H  */
    'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 
  /* I    J    K    L    M    N    O    P  */
    'N', 'N', 'M', 'N', 'K', 'N', 'N', 'N',
  /* Q    R    S    T    U    V    W    X  */
    'N', 'Y', 'W', 'A', 'A', 'B', 'S', 'N', 
  /* Y    Z */
    'R', 'T', 'N', 'N', 'N', 'N', 'N', 'N',
  /* a    b    c    d    e    f    g    h  */
    't', 'v', 'g', 'h', 'n', 'n', 'c', 'd', 
  /* i    j    k    l    m    n    o    p  */
    'n', 'n', 'm', 'n', 'k', 'n', 'n', 'n',
  /* q    r    s    t    u    v    w    x  */
    'n', 'y', 'w', 'a', 'a', 'b', 's', 'n', 
  /* y    z */
    'r', 't'
};
std::map<char, int> base2off = {
    { 'T', 0 }, { 'C', 1 }, { 'A', 2 }, { 'G', 3 },
    { 't', 0 }, { 'c', 1 }, { 'a', 2 }, { 'g', 3 }
};
/*  codon to amino acid map */
char trans_tbl[64] = {
    'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
    'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
    'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
    'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 
    'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
    'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 
    'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
    'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};

float bio_util::gc_count(
    const char *pstr,
    const size_t len
) {
    float count = 0.0F;
    for (size_t i = 0; i < len; i ++) {
        char c = std::tolower(pstr[i]);
        count += (int)(c == 'g' || c == 'c');
    }
    return count;
}

float bio_util::gc_fraction(
    const char *pstr,
    const size_t len
) { 
    return gc_count(pstr, len) / len; 
}

char *bio_util::get_complement(
    const std::string &origin
) {
    size_t len = origin.size();
    char *comp = NEW char[len+1];
    if (!comp) return nullptr;
    for (size_t i = 0; i < len; i ++)
        comp[i] = COMP[origin.at(len - i - 1)];
    comp[len] = '\0';
    return comp;
}

int bio_util::match_codon(
    const char *origin,
    const str_array &codons
) {
    int num_codons = (int) codons.size();
    for (int index = 0; index < num_codons; index ++) {
        const char *pstr = codons.at(index).c_str();
        if (!std::strncmp(pstr, origin, 3))
            return (int) index;
    }
    return -1;
}

char *bio_util::gene2protein(bio::orf &gene, int prolen) {
    char *protein = NEW char[prolen+1]; 
    if (!protein) return nullptr;
    try {
        for (int i = 0; i < prolen; i ++) {
            if (!(gene.partial5 || i)) {
                protein[0] = 'M';
                continue;
            }
            int f = base2off.at(gene.seq[i*3+0]);
            int s = base2off.at(gene.seq[i*3+1]);
            int t = base2off.at(gene.seq[i*3+2]);
            protein[i] = trans_tbl[f*16+s*4+t];
        }
    } catch (const std::out_of_range& e) {
        delete[] protein;
        return nullptr;
    }
    protein[prolen] = '\0';
    return protein;
}

static int refine_start(int_array &start_locs, int_array &start_types) {
    int t_start = start_locs.at(0);
    for (int k = 0; k < start_locs.size(); k ++) {
        int alter = start_locs.at(k);
        if (start_types[k] == 0 && (alter-t_start)<=75) {
            t_start = alter;
            break;
        }
    }
    return t_start;
}

op_type bio_util::check_overprint(const bio::orf &a, const bio::orf &b, float ratio) {
    if (!strcmp(a.host, b.host) && a.strand == b.strand) {
        int a_end = a.end, b_end = b.end;
        // 我是德布罗意的六世徒孙
        if (a_end < a.t_start) a_end += a.host_len;
        if (b_end < b.t_start) b_end += b.host_len;
        int olen = 0;
        if (a_end > b.t_start && b_end > a.t_start)
            olen = std::min(a_end, b_end) - std::max(a.t_start, b.t_start);
        if (olen == a.len || olen == b.len) 
            return op_type::INCLUDE;
        else if ((olen / (float)std::min(a.len, b.len)) >= ratio) 
            return op_type::INTERSECT;
    }
    return op_type::DISJOINT;
}
/**
 * @brief   splice the edge of an ORF.
 * 
 * @param orf3  The ORF to be spliced.
 * @param length    The length of the scaffold.
 * @param genom     The genome sequence.
 */
static void splice_edge(bio::orf &orf3, int length, char *genom) {
    int lstart = orf3.starts.at(0);
    int cache_len = length - lstart + orf3.end;
    orf3.seq = NEW char[cache_len+1];
    std::strcpy(orf3.seq, genom + lstart);
    std::strncat(orf3.seq, genom, orf3.end);
    orf3.seq[cache_len] = '\0';
    orf3.seq += orf3.t_start - lstart;
}

/* Error handling */
void bio_util::get_orfs(
    bio::record &scaffold,
    const str_array &starts,
    const str_array &stops,
    const int minlen,
    const bool circ,
    bio::orf_array &orfs
) {
    char *genome = (char *) scaffold.sequence.c_str();
    int length = (int) scaffold.sequence.size();
    char *hostname = (char*)scaffold.name.c_str();
    for (int neg_strand = 0; neg_strand < 2; neg_strand ++) {
        if (neg_strand) {
            genome = get_complement(scaffold.sequence);
            scaffold.complement = genome;
        }
        std::vector<int> pf_orfs(3);
        std::vector<int> fd_orfs(0); 
        for (int phase = 0; phase < 3; phase ++) {
            int_array slocs(0), types(0);
            int ps;
            pf_orfs.at(phase) = orfs.size();

            // search for standard ORFs
            for (ps = phase; ps < length; ps += 3) {
                int match = match_codon(genome+ps, starts);
                if (match > -1) {
                    slocs.push_back(ps);
                    types.push_back(match);
                } else if (slocs.size() && match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3;
                    int t_start = refine_start(slocs, types);
                    int seqlen = end - t_start;
                    if (seqlen < minlen) {
                        t_start = slocs.at(0);
                        seqlen = end - t_start;
                    }
                    if (seqlen >= minlen) {
                        char *pstr = genome + t_start;
                        float gc_frac = gc_fraction(pstr, seqlen);
                        char strand = neg_strand ? '-' : '+';
                        orfs.emplace_back(
                            hostname, length, std::move(slocs), std::move(types), 
                            end, t_start, seqlen, strand, pstr, gc_frac
                        );
                    }
                    slocs.clear(), types.clear();
                }
            }

            if (pf_orfs.at(phase) == orfs.size()) {
                pf_orfs.at(phase) = -1;
            }
            
            bool partial3 = (bool) slocs.size();
            // search for partial 3'-end ORFs
            if (partial3) {
                int end = ps; if (end > length) end = length;
                int t_start = refine_start(slocs, types);
                int seqlen = end - t_start;
                if (seqlen < minlen) {
                    t_start = slocs.at(0);
                    seqlen = end - t_start;
                }
                if (seqlen >= minlen || circ) {
                    char *pstr = genome + t_start;
                    float gc_frac = gc_fraction(pstr, seqlen);
                    char strand = neg_strand ? '-' : '+';
                    seqlen = seqlen / 3 * 3;
                    orfs.emplace_back(
                        hostname, length, std::move(slocs), std::move(types), 
                        end, t_start, seqlen, strand, pstr, gc_frac
                    );
                    orfs.at(orfs.size()-1).partial3 = true;
                }
            }

            // search for partial 5'-end ORFs
            int start = -1;
            int ps0 = (phase + ((3 - length % 3) % 3)) % 3;
            for (int ps = ps0; ps < length; ps += 3) {
                if (start == -1 && match_codon(genome+ps, starts) > -1) {
                    start = ps;
                    continue;
                }
                if (match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3;
                    if (circ && partial3) {
                        bio::orf &orf3 = orfs.at(orfs.size()-1);
                        int t_start = refine_start(orf3.starts, orf3.types);
                        int seqlen = length - orf3.t_start + end;
                        if (seqlen < minlen) {
                            t_start = orf3.starts.at(0);
                            seqlen = length - t_start + end;
                        }
                        if (seqlen >= minlen) {
                            orf3.pstr = genome + t_start;
                            orf3.t_start = t_start;
                            orf3.end = end;
                            orf3.len = seqlen;
                            splice_edge(orf3, length, genome);
                            orf3.gc_frac = gc_fraction(orf3.seq, seqlen);
                            orf3.partial3 = false;
                            if (start >= 0 && (end - start) >= minlen) {
                                fd_orfs.push_back(ps0);
                            }
                        } else orfs.pop_back();
                    } else {
                        int t_start = ps0;
                        int seqlen = end - t_start;
                        if (start < 0 && seqlen >= minlen) {
                            char *pstr = genome + t_start;
                            float gc_frac = gc_fraction(pstr, seqlen);
                            char strand = neg_strand ? '-' : '+';
                            orfs.emplace_back(
                                hostname, length, int_array(0), int_array(0), 
                                end, t_start, seqlen, strand, pstr, gc_frac
                            );
                            orfs.at(orfs.size()-1).partial5 = true;
                        }
                    }
                    break;
                }
            }
        }

        std::vector<size_t> to_del;
        for (auto index : fd_orfs) {
            if (pf_orfs.at(index) >= 0) {
                to_del.push_back(pf_orfs.at(index));
            }
        }
        std::sort(to_del.begin(), to_del.end(), std::greater<size_t>());
        for (auto idx : to_del) orfs.erase(orfs.begin() + idx);
    }
}