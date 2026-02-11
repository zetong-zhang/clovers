/**
 * @brief   Utility functions for bioinformatics.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     1.0.0
 * @date        2025-11-30
 * @modified    2026-02-10
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef BIO_UTIL
#define BIO_UTIL

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include "Encoding.hpp"

/*  codon to amino acid map */
extern char trans_tbl[64];

namespace bio_util {
    /**
     * @brief   count GC number in a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @param   len length of the sequence.
     * @return  float GC number.
     */
    double  gc_count(const char *seq, const size_t len);
    /**
     * @brief   calculate GC fraction in a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @param   len length of the sequence.
     * @return  float GC fraction.
     */
    double  gc_fraction(const char *seq, const size_t len);
    /**
     * @brief   get the complement sequence of a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @return  char* complement sequence.
     */
    char *  get_complement(const std::string &seq);
    /**
     * @brief   get ORFs in a nucleotide sequence.
     * 
     * @param   record    scaffold record.
     * @param   starts    start codon array.
     * @param   stops     stop codon array.
     * @param   circ      treat as circular.
     * @param   minlen    minimum length of the ORF.
     * @param   orfs      ORF array to be written.
     * @param   max_alt   maxinum number of alter start codons.
     */
    void    get_orfs(bio::record &record, const str_array &starts, const str_array &stops, 
        const int minlen, const bool circ, bio::orf_array &orfs, int max_alt);
    /**
     * @brief   translate a nucleotide sequence to a protein sequence.
     * 
     * @param   orf    ORF to be translated.
     * @param   prolen length of the protein sequence.
     * @return  char* protein sequence.
     */
    char *  gene2protein(bio::orf &orf, int prolen);
    /**
     * @brief   get the consensus sequence of a set of flanking sequences.
     * 
     * @param   flankings flanking sequence array.
     * @param   start     start position of the consensus.
     * @param   end       end position of the consensus.
     * @return  char* consensus sequence.
     */
    char *  get_consensus(pch_array flankings, int start, int end);
    /**
     * @brief   check if two ORFs overlap.
     * 
     * @param   a ORF a.
     * @param   b ORF b.
     * @param   min_olen minimum overlap length.
     * @return  op_type overlap type.
     */
    op_type check_overprint(const bio::orf &a, const bio::orf &b, double min_olen=0.6);
}

#endif