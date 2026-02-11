/** 
 * @brief       Bioinformatics I/O operations. 
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     1.0.0
 * @date        2025-11-30
 * @modified    2026-02-10
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */

#ifndef BIO_IO
#define BIO_IO
/* optional: zlib */
#ifdef ZLIB
#include <zlib.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <assert.h>
#include <iomanip>
#include "Model.hpp"

#include "BioUtil.hpp"
/* The buffer size for file reading. */
#define BUFF_SIZE 65536
/* The version info of the software. */
#define VERSION  "CLOVERS_v1.0.0"
/* The namespace for bioinformatics I/O operations. */
namespace bio_io {
    /**
     * @brief   Check if a file exists.
     * 
     * @param path  The path to the file.
     * @return true If the file exists.
     *         false If the file does not exist.
     */
    bool file_exists(std::string &path);
    /** 
     * @brief   Read a FASTA or GenBank file from a stream.
     * 
     * @param filename  The path to the file.
     * @param scaffolds The record array to store the scaffolds.
     * @return true If the file is successfully read.
     *         false If the file is not successfully read.
     */
    bool read_source(const std::string &filename, bio::record_array &scaffolds);
    /** 
     * @brief   Write the ORFs to a result file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @param format    The format of the output file.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_result(bio::orf_array &orfs, const std::string &filename, const std::string &format);
    /** 
     * @brief   Write the ORFs to a FASTA file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_faa(bio::orf_array &orfs, const std::string &filename);
    /** 
     * @brief   Write the ORFs to a FASTA file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_fna(bio::orf_array &orfs, const std::string &filename);
    /** 
     * @brief   Write overprinted genes to a GFF3 file.
     * 
     * @param orfs      The overprinted genes.
     * @param filename  The path to the output file.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_overprint(bio::orf_array &orfs, const std::string &filename);
    /** 
     * @brief   Write the svm model to a file.
     * 
     * @param model     The svm model.
     * @param mins      The minimum values of the scaler.
     * @param maxs      The maximum values of the scaler.
     * @param filename  The path to the output file.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_model(svm_model* model, double *mins, double *maxs, const std::string &filename);
    /** 
     * @brief   Read the svm model from a file.
     * 
     * @param model     The svm model.
     * @param mins      The minimum values of the scaler.
     * @param maxs      The maximum values of the scaler.
     * @param filename  The path to the input file.
     * @return true If the file is successfully read.
     *         false If the file is not successfully read.
     */
    bool read_model(svm_model *model, double *mins, double *maxs,const std::string &filename);
}

#endif