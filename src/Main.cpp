/*** * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
 *    PROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS    *
 *                                                         *    
 *   @copyright: (C) 2003-2026 TUBIC, Tianjin University   *
 *   @author:    Zetong Zhang, Yan Lin, Feng Gao           *
 *   @version:   1.0.0                                     *
 *   @date:      2025-11-30                                *
 *   @license:   GNU GPLv3                                 *
 *   @contact:   ylin@tju.edu.cn | fgao@tju.edu.cn         *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <omp.h>
#include <chrono>
#include <sys/stat.h>

#include "cxxopts.hpp"
// #include <thread>
#include "BioIO.hpp"

/* run the program in quiet mode */
static bool      QUIET = false;
/* default start codon types     */
static str_array STARTS { "ATG", "GTG", "TTG" };
/* default stop codon types      */
static str_array STOPS  { "TAA", "TAG", "TGA" };

int main(int argc, char *argv[]) {
    std::srand(32522);
    auto start_t = std::chrono::system_clock::now();  // program start time
    std::string exe_name;  // executable file name
    {
        std::string argv0(argv[0]);
        size_t pos = argv0.find_last_of("/\\");
        if (pos == std::string::npos) exe_name = argv0;
        else exe_name = argv0.substr(pos + 1);
    }
    
    /* set and parse parameters */
    cxxopts::Options options(exe_name);

    /* general parameters */
    options.add_options("General")
        ("h,help",     "Print help menu and exit.")

        ("q,quiet",    "Run quietly with no stderr output.")

        ("T,threads",  "Number of threads to use. (default: all)",
         cxxopts::value<uint32_t>());
    
    /* input/output parameters */
    options.add_options("Input/Output")
        ("i,input",    "Specify FASTA/Genbank input file or their compressed versions (gzip). (default: stdin)",
         cxxopts::value<std::string>())
        
        ("o,output",   "Specify output file. (default: stdout)",
         cxxopts::value<std::string>())
        
        ("f,format",   "Select output format (gff, gbk, med).",
         cxxopts::value<std::string>()->default_value("gff"))
        
        ("a,faa",      "Write protein translations to the selected file.",
         cxxopts::value<std::string>())
 
        ("d,fna",      "Write nucleotide sequences of genes to the selected file.",
         cxxopts::value<std::string>());
    
    /* clovers parameters */
    options.add_options("CLOVERS")
        ("g,table",    "Specify a translation table to use.",
         cxxopts::value<uint32_t>()->default_value("11"))

        ("l,minlen",   "Specify the mininum length of ORFs.",
         cxxopts::value<uint32_t>()->default_value("90"))

        ("c,circ",     "Treat topology as circular.")
    
        ("p,proc",     "Select procedure (single or meta).",
         cxxopts::value<std::string>()->default_value("single"))

        ("s,thres",    "Specify putative gene score threshold.",
         cxxopts::value<float>()->default_value("0.499"))
        
        ("t,train",    "Write (if none exists) or use the specified training file.",
         cxxopts::value<std::string>());
    
    /* tri-tisa+ parameters */
    options.add_options("TriTISA+")
        ("n,bypass",   "Bypass TriTISA+ and output longest ORFs.")

        ("M,maxiter",  "Max iteration times for TIS revision.",
         cxxopts::value<uint32_t>()->default_value("20"))
        
        ("I,tis",      "Write (if none exists) or use the specified TIS model file.",
         cxxopts::value<std::string>());
    
    /* gol-reporter parameters */
    options.add_options("GOP-Reporter")
        ("O,overprint","Write overprinted genes to the selected file (GFF3 format).",
         cxxopts::value<std::string>())

        ("A,amino",    "Write protein translations of overprinted genes to the selected file.",
         cxxopts::value<std::string>())

        ("D,nucl",     "Write nucleotide sequences of overprinted genes to the selected file.",
         cxxopts::value<std::string>())
        
        ("R,ratio",    "Specify minimum overlap ratio for overprinted genes.",
         cxxopts::value<float>()->default_value("0.6"));
    
    cxxopts::ParseResult args;
    try {
        args = options.parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nUse '-h' or '--help' to show help info\n";
        return 1;
    }

    /* show help information and exit with code 0 */
    if (argc <= 1 || args.count("help")) {
        std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
                  << "PROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS 1.0.1\n\n"
                  << "Copyright:  (C) 2003-2026 TUBIC,Tianjin University     \n"
                  << "Authors:    Zetong Zhang, Yan Lin*, Feng Gao*          \n"
                  << "Date:       November 30, 2025                          \n"
                  << "Contact:    ylin@tju.edu.cn | fgao@tju.edu.cn          \n"
                  << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -"
                  << options.help() << "\nExample: " << exe_name << ' '
                  << "-i example.fa -o example.gff -c -f gff               \n";
        if (argc <= 1) {
            std::cerr << "\nPress Enter or Ctrl+C to exit the help ... ";
            std::cin.get();
        }
        return 0;
    }

    /* set using quiet mode or not */
    if (args.count("quiet")) QUIET = true;
    if (!QUIET) std::cerr << "\nPROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS\n\n";
    
    /* convert time to timestamp for the calculation of used seconds */
    auto start_time_t = std::chrono::system_clock::to_time_t(start_t);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&start_time_t), "%Y-%m-%d %H:%M:%S");
    if (!QUIET) std::cerr << "Program Start Time: " << std::setw(29) << oss.str() << '\n';

    /* set the number of cpus for multi-thread process (openmp) */
    auto num_threads = args.count("threads") ? 
    args["threads"].as<uint32_t>() : omp_get_max_threads();
    omp_set_num_threads(num_threads);
    if (!QUIET) std::cerr << "Threads Used:" << std::setw(36) << num_threads << '\n';
    
    /* initialize translation table and start/stop codon list */
    uint32_t table = args["table"].as<uint32_t>();
    if (table == 1) {
        STARTS.assign({ "ATG "});
    } else if (table == 4) {
        STOPS.assign({ "TAA", "TAG" });
        trans_tbl[14] = 'W';
    } else if (table != 11) {
        std::cerr << "\nError: unsupported translation table " << table 
                  << " (options: 1 (standard), 4 (mycoplasma, spiroplasma), 11 (bacteria, archaea))\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Translation Table:" << std::setw(31) << table << '\n';

    /* set mininum gene length */
    auto minlen = args["minlen"].as<uint32_t>();
    if (minlen < 30) {
        std::cerr << "\nError: mininum gene length shorter than 30\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Mininum Gene Length:" << std::setw(26) << minlen << " nt\n";

    /* check output format */
    auto format = args["format"].as<std::string>();
    if (format != "gff" && format != "gbk" && format != "med") {
        std::cerr << "\nError: unknown output format (options: gff (GFF3), gbk (GenBank), med (MED))\n";
        return 1;
    }

    /* handle input */
    std::string input = "-";
    if (args.count("input")) input = args["input"].as<std::string>();
    bio::record_array scaffolds(0);
    if (!bio_io::read_source(input, scaffolds)) return 1;
    if (!scaffolds.size()) {
        std::cerr << "\nError: no valid sequence was read\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Number of Scaffolds: " << std::setw(28) << scaffolds.size() << '\n';

    size_t total_len = 0;
    float gc_cont = 0.0;

    /* extract all the orfs */
    bool circ = (bool) args.count("circ");
    bio::orf_array orfs;
    for (int i = 0; i < scaffolds.size(); i ++) {
        const auto sublen = scaffolds[i].sequence.length();
        total_len += sublen;
        if (total_len >= 100000000) {
            std::cerr << "\nError: sequence exceeds 100000000 bp and cannot be processed.\n";
            return 1;
        }
        gc_cont += bio_util::gc_count(scaffolds[i].sequence.c_str(), sublen);

        if (sublen >= minlen) {
            bio_util::get_orfs(scaffolds[i], STARTS, STOPS, minlen, circ, orfs);
            if (orfs.size() >= 10000000) {
                std::cerr << "\nError: total ORF exceeds 10000000 and cannot be processed.\n";
                return 1;
            }
        } else if (!QUIET) {
            std::cerr << "Warning: skip " << scaffolds[i].name << " (too short)\n";
        }
    }
    const uint32_t n_orfs = (uint32_t) orfs.size();
    if (n_orfs <= 0) {
        std::cerr << "\nError: no ORF matching the definition was found.\n";
        return 1;
    }

    /* calculate basic genome stats */
    gc_cont /= total_len;
    if (!QUIET) {
        std::cerr << "Genome Size: " << std::setw(33) << total_len << " bp\n";
        std::cerr << "G+C Content: " << std::setw(34) << std::fixed << std::setprecision(2) 
                  << gc_cont * 100 << " % \n";
        std::cerr << "Number of ORFs: " << std::setw(33) << n_orfs << '\n';
    }

    /* sort orfs by G+C content */
    std::sort(orfs.begin(), orfs.end(), [](bio::orf& a, bio::orf& b) { return a.gc_frac < b.gc_frac; });
    int gc_intv_count[N_MODELS+1] = {0};
    {
        float max_gc = 0.20;
        for (int i = 0, j = 0; i < n_orfs; i ++) {
            if (orfs.at(i).gc_frac > max_gc) {
                max_gc += 0.01;
                if ((++j) >= (N_MODELS+1)) break;
            }
            gc_intv_count[j] ++;
        }
    }

    /* convert orfs to zcurve params */
    float *params = NEW float[DIM_S*n_orfs];
    if (params == nullptr) {
        std::cerr << MEM_ERR_INFO;
        return 1;
    }
    encoding::encode_orfs(orfs, params, 3);

    /* load pre-trained models and calculate scores */
    if (!model::init_models()) return 1;
    if(!QUIET) std::cerr << "Loaded Model:" << std::setw(37) << "Prokaryota\n";
    bool training = args["proc"].as<std::string>() != "meta";
    float *i_scores = NEW float[n_orfs]();
    if (i_scores == nullptr) {
        std::cerr << MEM_ERR_INFO;
        return 1;
    }
    int off = 0;
    if (!QUIET) std::cerr << "Initialization:" << std::setw(34) << "0 %";
    for (int i = 0; i < N_MODELS; i ++) {
        int size = gc_intv_count[i];
        model::mlp_predict(i, params+off*DIM_S, size, i_scores+off);
        if (!QUIET) std::cerr << "\rInitialization:" << std::setw(32) << (int)(i*1.67) << " %";
        off += gc_intv_count[i];
    }

    /* train rbf-svm model */
    if (!QUIET) std::cerr << "\nTraining CDS Model ...";
    float *d_scores = NEW float[n_orfs]();
    if (d_scores == nullptr) {
        std::cerr << MEM_ERR_INFO;
        return 1;
    }
    float mins[DIM_S], maxs[DIM_S];
    svm_model* cds_model = nullptr;
    if (training) {
        std::string model_file;
        if (args.count("train")) {
            model_file = args["train"].as<std::string>();
            if (bio_io::file_exists(model_file)) {
                cds_model = NEW svm_model();
                if (!cds_model) {
                    std::cerr << "\n\n" << MEM_ERR_INFO;
                    return 1;
                }
                if (!bio_io::read_model(cds_model, mins, maxs, model_file)) {
                    delete cds_model;
                    cds_model = nullptr;
                } else {
                    encoding::minmax_scale(params, n_orfs, DIM_S, mins, maxs);
                }
            }
        }
        if (cds_model == nullptr) {
            cds_model = model::rbf_train(params, n_orfs, DIM_S, i_scores, mins, maxs);
        }
        // predict scores
        if (cds_model) {
            #ifdef _OPENMP
                #pragma omp parallel for
            #endif
            for (int i = 0; i < n_orfs; i ++) {
                d_scores[i] = svm_predict_score(cds_model, params + i*DIM_S, DIM_S);
            }
            if (!model_file.empty()) {
                bio_io::write_model(cds_model, mins, maxs, model_file);
            }
        }
        if (!QUIET) std::cerr << std::setw(28) << (cds_model ? "Done\n" : "Skipped\n");
    } else if (!QUIET) std::cerr << std::setw(28) << "Bypassed\n";

    /* classifying orfs and selection of seed orfs */
    auto thres = args["thres"].as<float>();
    bio::orf_array putative;
    for (int i = 0, j = 0; i < n_orfs; i ++) {
        orfs[i].score = std::max(i_scores[i]*W, d_scores[i]);
        if (orfs[i].score > thres) putative.push_back(orfs[i]);
    }

    /* revise gene starts*/
    bool flag = !(bool) args.count("bypass");
    if (flag) {
        if (!QUIET) std::cerr << "Revising TIS Model ..." << std::setw(27) << "Round #0";
        int n_seeds = 0, max_alter;
        for (int i = 0; i < putative.size(); i ++) 
            if (putative[i].len >= 300) n_seeds ++;
        float *params = nullptr, pFU, pFD;
        auto maxiter = args["maxiter"].as<uint32_t>();
        std::string model_file;
        if (args.count("tis")) {
            model_file = args["tis"].as<std::string>();
            if (bio_io::file_exists(model_file)) {
                params = NEW float[TIS_S];
                if (!bio_io::read_model(params, max_alter, pFU, pFD, model_file)) {
                    delete[] params;
                    params = nullptr;
                }
            }
        }
        if (params == nullptr && n_seeds >= MIN_MARKOV_SET) {
            params = NEW float[TIS_S];
            if (params == nullptr) {
                std::cerr << MEM_ERR_INFO << '\n';
                return 1;
            }
            int round = 0;
            for (int order = 1; order < 3; order ++) {
                for (int iter = 0; iter < maxiter; iter ++) {
                    round += 1;
                    if (!QUIET) std::cerr << "\rRevising TIS Model ..." << std::setw(27) 
                                          << "Round #" + std::to_string(round);
                    model::mm_train(putative, order, params, STARTS, table, pFU, pFD, max_alter);
                    float ratio = model::mm_revise(putative, order, params, pFU, pFD, max_alter);
                    if (ratio > 0.9) break;
                }
            }
        } else if (params != nullptr) {
            model::mm_revise(putative, 2, params, pFU, pFD, max_alter);
        } else {
            model::mm_revise(putative, 2, tis_params, 0.855970, 5.144030, 7);
            if (!QUIET) std::cerr << "\rRevising TIS Model ..." << std::setw(27) << "Skipped";
        }
        if (params && !model_file.empty()) {
            bio_io::write_model(params, max_alter, pFU, pFD, model_file);
        }
        delete[] params;
    } else if (!QUIET) std::cerr << "Revising TIS Model ..." << std::setw(27) << "Bybassed";

    if (cds_model) svm_free_model_content(cds_model);
    int num_putative = (int) putative.size();
    if (!QUIET) std::cerr << "\nNumber of Putative Genes:" << std::setw(24) << num_putative << "\n";

    delete[] d_scores;
    delete[] i_scores;
    delete[] params;

    /* write results to output */
    std::string output = "-";
    if (args.count("output")) output = args["output"].as<std::string>();
    std::sort(putative.begin(), putative.end(), [](bio::orf& a, bio::orf& b) {
        int host_cmp = std::strcmp(a.host, b.host);
        if (host_cmp < 0) return true;
        else if (host_cmp > 0) return false;
        else {
            int a_end = (a.strand=='+')?a.end:a.host_len-a.t_start;
            int b_end = (b.strand=='+')?b.end:b.host_len-b.t_start;
            return a_end < b_end;
        }
    });
    bio_io::write_result(putative, output, format);

    /* write protein sequences */
    if (args.count("faa")) {
        auto faa = args["faa"].as<std::string>();
        if(!bio_io::write_faa(putative, faa)) return 1;
    }

    /* write nucleotide sequences */
    if (args.count("fna")) {
        auto fna = args["fna"].as<std::string>();
        if(!bio_io::write_fna(putative, fna)) return 1;
    }

    /* check overprinted genes */
    float min_ratio = args["ratio"].as<float>();
    if (min_ratio < 0.0 || min_ratio > 1.0) {
        std::cerr << "\nError: invalid minimum overlap ratio (range: 0.0-1.0)\n";
        return 1;
    }
    bio::orf_array op_genes;
    for (int i = 0; i < num_putative; i ++) {
        for (int j = i + 1; j < num_putative; j ++) {
            bio::orf &gene_1 = putative[i];
            bio::orf &gene_2 = putative[j];
            op_type type = bio_util::check_overprint(gene_1, gene_2, min_ratio);
            if (type != op_type::DISJOINT) {
                op_genes.push_back(gene_1);
                op_genes.push_back(gene_2);
            }
        }
    }

    /* write overprinted gene coords */
    if (args.count("overprint")) {
        auto overprint = args["overprint"].as<std::string>();
        if(!bio_io::write_overprint(op_genes, overprint)) return 1;
    }

    /* write overprinted gene proteins */
    if (args.count("amino")) {
        auto amino = args["amino"].as<std::string>();
        if(!bio_io::write_faa(op_genes, amino)) return 1;
    }

    /* write overprinted gene nucleotide sequences */
    if (args.count("nucl")) {
        auto nucl = args["nucl"].as<std::string>();
        if(!bio_io::write_fna(op_genes, nucl)) return 1;
    }

    if (!QUIET) {
        auto end_t = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t);
        auto seconds = duration.count() / 1000.0;
        std::cerr << "\nFinished in " << std::fixed << std::setprecision(3) << seconds << " s\n";
    }
}