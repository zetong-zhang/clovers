/*** * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
 *    PROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS    *
 *                                                         *    
 *   @copyright: (C) 2003-2026 TUBIC, Tianjin University   *
 *   @author:    Zetong Zhang, Yan Lin, Feng Gao           *
 *   @version:   0.1.0                                     *
 *   @date:      2025-11-30                                *
 *   @license:   GNU GPLv3                                 *
 *   @contact:   ylin@tju.edu.cn | fgao@tju.edu.cn         *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <omp.h>
#include <chrono>
#include <sys/stat.h>

#include "cxxopts.hpp"
#include "BioIO.hpp"

/* run the program in quiet mode */
static bool   QUIET = false;

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
        ("i,input",    "Specify FASTA/Genbank input file. (default: stdin)",
         cxxopts::value<std::string>())
        
        ("o,output",   "Specify output file. (default: stdout)",
         cxxopts::value<std::string>())
        
        ("f,format",   "Select output format (gff, gbk).",
         cxxopts::value<std::string>()->default_value("gff"))
        
        ("a,faa",      "Write protein translations to the selected file.",
         cxxopts::value<std::string>())
 
        ("d,fna",      "Write nucleotide sequences of genes to the selected file.",
         cxxopts::value<std::string>());
    
    /* cds-finder parameters */
    options.add_options("CDS-Finder")
        ("g,table",    "Specify a translation table to use.",
         cxxopts::value<uint32_t>()->default_value("11"))

        ("l,minlen",   "Specify the mininum length of ORFs.",
         cxxopts::value<uint32_t>()->default_value("90"))

        ("c,circ",     "Treat topology as circular.")
    
        ("b,bypass",   "Bypass semi-supervised SVM training.")

        ("s,thres",    "Specify putative gene score threshold.",
         cxxopts::value<double>()->default_value("0.5"))
        
        ("t,train",    "Write (if none exists) or use the specified training file.",
         cxxopts::value<std::string>());
    
    /* tis-finder parameters */
    options.add_options("TIS-Finder")
        ("L,longest",  "Bypass TIS-Finder and output longest ORFs")

        ("E,alter",    "Maxinum number of alternative start sites in an ORF",
         cxxopts::value<uint32_t>()->default_value("6"))
        
        ("R,rbs",      "Write RBS model to the selected file.",
         cxxopts::value<std::string>());
    
    /* gol-reporter parameters */
    options.add_options("GOP-Reporter")
        ("O,overprint","Write overprinted genes to the selected file (GFF3 format).",
         cxxopts::value<std::string>())

        ("A,amino",    "Write protein translations of overprinted genes to the selected file.",
         cxxopts::value<std::string>())

        ("D,nucl",     "Write nucleotide sequences of overprinted genes to the selected file.",
         cxxopts::value<std::string>())
        
        ("M,ratio",    "Specify minimum overlap ratio for overprinted genes.",
         cxxopts::value<double>()->default_value("0.6"));
    
    cxxopts::ParseResult args;
    try {
        args = options.parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nUse '-h' or '--help' to show help info\n";
        return 1;
    }

    /* show help information and exit with code 0 */
    if (argc <= 1 || args.count("help")) {
        std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - - -\n"
                  << "PROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS\n\n"
                  << "Copyright: (C) 2003-2026 TUBIC,Tianjin University\n"
                  << "Authors:   Zetong Zhang, Yan Lin, Feng Gao       \n"
                  << "Date:      November 30, 2025                     \n"
                  << "Contact:   ylin@tju.edu.cn | fgao@tju.edu.cn     \n"
                  << "- - - - - - - - - - - - - - - - - - - - - - - - -"
                  << options.help() << "\nExample: " << exe_name << ' '
                  << "-i example.fa -o example.gff -c -f gff               \n";
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
    uint32_t table_code = args["table"].as<uint32_t>();
    str_array starts{ "ATG", "GTG" };
    str_array stops{ "TAA", "TAG", "TGA" };

    if (table_code == 1) {
        starts.assign({ "ATG "});
    } else if (table_code == 4) {
        stops.assign({ "TAA", "TAG" });
        trans_tbl[14] = 'W';
    } else if (table_code != 11) {
        std::cerr << "\nError: unsupported translation table " << table_code 
                  << " (options: 1 (standard), 4 (mycoplasma, spiroplasma), 11 (bacteria, archaea))\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Translation Table:" << std::setw(31) << table_code << '\n';

    /* set mininum gene length */
    auto minlen = args["minlen"].as<uint32_t>();
    if (minlen < 30) {
        std::cerr << "\nError: mininum gene length shorter than 30\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Mininum Gene Length:" << std::setw(26) << minlen << " nt\n";

    /* check output format */
    auto format = args["format"].as<std::string>();
    if (format != "gff" && format != "gbk") {
        std::cerr << "\nError: unknown output format (options: gff (GFF3), gbk (GenBank))\n";
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
    double gc_cont = 0.0;

    /* extract all the orfs */
    bool circ = (bool) args.count("circ");
    int max_alt = (int) args["alter"].as<uint32_t>();
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
            bio_util::get_orfs(scaffolds[i], starts, stops, minlen, circ, orfs, max_alt);
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
        double max_gc = 0.20;
        for (int i = 0, j = 0; i < n_orfs; i ++) {
            if (orfs.at(i).gc_frac > max_gc) {
                max_gc += 0.01;
                if ((++j) >= (N_MODELS+1)) break;
            }
            gc_intv_count[j] ++;
        }
    }

    /* convert orfs to zcurve params */
    double *params = NEW double[DIM_S*n_orfs];
    if (params == nullptr) {
        std::cerr << MEM_ERR_INFO;
        return 1;
    }
    encoding::encode_orfs(orfs, params, 3);

    /* load pre-trained models and calculate scores */
    if (!model::init_models()) return 1;
    if(!QUIET) std::cerr << "Loaded Model:" << std::setw(37) << "Prokaryota\n";
    double *i_scores = NEW double[n_orfs]();
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

    /* train rbs-svm model */
    if (!QUIET) std::cerr << "\nTraining CDS Model ...";
    double *d_scores = NEW double[n_orfs]();
    if (d_scores == nullptr) {
        std::cerr << MEM_ERR_INFO;
        return 1;
    }
    bool training = !((bool) args.count("bypass"));
    double mins[DIM_S], maxs[DIM_S];
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
                #pragma omp parallel for schedule(guided)
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
    auto thres = args["thres"].as<double>();
    bio::orf_array putative;
    int_array init_seeds;
    for (int i = 0, j = 0; i < n_orfs; i ++) {
        orfs[i].score = std::max(i_scores[i], d_scores[i]);
        if (orfs[i].score > thres) {
            putative.push_back(orfs[i]);
            j ++;
            if (!orfs[i].starts.size()) continue;
            int up_start = orfs[i].starts.at(0) - minlen;
            int up_end = orfs[i].starts.back() + minlen;
            if (up_start < 0 || up_end > orfs[i].host_len) continue;
            init_seeds.push_back(j-1);
        }
    }
    int n_init_seeds = (int) init_seeds.size();

    if (n_init_seeds > 0) {
        /* plot z-curves for the upstream flanking regions */
        int plot_range = minlen*2;
        pch_array r_flankings, p_flankings, q_flankings;
        int_array s_labels; // labels for start sites
        flt_array exp_lens; // exp length ratios of orfs
        double *sum_curve = NEW double[plot_range*3]();
        if (sum_curve == nullptr) {
            std::cerr << MEM_ERR_INFO;
            return 1;
        }
        int n_alters = 0;
        for (int i = 0; i < n_init_seeds; i ++) {
            bio::orf& pi = putative[init_seeds[i]];
            char *pstr = pi.pstr - minlen;
            r_flankings.push_back(pstr);
            double *curve = NEW double[plot_range*3]();
            if (curve == nullptr) {
                std::cerr << MEM_ERR_INFO;
                return 1;
            }
            encoding::z_curve(pstr, plot_range, curve);
            for (int i = 0; i < plot_range*3; i ++) {
                sum_curve[i] += curve[i];
            }
            delete[] curve;
            int l0 = pi.end - pi.starts.at(0);
            for (int j = 0; j < pi.starts.size(); j ++) {
                int cs = pi.starts[j], ts = pi.t_start;
                char *pstr = pi.pstr + (cs - ts);
                p_flankings.push_back(pstr - minlen);
                q_flankings.push_back(pstr);
                pi.scores.push_back(n_alters); // index of score
                s_labels.push_back(cs==ts ? 1 : -1); // label
                exp_lens.push_back(std::exp(-(double) (pi.end - cs) / l0));
                n_alters ++;
            }
        }
        for (int i = 0; i < plot_range*3; i ++) sum_curve[i] /= (double) n_init_seeds;

        /* detect RBS region and spacer */
        bio::region *rbs_region = nullptr;
        int spacer_len, rbs_len;
        for (int i = 0; i < 3; i ++) {
            bio::region *root = NEW bio::region(0, -1);
            int c = encoding::find_island(sum_curve+i*plot_range, minlen, 6, 0.98, root);
            if (root == nullptr || c < 0) {
                std::cerr << MEM_ERR_INFO;
                return 1;
            }
            rbs_region = root->next;
            for (int i = 0; i < (c - 1); i++) {
                bio::region *ths = rbs_region;
                rbs_region = rbs_region->next;
                delete ths;
            }
            delete root;
            if (rbs_region) {
                spacer_len = minlen - rbs_region->end + 1;
                rbs_len = rbs_region->end - rbs_region->start;
                if ((spacer_len >= 4 && spacer_len <= 13) &&
                    (rbs_len    >= 3 && rbs_len    <= 10)) break;
                delete rbs_region;
                rbs_region = nullptr;
            }
        }
        if (rbs_region == nullptr) {
            if (!QUIET) std::cerr << "Warning: detect no RBS consensus, using default\n";
            rbs_region = NEW bio::region(minlen-13, minlen-7);
        } else {
            char *cons = bio_util::get_consensus(r_flankings, rbs_region->start, rbs_region->end);
            if (!QUIET) std::cerr << "RBS Consensus: " << std::setw(35) << ("5'-" + std::string(cons) + "-3'\n")
                                  << "RBS Spacer Length: " << std::setw(27) << spacer_len << " bp\n"
                                  << "Number of Alter Starts: " << std::setw(25) << n_alters << "\n";
            delete[] cons;
        }
        
        bool flag = !(bool) args.count("longest");
        
        if (!QUIET) std::cerr << "Training RBS Model ...";
        if (flag) {
            /* encode flankings to zcurve params */
            params = NEW double[DIM_S*n_alters*2];
            if (params == nullptr) {
                std::cerr << MEM_ERR_INFO;
                return 1;
            }
            encoding::encode_seqs(p_flankings, minlen, params, 3);
            encoding::encode_seqs(q_flankings, minlen, params+n_alters*DIM_S, 3);
            encoding::minmax_scale(params, n_alters*2, DIM_S, mins, maxs);

            /* initiation of RBS model training */
            double*  rbs_data = model::build_rbs_data(p_flankings, cds_model, params, 
                                                      rbs_region, plot_range, exp_lens);
            double*  rbs_model = NEW double[8]();
            double** p_data = NEW double*[n_alters];
            double*  rbs_scores = NEW double[n_alters];
            if (rbs_data == nullptr || rbs_model == nullptr || p_data == nullptr || rbs_scores == nullptr) {
                std::cerr << MEM_ERR_INFO;
                return 1;
            }
            int ns_alt = n_alters;
            for (int i = 0; i < ns_alt; i ++) p_data[i] = rbs_data + 7*i;
            
            /* train fisher model */
            flag = model::fisher_train(p_data, ns_alt, 7, s_labels.data(), rbs_model);
            if (flag) {
                /* assign fisher scores to all alter starts */
                model::fisher_predict(rbs_data, n_alters, 7, rbs_model, rbs_scores);
                ns_alt = 0;
                for (int i = 0; i < n_init_seeds; i ++) {
                    bio::orf& pi = putative[init_seeds[i]];
                    /* find the best alter start */
                    int max_idx = -1;
                    double max_score = std::numeric_limits<double>::lowest();
                    for (int j = 0; j < pi.starts.size(); j ++) {
                        int score_idx = pi.scores[j];
                        if (rbs_scores[score_idx] > max_score) {
                            max_score = rbs_scores[score_idx];
                            max_idx = j;
                        }
                    }
                    int b_start = pi.starts[max_idx];
                    if (max_idx> -1 && max_score > 0) {
                        pi.seq = pi.seq + (b_start - pi.t_start);
                        pi.len = pi.end - b_start;
                        pi.t_start = b_start;
                    }
                }
            }

            if (!QUIET) std::cerr << std::setw(28) << (flag ? "Done\n" : "Skipped\n");
        } else if (!QUIET) std::cerr << std::setw(28) << "Bybassed\n";
    }

    if (cds_model) svm_free_model_content(cds_model);
    int num_putative = (int) putative.size();
    if (!QUIET) std::cerr << "Number of Putative Genes:" << std::setw(24) << num_putative << "\n";

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
    double min_ratio = args["ratio"].as<double>();
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