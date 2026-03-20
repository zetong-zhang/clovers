#include "BioIO.hpp"

/**
 * @brief   Check if a file is gzipped.
 * 
 * @param filename  The path to the file.
 * @return true If the file is gzipped.
 *         false If the file is not gzipped.
 */
static bool is_gzip(const std::string &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) return false;
    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    // 1st and 2nd magic of gzip are 0x1F and 0x8B
    return (magic[0] == 0x1F && magic[1] == 0x8B);
}

bool bio_io::file_exists(std::string &path) noexcept {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

/**
 * @brief   Read a FASTA or GenBank file from a stream.
 * 
 * @param handle    The input stream.
 * @param scaffolds The record array to store the scaffolds.
 * @return true If the file is successfully read.
 *         false If the file is not successfully read.
 */
static bool read_stream(
    std::istream &handle,
    bio::record_array &scaffolds
) {
    std::string header;
    // skip empty lines
    while(std::getline(handle, header))
        if (!header.empty()) break;
    if (header.empty()) {
        std::cerr << "\nError: empty content \n";
        return false;
    }
    std::string line, buffer;
    // FASTA files should start with '>'
    if (header.at(0) == '>') {
        while (std::getline(handle, line)) {
            if (line.empty()) continue;
            else if (line.at(0) == '>') {
                size_t sep = header.find(' ');
                std::string entry = sep>1 ? header.substr(1, sep-1):header.substr(1);
                if (!entry.empty()) {
                    if (entry.back() == '\r') entry.pop_back();
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    header = std::move(line);
                    buffer.clear();
                } else {
                    std::cerr << "\nError: invalid FASTA header \n";
                    return false;
                }
            } else {
                if (!std::isalpha(line.back())) line.pop_back();
                buffer.append(std::move(line));
            }
            line.clear();
        }
        if (!buffer.empty()) {
            size_t sep = header.find(' ');
            std::string entry = sep>1 ? header.substr(1, sep-1):header.substr(1);
            scaffolds.emplace_back(std::move(entry), std::move(buffer));
        }
        return true;
    // GenBank files should start with 'LOCUS ... '
    } else if (!header.find("LOCUS")) {
        bool in_origin = false;
        std::string entry;
        std::istringstream header_stream(header);
        header_stream >> entry;
        if (!(header_stream >> entry)) {
            std::cerr << "\nError: invalid LOCUS line \n";
            return false;
        }
        while (std::getline(handle, line)) {
            if (line.empty()) continue;
            else if (!line.find("ORIGIN")) in_origin = true;
            else if (!line.find("//")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_origin = false;
            } else if (!line.find("LOCUS")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_origin = false; 

                std::istringstream header_stream(line);
                header_stream >> entry;
                if (!(header_stream >> entry)) {
                    std::cerr << "\nError: invalid LOCUS line \n";
                    return false;
                }
            } else if (in_origin) {
                for (char c : line) if (std::isalpha(c)) buffer += std::toupper(c);
            }
        } 
        if (!buffer.empty()) {
            scaffolds.emplace_back(std::move(entry), std::move(buffer));
        }
        return true;
    // EMBL files should start with 'ID   '
    } else if (!header.find("ID   ")) {
        bool in_sequence = false;
        std::string entry;
        std::istringstream header_stream(header);
        header_stream >> entry;
        if (!(header_stream >> entry)) {
            std::cerr << "\nError: invalid ID line \n";
            return false;
        }
        while (std::getline(handle, line)) {
            if (line.empty()) continue;
            else if (!line.find("SQ   ")) in_sequence = true;
            else if (!line.find("//")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_sequence = false;
            } else if (!line.find("ID   ")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_sequence = false;

                std::istringstream header_stream(line);
                header_stream >> entry;
                if (!(header_stream >> entry)) {
                    std::cerr << "\nError: invalid ID line \n";
                    return false;
                }
            } else if (in_sequence) {
                for (char c : line) if (std::isalpha(c)) buffer += std::toupper(c);
            }
        }
        if (!buffer.empty()) {
            scaffolds.emplace_back(std::move(entry), std::move(buffer));
        }
        return true;
    } else {
        std::cerr << "\nError: unsupported file type (options: GenBank, EMBL, FASTA)\n";
        return false;
    }
}
/**
 * @brief   Write the ORFs to a result file.
 * 
 * @param handle    The output stream.
 * @param orfs      The ORF array to be written.
 * @param is_circ   Whether the scaffold is circular.
 * @param format    The format of the output file.
 * @return true If the file is successfully written.
 *         false If the file is not successfully written.
 */
static bool write_stream(
    std::ostream &handle,
    bio::orf_array &orfs,
    const std::string &date,
    bool is_circ,
    const std::string &format
) {
    int count = (int) orfs.size();
    if (format == "gff") {
        handle << "##gff-version 3 \n";
        char *last_scaffold = nullptr;
        for (int i = 0, j = 0; i < count; i ++) {
            // seqid + source + type
            if (last_scaffold != orfs[i].host) {
                handle << "# " << orfs[i].host << '\t' << orfs[i].host_len << " bp\t"
                       << (is_circ ? "circular" : "linear") << "\tUNA\t" << date << '\n';
                last_scaffold = orfs[i].host;
            }
            handle << orfs[i].host << '\t' << VERSION << "\tCDS\t";
            int rstart, rend, host_len = orfs[i].host_len;
            if (orfs[i].strand == '+') {
                rstart = orfs[i].t_start + 1;
                rend = orfs[i].end;
            } else {
                rstart = host_len - orfs[i].end + 1;
                rend = host_len - orfs[i].t_start;
            }
            // start + end
            handle << rstart << '\t' << rend << '\t';
            // score
            handle << std::fixed << std::setprecision(3) 
                   << orfs[i].score << '\t';
            // strand + phase
            handle << orfs[i].strand << "\t0\tID=orf";
            // attributes
            handle << std::setw(6) << std::setfill('0') << (++j);
            if (orfs[i].partial3) handle << ";partial3=true\n";
            else if (orfs[i].partial5) handle << ";partial5=true\n";
            else handle << "\n";
        }
        return true;
    } else if (format == "gbk") {
        char *last_scaffold = nullptr;
        for (int i = 0, j = 0; i < count; i ++) {
            if (last_scaffold != orfs[i].host) {
                if (last_scaffold) handle << "ORIGIN\n//\n";
                handle << "LOCUS       " << orfs[i].host << ' ' << std::setw(27-(int)std::strlen(orfs[i].host)) 
                       << orfs[i].host_len << " bp" << "    DNA" << std::setw(13) 
                       << (is_circ ? "circular" : " linear") << " UNA " << date << '\n';
                handle << "DEFINITION  " << orfs[i].host << "\n";
                last_scaffold = orfs[i].host;
                handle << "FEATURES             Location/Qualifiers\n";
            }
            bool neg_strand = (bool)(orfs[i].strand == '-');
            handle << "     CDS             ";
            int rstart, rend, host_len = orfs[i].host_len;
            if (!neg_strand) {
                rstart = orfs[i].t_start + 1;
                rend = orfs[i].end;
            } else {
                rstart = host_len - orfs[i].end + 1;
                rend = host_len - orfs[i].t_start;
                handle << "complement(";
            }
            if (rstart > rend) {
                handle << "join(" << rstart << ".." 
                << orfs[i].host_len << ",1.." << rend << ')';
            } else {
                if (orfs[i].partial5 && !neg_strand || 
                    orfs[i].partial3 && neg_strand) handle << '<';
                handle << rstart << "..";
                if (orfs[i].partial3 && !neg_strand || 
                    orfs[i].partial5 && neg_strand) handle << '>';
                handle << rend;
            }
            if (neg_strand) handle << ")";
            handle << "\n                     /note=\"version=" << VERSION 
                   << ";ID=orf" << std::setw(6) << std::setfill('0') << (++j)
                   << ";score=" << std::fixed << std::setprecision(3) 
                   << orfs[i].score << "\"\n";
        }
        if (count > 0) handle << "ORIGIN\n//\n";
        return true;
    } else if (format == "med") {
        char *last_scaffold = nullptr;
        for (int i = 0, j = 0; i < count; i ++) {
            if (last_scaffold != orfs[i].host) {
                handle << "# " << orfs[i].host << '\t' << orfs[i].host_len << " bp\t"
                       << (is_circ ? "circular" : "linear") << "\tUNA\t" << date << '\n';
                last_scaffold = orfs[i].host;
            }
            bool neg_strand = (bool)(orfs[i].strand == '-');
            int rstart, rend, host_len = orfs[i].host_len;
            if (!neg_strand) {
                rstart = orfs[i].t_start + 1;
                rend = orfs[i].end;
            } else {
                rstart = host_len - orfs[i].end + 1;
                rend = host_len - orfs[i].t_start;
            }
            handle << rstart << ' ' << rend << '\t' << orfs[i].strand << "\n";
        }
        return true;
    }
    std::cerr << "Error: unsupported output format '" 
              << format << "' (options: gff, gbk, med)\n";
    return false;
}
// API functions 

bool bio_io::read_source(
    const std::string &filename,
    bio::record_array &scaffolds
) {
    if (filename == "-") return read_stream(std::cin, scaffolds);
#ifdef ZLIB
    else if (is_gzip(filename)) {
        gzFile file = gzopen(filename.c_str(), "rb");
        assert(file != nullptr);
        char buff[BUFF_SIZE];
        int bytes_read, err_code;
        std::string content;
        while ((bytes_read = gzread(file, buff, BUFF_SIZE)))
            content.append(buff, bytes_read);
        const char *err_msg = gzerror(file, &err_code);
        if (err_code != Z_OK && err_code != Z_STREAM_END) {
            gzclose(file);
            std::cerr << "\nError: failed to decompress " << filename << "\n";
            return false;
        }
        gzclose(file);
        std::istringstream handle(std::move(content));
        return read_stream(handle, scaffolds);
    }  
#endif
    else {
        std::ifstream handle(filename);
        if (!handle.is_open()) {
            if (filename != "example.fa") 
                std::cerr << "\nError: cannot open " << filename << "\n";
            return false;
        }
        return read_stream(handle, scaffolds);
    }
}

bool bio_io::write_result(
    bio::orf_array &orfs, 
    const std::string &date,
    bool is_circ,
    const std::string &filename,
    const std::string &format
) {
    if (filename == "-") return write_stream(std::cout, orfs, date, is_circ, format);
    else {
        std::ofstream handle(filename);
        if (!handle.is_open()) {
            std::cerr << "\nError: failed to open " << filename << '\n';
            return false;
        }
        return write_stream(handle, orfs, date, is_circ, format);
    }
}

bool bio_io::write_faa(
    bio::orf_array &orfs, 
    const std::string &filename
) {
    std::ofstream handle(filename);
    if (!handle.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    int count = (int) orfs.size();
    for (int i = 0; i < count; i ++) {
        int rstart, rend, host_len = orfs[i].host_len;
        if (orfs[i].strand == '+') {
            rstart = orfs[i].t_start + 1;
            rend = orfs[i].end;
        } else {
            rstart = host_len - orfs[i].end + 1;
            rend = host_len - orfs[i].t_start;
        }
        handle << '>' << orfs[i].host << ':' << rstart << ".." << rend << '(' 
               << orfs[i].strand << ") score=" << std::fixed << std::setprecision(3) 
               << orfs[i].score << '\n';
        int prolen = (int) (orfs[i].len / 3);
        char *protein = bio_util::gene2protein(orfs[i], prolen);
        if (!protein) {
            std::cerr << "Warning: invalid characters encountered in "
                      << orfs[i].host << ':' << rstart 
                      << ".." << rend << '(' << orfs[i].strand << ")\n";
            continue;
        }
        for (int j = 0; j < prolen; j += 80) {
            int line_length = std::min(80, prolen - j);
            handle.write(protein + j, line_length);
            handle << '\n';
        }
        delete[] protein;
    }
    handle.close();
    return true;
}

bool bio_io::write_fna(
    bio::orf_array &orfs, 
    const std::string &filename
) {
    std::ofstream handle(filename);
    if (!handle.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    int count = (int) orfs.size();
    for (int i = 0; i < count; i ++) {
        int rstart, rend, host_len = orfs[i].host_len;;
        if (orfs[i].strand == '+') {
            rstart = orfs[i].t_start + 1;
            rend = orfs[i].end;
        } else {
            rstart = host_len - orfs[i].end + 1;
            rend = host_len - orfs[i].t_start;
        }
        handle << '>' << orfs[i].host << ':' << rstart << ".." << rend << '(' 
               << orfs[i].strand << ") score=" << std::fixed << std::setprecision(3) 
               << orfs[i].score << '\n';
        for (int j = 0; j < orfs[i].len; j += 80) {
            int line_length = std::min(80, orfs[i].len - j);
            handle.write(orfs[i].seq + j, line_length);
            handle << '\n';
        }
    }
    handle.close();
    return true;
}

bool bio_io::write_model(
    svm_model *model,
    float *mins, float *maxs,
    const std::string &filename
) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    try {
        // write scaler
        outfile.write((char *) mins, DIM_S * sizeof(float));
        outfile.write((char *) maxs, DIM_S * sizeof(float));
        // write gamma
        float gamma = model->param.gamma;
        outfile.write((char*) &gamma, sizeof(float));
        // write number of support vectors
        int n_sv = model->l;
        outfile.write((char*) &n_sv, sizeof(int));
        // write support vectors
        for (int i = 0; i < n_sv; i ++) {
            float *sv = model->SV[i];
            outfile.write((char*) sv, DIM_S * sizeof(float));
        }
        // write coefficients
        for (int i = 0; i < n_sv; i ++) {
            float coef = model->sv_coef[0][i];
            outfile.write((char*) &coef, sizeof(float));
        }
        // write intercept
        float intercept = model->rho[0];
        outfile.write((char*) &intercept, sizeof(float));
    } catch (const std::ios_base::failure& e) {
        std::cerr << "\nError: failed to write svm model to " << filename << '\n';
        return false;
    }
    return true;
}

bool bio_io::write_model(
    float *params, int max_alter,
    float pFU, float pFD,
    const std::string &filename
) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    try {
        outfile.write((char*) params, TIS_S * sizeof(float));
        outfile.write((char*) &max_alter, sizeof(int));
        outfile.write((char*) &pFU, sizeof(float));
        outfile.write((char*) &pFD, sizeof(float));
    } catch (const std::ios_base::failure& e) {
        std::cerr << "\nError: failed to write TriTISA+ model to " << filename << '\n';
        return false;
    }
    return true;
}

bool bio_io::read_model(
    svm_model *model,
    float *mins, float *maxs,
    const std::string &filename
) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    try {
        // read scaler
        infile.read((char *) mins, DIM_S * sizeof(float));
        infile.read((char *) maxs, DIM_S * sizeof(float));
        // read gamma
        float gamma;
        infile.read((char*) &gamma, sizeof(float));
        model->param.gamma = gamma;
        // read number of support vectors
        int n_sv;
        infile.read((char*) &n_sv, sizeof(int));
        model->l = n_sv;
        // read support vectors
        float *sv_cache = NEW float[n_sv * DIM_S];
        if (!sv_cache) return false;
        infile.read((char*) sv_cache, n_sv * DIM_S * sizeof(float));
        model->SV = NEW float*[n_sv];
        if (!model->SV) return false;
        for (int i = 0; i < n_sv; i ++) {
            model->SV[i] = sv_cache + i * DIM_S;
        }
        // read coefficients
        model->sv_coef = NEW float*[1];
        model->sv_coef[0] = NEW float[n_sv];
        if (!model->sv_coef[0]) return false;
        infile.read((char*) model->sv_coef[0], n_sv * sizeof(float));
        // read intercept
        model->rho = NEW float[1];
        infile.read((char*) &model->rho[0], sizeof(float));
    } catch (const std::ios_base::failure& e) {
        std::cerr << "\nError: failed to read svm model from " << filename << '\n';
        return false;
    }
    return model;
}

bool bio_io::read_model(
    float *params, int &max_alter,
    float &pFU, float &pFD,
    const std::string &filename
) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    try {
        infile.read((char*) params, TIS_S * sizeof(float));
        infile.read((char*) &max_alter, sizeof(int));
        infile.read((char*) &pFU, sizeof(float));
        infile.read((char*) &pFD, sizeof(float));
    } catch (const std::ios_base::failure& e) {
        std::cerr << "\nError: failed to read TriTISA+ model from " << filename << '\n';
        return false;
    }
    return true;
}

bool bio_io::write_overprint(
    bio::orf_array &orfs,
    const float ratio,
    const int min_olen,
    const std::string &filename
) {
    std::ofstream handle(filename);
    if (!handle.is_open()) {
        std::cerr << "\nError: failed to open " << filename << '\n';
        return false;
    }
    handle << "##gff-version 3\n";
    int count = (int) orfs.size() / 2;
    for (int i = 0; i < count; i ++) {
        bio::orf &gene_1 = orfs[2*i], &gene_2 = orfs[2*i+1];
        op_type type = bio_util::check_overprint(gene_1, gene_2, ratio, min_olen);
        handle << "# index=" << std::setw(4) << std::setfill('0') << (i+1) << ";type=" 
               << (type == op_type::INTERSECT ? "INTERSECT\n" : "INCLUDE\n");
        int rstart, rend, host_len;

        host_len = gene_1.host_len;
        if (gene_1.strand == '+') {
            rstart = gene_1.t_start + 1;
            rend = gene_1.end;
        } else {
            rstart = host_len - gene_1.end + 1;
            rend = host_len - gene_1.t_start;
        }
        handle << gene_1.host << '\t' << VERSION << "\tCDS\t";
        // start + end
        handle << rstart << '\t' << rend << '\t';
        // score
        handle << std::fixed << std::setprecision(3) 
               << gene_1.score << '\t';
        // strand + phase
        handle << gene_1.strand << "\t0\t.\n";

        host_len = gene_2.host_len;
        if (gene_2.strand == '+') {
            rstart = gene_2.t_start + 1;
            rend = gene_2.end;
        } else {
            rstart = host_len - gene_2.end + 1;
            rend = host_len - gene_2.t_start;
        }
        // seqid + source + type
        handle << gene_2.host << '\t' << VERSION << "\tCDS\t";
        // start + end
        handle << rstart << '\t' << rend << '\t';
        // score
        handle << std::fixed << std::setprecision(3)
               << gene_2.score << '\t';
        // strand + phase
        handle << gene_2.strand << "\t0\t.\n";
    }
    handle.close();
    return true;
}