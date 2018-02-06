/*
 * Program: "cb-np_count-isoforms"
 * Requires: C++14, OpenMP, TCLAP
 * Purpose: Quickly count different transcript splicing isoforms present in the
 *          given alignment.
 * Author: Matthew Bauer
 */

#include <tclap/CmdLine.h>
#include <regex>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <stdint.h>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <omp.h>

double cutoff;
std::vector<std::pair<std::regex, double>> cutoff_groups;
std::ifstream alig_fin, anno_fin;

template<class A, class B>
B atob(A const& a) {
    std::stringstream ss;
    ss << a;
    B b;
    ss >> b;
    return b;
}

std::unordered_map<std::string, double> cutoff_table;
void register_feature_cutoff(std::string const& name) {
    double c = cutoff;
    for(auto const& group : cutoff_groups) {
        if(std::regex_match(name, group.first))
        {
            c = group.second;
            break;
        }
    }
    cutoff_table[name] = c;
}

struct Interval {
    uint64_t start;
    uint64_t end;
    Interval() {}
    Interval(uint64_t s, uint64_t e) : start(s), end(e) {
        if(start >= end) {
            std::stringstream ss;
            ss << "Invalid interval: ("
               << atob<uint64_t, std::string>(start) << ", "
               << atob<uint64_t, std::string>(end) << ")";
            throw std::invalid_argument(ss.str());
        }
    }
};
struct NamedInterval : public Interval {
    std::string name;
    NamedInterval() {}
    NamedInterval(uint64_t s, uint64_t e, std::string const& n)
        : Interval(s, e), name(n) {}
    NamedInterval(Interval const& ivl)
        : NamedInterval(ivl.start, ivl.end, "") {}
    NamedInterval(NamedInterval const& ivl)
        : NamedInterval(ivl.start, ivl.end, ivl.name) {}
};
bool operator<(Interval const& a, Interval const& b) {
    return a.end <= b.start;
}
int cmp_ivls(Interval const& a, Interval const& b) {
    if(a.end <= b.start)
        return -1;
    if(a.start >= b.end)
        return 1;
    return 0;
}
uint64_t get_overlap(Interval const& a, Interval const& b) {
    if(cmp_ivls(a, b) != 0)
        return 0;
    return std::min(a.end, b.end) - std::max(a.start, b.start);
}

std::vector<std::string> strsplit(
        std::string const& s, char d='\t', bool prune=false) {
    std::vector<std::string> ret;
    std::string growing;
    for(char ch : s) {
        if(ch != d)
            growing += ch;
        else
        {
            ret.push_back(growing);
            growing = "";
        }
    }
    if(!prune || !growing.empty())
        ret.push_back(growing);
    return ret;
}

void parse_args(int argc, char **argv)
{
    try {
        TCLAP::CmdLine cmd(
                ("Quickly count different transcript splicing isoforms present "
                 "in the alignment."),
                ' ', "1.0.0");
        TCLAP::UnlabeledValueArg<std::string> anno_arg(
                "annotation", "The BED6 annotation.", true, "", "file path");
        TCLAP::UnlabeledValueArg<std::string> alig_arg(
                "alignment", "The BED12 alignment.", true, "", "file path");
        TCLAP::ValueArg<double> cutoff_arg(
                "c", "cutoff", "The default cutoff for feature inclusion.",
                false, 50.0, "percentage");
        TCLAP::MultiArg<std::string> cutoffgroups_arg(
                "C", "cutoff-group",
                ("A group of features with a specific inclusion cutoff. "
                 "Format: 'R C', where R is a regex pattern that matches "
                 "the group in question and C is the cutoff that should be "
                 "used for the group. R and C are separated by a space."),
                false, "spec");
        cmd.add(anno_arg);
        cmd.add(alig_arg);
        cmd.add(cutoff_arg);
        cmd.add(cutoffgroups_arg);
        cmd.parse(argc, argv);

        // try to open files
        anno_fin.open(anno_arg.getValue());
        alig_fin.open(alig_arg.getValue());

        // check cutoff validity
        cutoff = cutoff_arg.getValue();
        if(cutoff < 0.0 || cutoff > 100.0)
            throw std::runtime_error(
                    std::string("Invalid cutoff: ")
                    + atob<double, std::string>(cutoff));

        for(std::string const& spec : cutoffgroups_arg.getValue()) {
            auto delim_pos = spec.find_first_of(' ');
            // check spec validity
            if(delim_pos >= spec.size() || delim_pos != spec.find_last_of(' '))
                throw std::runtime_error(
                        std::string("Invalid cutoff group specifier: ") + spec);
            std::string regex_part = std::string(
                    spec.begin(), spec.begin()+delim_pos);
            std::string cutoff_part = std::string(
                    spec.begin()+delim_pos+1, spec.end());
            try {
                // try to create <regex, float> pair
                std::regex r(regex_part);
                float c = std::stod(cutoff_part);
                if(c < 0.0 || c > 100.0)
                    throw std::invalid_argument("");
                cutoff_groups.emplace_back(r, c);
            }
            catch(std::regex_error const& err) {
                throw std::runtime_error(
                        std::string("Invalid regex: ") + regex_part);
            }
            catch(std::invalid_argument const& err) {
                throw std::runtime_error(
                        std::string("Invalid cutoff: ") + cutoff_part);
            }
        }
    }
    catch(TCLAP::ArgException& err) {
        std::stringstream ss;
        ss
            << "E: " << err.error()
            << " for arg " << err.argId();
        throw std::runtime_error(ss.str());
    }
}

class CoveredRegion
{
    Interval mivl;
    std::multiset<Interval> ivls;

public:
    CoveredRegion() {}
    CoveredRegion(Interval const& main_ivl) : mivl(main_ivl) {}
    void add(Interval ivl) {
        if(cmp_ivls(ivl, mivl) != 0)
            return; // new interval outside of main interval
        // cut down to size
        ivl = Interval(
                std::max(ivl.start, mivl.start),
                std::min(ivl.end, mivl.end));
        // merge into set
        uint64_t min_start = ivl.start;
        uint64_t max_end = ivl.end;
        auto overlapping = ivls.equal_range(ivl);
        for(auto it = overlapping.first; it != overlapping.second; it++) {
            min_start = std::min(min_start, it->start);
            max_end = std::max(max_end, it->end);
        }
        ivls.erase(overlapping.first, overlapping.second);
        ivls.emplace(min_start, max_end);
    }
    uint64_t get_raw_coverage() const {
        uint64_t ret = 0;
        for(Interval const& ivl : ivls)
            ret += ivl.end-ivl.start;
        return ret;
    }
    double get_percent_coverage() const {
        return 100.0*get_raw_coverage()/((double)(mivl.end-mivl.start));
    }
};

std::multiset<NamedInterval> annotation;
typedef std::unordered_map<std::string, CoveredRegion> CoverageCounter;
CoverageCounter empty_coverage_counter;
void load_annotation() {
    // load from file
    std::string line;
    int linenum = 0;
    while(linenum++, getline(anno_fin, line)) {
        if(line.empty())
            continue;
        std::vector<std::string> fields = strsplit(line);
        if(fields.size() < 4)
        {
            std::stringstream ss;
            ss << "Invalid annotation line: " << linenum;
            throw std::runtime_error(ss.str());
        }
        try
        {
            uint64_t start = atob<std::string, uint64_t>(fields[1]);
            uint64_t end = atob<std::string, uint64_t>(fields[2]);
            std::string name = fields[3];
            annotation.emplace(start, end, name);
            if(empty_coverage_counter.count(name))
            {
                std::stringstream ss;
                ss << "Feature name used more than once: " << name;
                throw std::runtime_error(ss.str());
            }
            empty_coverage_counter[name] = CoveredRegion(Interval(start, end));
            register_feature_cutoff(name);
        }
        catch(std::ios::failure const& err) {
            std::stringstream ss;
            ss
                << "Invalid annotation line: " << linenum << "; "
                << err.what();
            throw std::runtime_error(ss.str());
        }
    }

    // validate
    if(annotation.empty())
        throw std::runtime_error("Annotation is empty!");
    NamedInterval prev = *(annotation.begin());
    for(auto it = std::next(annotation.begin()); it != annotation.end(); it++) {
        if(cmp_ivls(prev, *it) != -1)
            throw std::runtime_error("Annotation records overlap!");
        prev = *it;
    }
}

std::vector<Interval> process_bed12_record(std::string const& rec) {
    std::runtime_error err("Invalid BED12 format!");
    std::vector<Interval> ret;
    std::vector<std::string> fields = strsplit(rec);
    if(fields.size() != 12)
        throw err;
    int64_t start, end;
    try {
        start = atob<std::string, int64_t>(fields[1]);
        end = atob<std::string, int64_t>(fields[2]);
    }
    catch(std::ios::failure const& e) {
        throw err;
    }
    std::vector<std::string> blk_lens_strs = strsplit(
            fields[10], ',', true);
    std::vector<std::string> blk_starts_strs = strsplit(
            fields[11], ',', true);
    if(blk_lens_strs.size() != blk_starts_strs.size())
        throw err;
    try {
        for(int i = 0; i < blk_lens_strs.size(); i++) {
            int64_t blk_start = atob<std::string, int64_t>(blk_starts_strs[i]);
            int64_t blk_len = atob<std::string, int64_t>(blk_lens_strs[i]);
            ret.emplace_back(start+blk_start, start+blk_start+blk_len);
        }
    }
    catch(std::ios::failure const& e) {
        throw err;
    }
    return ret;
}

std::vector<std::vector<Interval>> reads;
void load_all_reads() {
    std::string line;
    int linenum = 0;
    while(linenum++, std::getline(alig_fin, line)) {
        try {
            reads.push_back(process_bed12_record(line));
        }
        catch(std::runtime_error const& err) {
            std::stringstream ss;
            ss
                << "Error on alignment line " << linenum << ": "
                << err.what();
            throw std::runtime_error(ss.str());
        }
    }
}

std::string serialize_coverage_counter(CoverageCounter const& cc) {
    std::vector<std::string> included;
    for(auto const& feat : cc)
        if(feat.second.get_percent_coverage() >= cutoff_table[feat.first])
            included.push_back(feat.first);

    std::sort(included.begin(), included.end());
    std::string ret;
    if(!included.empty()) {
        ret += included[0];
        for(auto it = included.begin()+1; it != included.end(); it++)
        {
            ret += ",";
            ret += *it;
        }
    }
    return ret;
}

int main(int argc, char **argv) {
    try {
        parse_args(argc, argv);

        load_annotation();
        anno_fin.close();

        load_all_reads();
        alig_fin.close();
    }
    catch(std::runtime_error const& err) {
        alig_fin.close();
        anno_fin.close();
        std::cerr << err.what() << std::endl;
        return -1;
    }

    int N = reads.size();
    std::unordered_map<std::string, uint32_t> final_count;
#pragma omp parallel
    {
        int i;
        std::string this_iso;
        CoverageCounter cc;
        std::unordered_map<std::string, uint32_t> intermediate_count;

        // (in parallel) process each read
#pragma omp for schedule(static) nowait
        for(i = 0; i < N; i++) {
            cc = empty_coverage_counter;
            for(Interval const& ivl : reads[i])
            {
                auto eqr = annotation.equal_range(ivl);
                for(auto it = eqr.first; it != eqr.second; it++)
                    cc[it->name].add(ivl);
            }
            this_iso = serialize_coverage_counter(cc);
            if(!this_iso.empty()) {
                if(!intermediate_count.count(this_iso))
                    intermediate_count[this_iso] = 0;
                intermediate_count[this_iso]++;
            }
        }

        // merge this thread's results with the final results
#pragma omp critical
        {
            for(auto const& p : intermediate_count) {
                if(!final_count.count(p.first))
                    final_count[p.first] = 0;
                final_count[p.first] += p.second;
            }
        }
    }

    // sort & print results
    std::vector<std::pair<std::string, uint32_t>> final_count_sorted(
            final_count.begin(), final_count.end());
    std::sort(
            final_count_sorted.begin(), final_count_sorted.end(),
            [](auto const& a, auto const& b){return a.second < b.second;});
    for(auto const& p : final_count_sorted)
        std::cout << p.first << "\t" << p.second << std::endl;

    return 0;
}
