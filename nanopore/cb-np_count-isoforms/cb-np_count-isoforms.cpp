/*
 * Program: "cb-np_count-isoforms"
 * Requires: C++17, TCLAP
 * Purpose: Quickly count different transcript splicing isoforms present in the
 *          given alignment.
 * Author: Matthew Bauer
 */

#include <tclap/CmdLine.h>
#include <string>
#include <vector>
#include <stdint.h>
#include <functional>
#include <algorithm>
#include <regex>
#include <unordered_map>
#include <set>
#include <fstream>
#include <iostream>

/* PROTOTYPES */

/** Split a string by the given delimiter.
 *
 * Single delimiters at the very end of the string will be ignored.
 */
template<char delim='\t'>
std::vector<std::string> split(std::string const& str);

/** Convert from type A to type B using a stringstream. */
template<class A, class B>
B ss_conv(A const& a);

/* Shortcut for creating error messages. */
/** Useful, flexible exception class. */
class MyErr {
protected:
    std::string msg;

public:
    /** Compose the error.
     *
     * Takes any number of parameters of any time, as long as they are
     * stringstream-compatible.
     */
    template<class... Args>
    explicit MyErr(Args const&... args);

    friend std::ostream& operator<<(std::ostream& out, MyErr const& err);

    /** Combine two errors in a stack-trace-like manner. */
    friend MyErr operator+(MyErr const& a, MyErr const& b);
};

/** Object to store and parse command-line options. */
struct Config {
    std::string path_anno, /** The annotation file. */
        path_alig, /** The alignment file. */
        path_distdir; /** The directory to save coverage distributions in. */

    std::string default_cutoff; /** Specifies the default cutoff. */
    std::vector<std::string> cutoff_groups; /** REGEX-Cutoff pairs. */

    /** Parse the given command-line args. */
    explicit Config(int argc, char **argv);
};

/** A half-open interval: [start, end). */
struct Interval {
    uint64_t start;
    uint64_t end;

    /** [0, 1). */
    Interval();

    /** [s, e). */
    Interval(uint64_t s, uint64_t e);

    /** Test if a is strictly (no overlap) less than b.
     *
     * This is convenient in multisets because it makes it very eassy to find
     * overlapping intervals. */
    friend bool operator<(Interval const& a, Interval const& b);
};

/* Represents an annotation feature. */
struct Feature : public Interval {
    std::string name;

    /** Nameless, default interval. */
    Feature();
    Feature(uint64_t s, uint64_t e, std::string const& n);

    /** Nameless copy. */
    Feature(Interval const& ivl);

    /** Parse a BED6 record. */
    Feature(std::string const& bed6_record);
};

/** Represents a single read/alignment.
 *
 * Use as an iterable of intervals. */
class Alignment : public std::vector<Interval> {
public:
    /** Parse a BED12 record. */
    explicit Alignment(std::string const& bed12_record);
};

/** Represents the entire annotation.
 *
 * Use as an iterable of Features. */
class Annotation : public std::multiset<Feature> {
public:
    /** Load from a BED6 file. */
    explicit Annotation(std::string const& fpath);

    /** Get all Features that overlap the given Interval. */
    std::vector<std::reference_wrapper<Feature const>> get_overlapping(
            Interval const& ivl) const;
};

/** Convenience class representing a percent value. */
struct Percent {
    double val;

    Percent(); /** 0% */
    Percent(double d); /** Throw if invalid percent. */

    operator double() const;
    friend std::istream& operator>>(std::istream& inp, Percent& p);
};

/** A two-tailed cutoff (stores both inclusion and exclusion cutoffs). */
struct Cutoff {
    Percent inclusion; /** Inclusion minimum (inclusive). */
    Percent exclusion; /** Exclusion maximum (inclusive). */

    Cutoff();
    Cutoff(Percent i, Percent e);
    explicit Cutoff(std::string const& s); /** Parse: "inclusion,exclusion". */

    /** Advice given to the caller of advise(). */
    enum class Action {
        INCLUDE, /** Inclusion cutoff met; the feature should be counted. */
        EXCLUDE, /** Exclusion cutoff met; the feature should be not counted. */
        REMOVE, /** Neither cutoff met; the read should be thrown out! */
    };

    /** Decide what to do with the given percent coverage. */
    Action advise(Percent coverage) const;
};

/** Abstract interface for generating cutoffs from features. */
class CutoffGenerator {
protected:
    CutoffGenerator();

public:
    virtual Cutoff get_cutoff(Feature const& f) const =0;
    Cutoff operator[](Feature const& f) const; /** Alias for get_cutoff(). */
};

/** REGEX-group cabable "raw" percentage cutoff generator. */
class PercentCG : public CutoffGenerator {
public:
    /** Represents a group of features with common cutoffs. */
    struct FeatureGroup {
        Cutoff cutoff;
        std::regex pattern;

        /** Parse: "pattern inclusion,exclusion". */
        explicit FeatureGroup(std::string const& spec);
    };

protected:
    /** The cutoff used for features that don't match any group. */
    Cutoff default_cutoff; 
    std::vector<FeatureGroup> groups;

public:
    PercentCG(Cutoff def, std::vector<std::string> const& specs);
    virtual Cutoff get_cutoff(Feature const& f) const;
};

/** A generic region that can be partially covered with any number of
 *  intervals. */
class CoveredInterval : public Interval {
protected:
    std::multiset<Interval> covered;

public:
    /** [0, 1). */
    CoveredInterval();
    /** [s, e). */
    CoveredInterval(uint64_t s, uint64_t e);

    /** Cover part of the interval. */
    void absorb(Interval const& ivl);

    /** The (integer) number of bases covered overall. */
    uint64_t get_raw_coverage() const;
    /** The percentage of bases covered. */
    double get_percent_coverage() const;
};

/** Generates/decides isoforms given the cutoff configuration, the annotation,
 * and the alignments. */
class IsoformGenerator {
public:
    /** An isoform.
     *
     * Includes info about both the included features and the raw coverage
     * information. */
    class Isoform : public std::unordered_map<std::string, CoveredInterval> {
    public:
        std::string repr;
        bool good;

        operator std::string() const;
        operator bool() const;
    };

protected:
    using CutoffTable = std::unordered_map<std::string, Cutoff>;
    CutoffTable cutoffs; /** Memo table. */
    Isoform blank_iso; /** Isoform tempalte. */
    Annotation const& anno; /** A reference to the annotation. */

public:
    /** Generate template isoform; fill in cutoff table. */
    IsoformGenerator(Annotation const& anno, CutoffGenerator const& cgen);

    /** Decide on an isoform for a raw alignment. */
    Isoform generate(Alignment const& alig);
};

/** Break an interval into bins and count hits inside each bin.
 *
 * Useful for generating histograms. */
class BinCounter : public std::vector<unsigned int> {
protected:
    double min, max;
    double binsz;

public:
    BinCounter();
    BinCounter(double min, double max, unsigned int n);

    /** Count the given value sensibly. */
    void add(double val);

    /** Combine two BinCounters.
     *
     * They must have the same range and number of bins! */
    friend BinCounter operator+(BinCounter const& a, BinCounter const& b);

    /** Write a (tsv) table. */
    friend std::ostream& operator<<(std::ostream& out, BinCounter const& bc);
};

/** Object for keeping track of the coverage distributions for all features. */
class CovDistCounter {
protected:
    std::unordered_map<std::string, BinCounter> dists;

public:
    CovDistCounter(Annotation const& anno);

    /** Incorporate the given Isoform. */
    void count(IsoformGenerator::Isoform const& iso);

    /** Save each coverage distribution in a separate file in the given
     * directory. */
    void save(std::string const& dir) const; 
};

/* IMPLEMENTATIONS */
template<char delim='\t'>
std::vector<std::string> split(std::string const& str) {
    std::vector<std::string> ret;
    std::string growing;
    for (char ch : str) {
        if (ch == delim) {
            ret.push_back(growing);
            growing.clear();
        }
        else {
            growing += ch;
        }
    }
    if (!growing.empty())
        ret.push_back(growing);
    return ret;
}

template<class A, class B>
B ss_conv(A const& a) {
    std::stringstream ss;
    ss << a;
    B b;
    try {
        ss >> b;
        if (!ss) {
            throw MyErr("failed to convert: ", a);
        }
        return b;
    }
    catch (std::ios_base::failure const& err) {
        throw MyErr("failed to convert: ", a);
    }
}

template<class... Args>
MyErr::MyErr(Args const&... args) {
    std::stringstream ss;
    (ss << ... << args);
    msg = ss.str();
}
std::ostream& operator<<(std::ostream& out, MyErr const& err) {
    out << err.msg;
    return out;
}
MyErr operator+(MyErr const& a, MyErr const& b) {
    return MyErr(a.msg + "\n\t" + b.msg);
}

Config::Config(int argc, char **argv) {
    using namespace TCLAP;
    try {
        CmdLine cmd(
                ("Quickly count different transcript splicing isoforms present "
                 "in the alignment."),
                ' ', "2.0.0");
        UnlabeledValueArg<std::string> anno_arg(
                "annotation", "The BED6 annotation.", true, "", "anno.bed");
        UnlabeledValueArg<std::string> alig_arg(
                "alignment", "The BED12 alignment.", true, "", "alig.bed");
        ValueArg<std::string> defcut_arg(
                "c", "cutoff",
                "The default cutoff. Format: inclusion,exclusion",
                true, "", "i,e");
        MultiArg<std::string> cutgroups_arg(
                "g", "cutoff-group",
                "A cutoff group. Format: 'regex inclusion,exclusion'",
                false, "pat i,e");
        ValueArg<std::string> covdistdir_arg(
                "d", "cov-dist-dir",
                "The directory to store coverage distribution tables in.",
                false, "", "dir");
        cmd.add(anno_arg);
        cmd.add(alig_arg);
        cmd.add(defcut_arg);
        cmd.add(cutgroups_arg);
        cmd.add(covdistdir_arg);
        cmd.parse(argc, argv);

        path_anno = anno_arg.getValue();
        path_alig = alig_arg.getValue();
        path_distdir = covdistdir_arg.getValue();

        default_cutoff = defcut_arg.getValue();
        cutoff_groups = cutgroups_arg.getValue();
    }
    catch (ArgException const& err) {
        throw MyErr(err.error(), " for arg ", err.argId());
    }
}

Interval::Interval() : Interval(0, 1) {}
Interval::Interval(uint64_t s, uint64_t e) {
    start = s;
    end = e;
    if (end <= start) {
        throw MyErr("invalid interval: (", s, ", ", e, ")");
    }
}
bool operator<(Interval const& a, Interval const& b) {
    return a.end <= b.start;
}

Feature::Feature() : Feature(0, 0, "") {}
Feature::Feature(uint64_t s, uint64_t e, std::string const& n)
    : Interval(s, e) {
    name = n;
}
Feature::Feature(Interval const& ivl)
    : Feature(ivl.start, ivl.end, "")
{}
Feature::Feature(std::string const& bed6_rec) {
    std::vector<std::string> fields = split(bed6_rec);
    if (fields.size() < 4) {
        throw MyErr("incomplete BED6 record: ", bed6_rec);
    }

    start = ss_conv<std::string, uint64_t>(fields[1]);
    end = ss_conv<std::string, uint64_t>(fields[2]);
    name = fields[3];

    if (end <= start) {
        throw MyErr("invalid interval: (", start, ", ", end, ")");
    }
}

Alignment::Alignment(std::string const& bed12_rec) {
    std::vector<std::string> fields = split(bed12_rec);
    if (fields.size() < 12) {
        throw MyErr("incomplete BED12 record: ", bed12_rec);
    }

    try {
        uint64_t start = ss_conv<std::string, uint64_t>(fields[1]);
        uint64_t end = ss_conv<std::string, uint64_t>(fields[2]);

        std::vector<std::string> block_starts = split<','>(fields[11]);
        std::vector<std::string> block_sizes = split<','>(fields[10]);
        if (block_starts.size() != block_sizes.size()) {
            throw MyErr("inconsistent number of blocks");
        }

        for (int i = 0; i < block_starts.size(); i++) {
            uint64_t b_start = start + ss_conv<std::string, uint64_t>(block_starts[i]);
            uint64_t b_end = b_start + ss_conv<std::string, uint64_t>(block_sizes[i]);
            push_back(Interval(b_start, b_end));
        }
    }
    catch (MyErr const& err) { // from ss_conv
        throw err + MyErr("in BED12 record: ", bed12_rec);
    }
}

Annotation::Annotation(std::string const& fpath) {
    std::ifstream fin(fpath);
    if (!fin) {
        throw MyErr("failed to open annotation: ", fpath);
    }

    std::string line;
    while (std::getline(fin, line), fin) {
        insert(Feature(line));
    }

    // check that no intervals overlap
    Feature const* last = nullptr;
    for (Feature const& ivl : *this) {
        if (last && !(*last < ivl)) {
            throw MyErr("annotation features overlap: ",
                    last->name, ", ", ivl.name);
        }
        last = &ivl;
    }
}
std::vector<std::reference_wrapper<Feature const>> Annotation::get_overlapping(
        Interval const& ivl) const {
    auto eqr = equal_range(ivl);
    std::vector<std::reference_wrapper<Feature const>> ret;
    for (auto it = eqr.first; it != eqr.second; it++) {
        ret.emplace_back(*it);
    }
    return ret;
}

Percent::Percent() : Percent(0.0) {}
Percent::Percent(double d) {
    if (d < 0.0 || d > 100.0) {
        throw MyErr("invalid percent: ", d);
    }
    val = d;
}
Percent::operator double() const {
    return val;
}
std::istream& operator>>(std::istream& inp, Percent& p) {
    inp >> p.val;
    if (p.val < 0.0 || p.val > 100.0) {
        throw MyErr("invalid percent: ", p.val);
    }
    return inp;
}

Cutoff::Cutoff() : Cutoff(0, 0) {}
Cutoff::Cutoff(Percent i, Percent e) {
    if (i < e) {
        throw MyErr("inclusion cutoff (", i,
                ") must be >= exclusion cutoff (", e, ")");
    }
    inclusion = i;
    exclusion = e;
}
Cutoff::Cutoff(std::string const& s) {
    std::vector<std::string> fields = split<','>(s);
    if (fields.size() != 2) {
        throw MyErr("improper number of fields in cutoff spec: ", s);
    }
    try {
        inclusion = ss_conv<std::string, Percent>(fields[0]);
        exclusion = ss_conv<std::string, Percent>(fields[1]);
        if (inclusion < exclusion) {
            throw MyErr("inclusion cutoff (", inclusion,
                    ") must be >= exclusion cutoff (", exclusion, ")");
        }
    }
    catch (MyErr const& err) {
        throw err + MyErr("in parsing of cutoff spec: ", s);
    }
}
Cutoff::Action Cutoff::advise(Percent val) const {
    if (val >= inclusion)
        return Action::INCLUDE;
    if (val <= exclusion)
        return Action::EXCLUDE;
    return Action::REMOVE;
}

CutoffGenerator::CutoffGenerator() {}
Cutoff CutoffGenerator::operator[](Feature const& f) const {
    return get_cutoff(f);
}

PercentCG::FeatureGroup::FeatureGroup(std::string const& spec) {
    try {
        std::vector<std::string> parts = split<' '>(spec);
        if (parts.size() != 2) {
            throw MyErr("improper specifier");
        }

        try {
            pattern = std::regex(parts[0], std::regex_constants::extended);
        }
        catch (std::regex_error const& err) {
            throw MyErr("invalid REGEX pattern: '", parts[0], "'");
        }

        cutoff = Cutoff(parts[1]);
    }
    catch (MyErr const& err) {
        throw err + MyErr("in registration of cutoff group with spec: '",
               spec, "'");
    }
}

PercentCG::PercentCG(Cutoff def, std::vector<std::string> const& specs) {
    default_cutoff = def;
    for (std::string const& spec : specs) {
        groups.push_back(FeatureGroup(spec));
    }
}
Cutoff PercentCG::get_cutoff(Feature const& f) const {
    for (FeatureGroup const& grp : groups) {
        if (std::regex_match(f.name, grp.pattern)) {
            return grp.cutoff;
        }
    }
    return default_cutoff;
}

CoveredInterval::CoveredInterval() : CoveredInterval(0, 1) {}
CoveredInterval::CoveredInterval(uint64_t s, uint64_t e) : Interval(s, e) {}
void CoveredInterval::absorb(Interval const& ivl) {
    if (ivl < *this || *this < ivl) { // ==> the entire interval is out of range
        return;
    }
    // calculate the coordinates of the new mega-interval created by this
    // interval combining any number of intervals together
    uint64_t min_start = std::max(start, ivl.start); // clip if extends beyond bounds
    uint64_t max_end = std::min(end, ivl.end); // clip if extends beyond bounds
    auto eqr = covered.equal_range(ivl);
    for (auto it = eqr.first; it != eqr.second; it++) {
        min_start = std::min(min_start, it->start);
        max_end = std::max(max_end, it->end);
    }

    // get rid of the old intervals (that have just been combined)
    covered.erase(eqr.first, eqr.second);
    // add back our new, combined interval
    covered.emplace(min_start, max_end);
}
uint64_t CoveredInterval::get_raw_coverage() const {
    uint64_t ret = 0;
    for (Interval const& ivl : covered) {
        ret += ivl.end-ivl.start;
    }
    return ret;
}
double CoveredInterval::get_percent_coverage() const {
    uint64_t sz = end-start;
    return 100.0*((double)get_raw_coverage())/((double)sz);
}

IsoformGenerator::Isoform::operator std::string() const {
    return repr;
}
IsoformGenerator::Isoform::operator bool() const {
    return good;
}

IsoformGenerator::IsoformGenerator(Annotation const& annot,
        CutoffGenerator const& cgen) : anno(annot) {
    // build our template isoform; fill in cutoff table
    for (Feature const& f : anno) {
        blank_iso[f.name] = CoveredInterval(f.start, f.end);
        cutoffs[f.name] = cgen.get_cutoff(f);
    }
}
IsoformGenerator::Isoform IsoformGenerator::generate(Alignment const& alig) {
    Isoform ret = blank_iso;
    for (Interval const& ivl : alig) {
        for (Feature const& f : anno.get_overlapping(ivl)) {
            ret[f.name].absorb(ivl);
        }
    }

    // serialize
    std::vector<std::string> included;
    for (auto const& p : ret) {
        Cutoff::Action adv = cutoffs[p.first].advise(
                ret[p.first].get_percent_coverage());
        if (adv == Cutoff::Action::INCLUDE) {
            included.push_back(p.first);
        }
        else if (adv == Cutoff::Action::REMOVE) {
            ret.good = false;
            return ret;
        }
    }
    if (included.empty()) {
        ret.good = false;
        return ret;
    }

    std::sort(included.begin(), included.end());
    for (int i = 0; i < included.size(); i++) {
        if (i != 0) {
            ret.repr += ",";
        }
        ret.repr += included[i];
    }
    ret.good = true;
    return ret;
}

BinCounter::BinCounter() : BinCounter(0, 1, 1) {}
BinCounter::BinCounter(double minval, double maxval, unsigned int n) {
    min = minval;
    max = maxval;
    if (max <= min) {
        throw MyErr("invalid bin range: (", min, ", ", max, ")");
    }

    binsz = (max-min)/((double)n);
    resize(n, 0);
}
void BinCounter::add(double val) {
    if (val < min || val > max) {
        throw MyErr("value out of bin range (", min, ", ", max, "): ", val);
    }
    unsigned int bin = (val-min)/binsz; // floor/truncate is correct
    if (bin == size()) {
        bin--;
    }
    at(bin)++;
}
BinCounter operator+(BinCounter const& a, BinCounter const& b) {
    if (a.min != b.min || a.max != b.max || a.size() != b.size()) {
        throw MyErr("incompatible bins: [(", a.min, ", ", b.max, "), n=", a.size(),
                "] and [(", b.min, ", ", b.max, "), n=", b.size(), "]");
    }
    BinCounter ret = a;
    for (int i = 0; i < b.size(); i++) {
        ret[i] += b[i];
    }
    return ret;
}
std::ostream& operator<<(std::ostream& out, BinCounter const& bc) {
    double lower, upper;
    for (int i = 0; i < bc.size(); i++) {
        lower = ((double)i)*bc.binsz;
        upper = ((double)(i+1))*bc.binsz;
        out << lower << "_to_" << upper << "\t" << bc[i] << std::endl;
    }
    return out;
}

CovDistCounter::CovDistCounter(Annotation const& anno) {
    for (Feature const& feat : anno) {
        dists[feat.name] = BinCounter(0.0, 100.0, 5);
    }
}
void CovDistCounter::count(IsoformGenerator::Isoform const& iso) {
    for (auto const& p : iso) {
        dists[p.first].add(p.second.get_percent_coverage());
    }
}
void CovDistCounter::save(std::string const& dir) const {
    for (auto const& p : dists) {
        std::string fpath = dir + "/" + p.first + ".tsv";
        std::ofstream fout(fpath);
        if (!fout) {
            throw MyErr("could not open file: ", fpath);
        }
        fout << p.second;
        fout.close();
    }
}

int main(int argc, char **argv) {
    try {
        Config cfg(argc, argv);
        PercentCG cgen(Cutoff(cfg.default_cutoff), cfg.cutoff_groups);

        Annotation anno(cfg.path_anno);

        IsoformGenerator isogen(anno, cgen);
        CovDistCounter dist_counter(anno);

        // process all reads
        std::ifstream alig_fin(cfg.path_alig);
        if (!alig_fin) {
            throw MyErr("could not open alignment file: ", cfg.path_alig);
        }
        std::string alig_line;
        std::unordered_map<std::string, unsigned int> iso_table;
        while (std::getline(alig_fin, alig_line), alig_fin) {
            if (alig_line.empty()) {
                continue;
            }
            Alignment alig(alig_line);
            IsoformGenerator::Isoform iso = isogen.generate(alig);

            // add to the table
            if (iso) {
                if (!iso_table.count(iso)) {
                    iso_table[iso] = 0;
                }
                iso_table[iso]++;
            }

            // save in distribution
            dist_counter.count(iso);
        }

        // output isoform table
        std::vector<std::pair<std::string, unsigned int>> iso_table_sorted(
                iso_table.begin(), iso_table.end());
        std::sort(iso_table_sorted.begin(), iso_table_sorted.end(),
                [](auto const& a, auto const& b){return a.second > b.second;});
        for (auto const& p : iso_table_sorted) {
            std::cout << p.first << '\t' << p.second << std::endl;
        }

        // output covdist
        if (!cfg.path_distdir.empty()) {
            dist_counter.save(cfg.path_distdir);
        }
    }
    catch (MyErr const& err) {
        std::cerr << "E: " << err << std::endl;
    }
}
