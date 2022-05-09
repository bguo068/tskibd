#include <algorithm>
#include <tskit.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <err.h>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>

using namespace std;

void
tsk_check_error(int val, int line_no)
{
    if (val < 0) {
        fprintf(stderr, "line %d: %s", line_no, tsk_strerror(val));
        exit(EXIT_FAILURE);
    }
}

struct IBDSegments {
    uint32_t id1;
    uint32_t id2;
    uint32_t start;
    uint32_t end;
    uint32_t mrca;
    double tmrca;
    uint8_t mut;
};

template <class T> class LowerTriangularMatrix
{
    vector<T> __vec; // convert float to int16_t by round(10xcM);
    uint32_t d;      // d = nsam if is not hap_pair matrix; otherwise d = 2 * nsam;

  public:
    // the constructor will initialize the array element to 0
    LowerTriangularMatrix(uint32_t dimension) : d(dimension)
    {
        size_t num_element = d;
        num_element = num_element * (num_element - 1) / 2;
        __vec.resize(num_element, 0);
    }
    LowerTriangularMatrix() = delete;

    void
    resize(uint32_t dimension, T val = 0)
    {
        size_t num_element = dimension;
        d = dimension;
        num_element = num_element * (num_element - 1) / 2;
        __vec.resize(num_element, 0);
    }

    // Lower trangular  row > col;
    //  -  -  -  -  -
    //  0  -  -  -  -
    //  1  2  -  -  -
    //  3  4  5  -  -
    //  6  7  8  9  -
    //
    size_t
    get_array_size()
    {
        size_t sz = d;
        return sz * (sz - 1) / 2;
    }

    size_t
    get_arr_index(size_t row, size_t col)
    {
        return row * (row - 1) / 2 + col;
    }

    uint32_t
    get_dimension()
    {
        return d;
    }

    void
    get_matrix_index(size_t arr_index, uint32_t &out_row, uint32_t &out_col)
    {
        // r* (r - 1) / 2 + c = y;
        // 0 <= c <= r-1
        // sqrt(2*y + 0.25) + 0.5 >= r >= sqrt(2*y + 2.25) - 0.5
        //

        out_row = sqrtl(2 * arr_index + 0.25) + 0.5;
        size_t tmp = out_row;
        out_col = arr_index - tmp * (tmp - 1) / 2;
    }

    T &
    at(size_t row, size_t col)
    {
        return __vec[row * (row - 1) / 2 + col];
    }

    T &
    at(size_t arr_index)
    {
        return __vec[arr_index];
    }

    static void
    test()
    {
        LowerTriangularMatrix<int> m(5UL);
        size_t num_element = m.get_array_size();
        for (size_t i = 0; i < num_element; i++) {
            m.at(i) = i;
            uint32_t r, c;
            m.get_matrix_index(i, r, c);
            assert(m.at(i) == m.at(r, c));
        }
        for (uint32_t i = 0; i < m.get_dimension(); i++) {
            for (uint32_t j = 0; j < m.get_dimension(); j++) {
                if (i > j)
                    cout << m.at(i, j) << '\t';
                else
                    cout << "-\t";
            }
            cout << '\n';
        }
    }
};

class Args
{
  public:
    int chrom;
    long bp_per_cm;
    long time_scale_factor;
    long true_IBD_sampling_window;
    bool hom_het;
    double delta; // tolerance when comparing ibd tmrca with segment tmrca.
    long min_length_bp;
    double mincm;
    // inferred data
    string treeseq_fn;
    int nsam;
    long nsegsites;
    int tree_sample_nodes_label_start_num; // macs starts with 0; msms starts
                                           // with 1;
    vector<long> positions;
    string haplotypes;
    vector<IBDSegments> ibdseg_arr;

    Args(int argc, char *argv[])
    {
        string u;
        u = "Usage : \n"
            "tsk_trueibd 1 1000000 10000 2.0 xxx.trees \n"
            "\n"
            "Positional parameters: \n"
            "1. chromN:                   Chromosome number (int). Use for \n"
            "                               name output file like 1.vcf, 1.ibd\n"
            "                               1.map file. It is also used as    \n"
            "                               chromosome contig name in vcf file\n"
            "2. bp_per_cm:                1000000 for human; 15000 for p.f.\n"
            "3. sampling_window:          1: check all trees; 1000: will check \n"
            "                               trees covers 0, 1000, 2000 \n"
            " 	                            suggestions: human->10000, pf->150; \n"
            "4. mincm:                    min length of ibd segment to keep\n\n"
            "5. treeseq_file:             tree sequence file that has finished \n"
            "                               coalescent\n\n";

        if (argc != 6) {
            cerr << u;
            exit(-1);
        } else {
            chrom = atoi(argv[1]);
            bp_per_cm = strtol(argv[2], NULL, 10);
            true_IBD_sampling_window = strtol(argv[3], NULL, 10);
            if (true_IBD_sampling_window == 0) {
                fprintf(stderr, "true_IBD_sampling_window parameter error!\n");
                exit(1);
            }
            mincm = strtod(argv[4], NULL);
            treeseq_fn = argv[5];
        }
        min_length_bp = mincm * bp_per_cm;
        nsam = -1;
    }

    void
    print(ostream &os)
    {
        cout << "====================================\n";
        cout << "chromN:               " << chrom << '\n';
        cout << "bp_per_cm:            " << bp_per_cm << '\n';
        cout << "time_scale_factor:    " << time_scale_factor << '\n';
        cout << "true_ibd_sample_win:  " << true_IBD_sampling_window << '\n';
        cout << "hom_het:              " << hom_het << '\n';
        if (nsam != -1)
            cout << "nsam:                 " << nsam << '\n';
    }
};

bool
operator<(const IBDSegments &seg1, const IBDSegments &seg2)
{
    return tie(seg1.id1, seg2.id2, seg1.start, seg1.end)
           < tie(seg2.id1, seg2.id2, seg2.start, seg2.end);
}

class Tsk
{
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_table_collection_t tables;
    int tree_iter;
    uint32_t tree_counter;
    uint32_t sampled_tree_counter;
    vector<tsk_id_t> stack;
    vector<tsk_id_t> leaves1;
    vector<tsk_id_t> leaves2;
    vector<tsk_id_t> internal_nodes;
    vector<tsk_id_t> children;
    // each element is the running ibd start bp value of a pair of haplotypes
    LowerTriangularMatrix<uint32_t> ibd_start_mat;
    // each element is the running ibd end bp value of a pair of haplotypes
    LowerTriangularMatrix<uint32_t> ibd_end_mat;
    // each element is mrca of a pair of haplotypes at the previous tree/ibd
    // segment time
    LowerTriangularMatrix<tsk_id_t> ibd_mrca_mat;
    // whether the ibd MRCA carry mutation
    LowerTriangularMatrix<uint8_t> ibd_mut_mat;
    // each element is mrca of a pair of haplotypes at the current tree
    LowerTriangularMatrix<tsk_id_t> cur_mrca_mat;
    LowerTriangularMatrix<uint8_t> cur_mut_mat;
    bool tree_contains_mutation;

  public:
    Tsk(string ts_fn)
        : ibd_start_mat(LowerTriangularMatrix<uint32_t>(1)),
          ibd_end_mat(LowerTriangularMatrix<uint32_t>(1)),
          ibd_mrca_mat(LowerTriangularMatrix<tsk_id_t>(1)),
          ibd_mut_mat(LowerTriangularMatrix<uint8_t>(1)),
          cur_mrca_mat(LowerTriangularMatrix<tsk_id_t>(1)),
          cur_mut_mat(LowerTriangularMatrix<uint8_t>(1))
    {
        int ret;
        ret = tsk_table_collection_load(&tables, ts_fn.c_str(), 0);
        tsk_check_error(ret, __LINE__);

        ret = tsk_table_collection_build_index(&tables, 0);
        tsk_check_error(ret, __LINE__);

        ret = tsk_treeseq_init(&ts, &tables, 0);
        tsk_check_error(ret, __LINE__);

        ret = tsk_tree_init(&tree, &ts, 0);

        tree_counter = 0;
        sampled_tree_counter = 0;

        ibd_start_mat.resize(ts.num_samples);
        ibd_end_mat.resize(ts.num_samples);
        ibd_mut_mat.resize(ts.num_samples);
        ibd_mrca_mat.resize(ts.num_samples);
        cur_mrca_mat.resize(ts.num_samples);
        cur_mut_mat.resize(ts.num_samples);
    }

    /*
        return 0 when reach the end
        return 1 when succefully move to the next tree
    */
    int
    iter_tree()
    {
        if (tree_counter == 0) {
            tree_iter = tsk_tree_first(&tree);
            tsk_check_error(tree_iter, __LINE__);
        } else {
            tree_iter = tsk_tree_next(&tree);
            tsk_check_error(tree_iter, __LINE__);
        }
        if (tree_iter == 1) {
            tree_counter++;
            assert(tsk_tree_get_num_roots(&tree) == 1);
        }

        tree_contains_mutation = is_mutation_in_tree_interval();

        return tree_iter;
    }

    bool
    is_multiple_in_tree_interval(uint32_t sampling_window_bp)
    {
        uint32_t left = lround(tree.interval.left);
        uint32_t right = lround(tree.interval.right);
        uint32_t multiple = right / sampling_window_bp * sampling_window_bp;
        bool isin = true;
        if (right == multiple || left > multiple) {
            isin = false;
        } else {
            sampled_tree_counter++;
        }
        return isin;
    }

    bool
    is_mutation_in_tree_interval()
    {
        bool yesno = false;
        for (tsk_size_t i = 0; i < tables.sites.num_rows; i++) {
            if (tables.sites.position[i] >= tree.interval.left
                && tables.sites.position[i] < tree.interval.right) {
                yesno = true;
                break;
            }
        }
        return yesno;
    }

    // private call once
    bool
    do_node_carry_mutation(tsk_id_t u)
    {
        bool yesno = false;
        tsk_id_t v;
        for (tsk_size_t i = 0; i < tables.mutations.num_rows; i++) {
            v = tables.mutations.node[i];
            if (tsk_tree_is_descendant(&tree, u, v)) {
                yesno = true;
                break;
            }
        }
        return yesno;
    }

    // public can call many times
    bool
    do_tree_contains_mutation()
    {
        return tree_contains_mutation;
    }
    void
    traverse_tree()
    {
        tsk_id_t u, v;
        stack.resize(0);

        stack.push_back(tree.virtual_root);
        while (stack.size() > 0) {
            u = stack.back();
            stack.pop_back();
            cout << "Visit node " << u << '\t';
            /* Put nodes on the stack right-to-left, so we visit in left-to-right */
            for (v = tree.right_child[u]; v != TSK_NULL; v = tree.left_sib[v]) {
                stack.push_back(v);
            }
        }
    }

    void
    find_leaves(tsk_id_t u, vector<tsk_id_t> &leaves_out)
    {
        tsk_id_t v;
        stack.resize(0);
        leaves_out.resize(0);

        stack.push_back(u);
        while (stack.size() > 0) {
            u = stack.back();
            stack.pop_back();
            // if u is a leave node
            if (tree.right_child[u] == TSK_NULL) {
                leaves_out.push_back(u);
            }
            /* Put nodes on the stack right-to-left, so we visit in left-to-right */
            for (v = tree.right_child[u]; v != TSK_NULL; v = tree.left_sib[v]) {
                stack.push_back(v);
            }
        }
    }

    void
    find_children(tsk_id_t p)
    {
        tsk_id_t v, u;
        children.resize(0);
        u = tree.left_child[p];
        v = tree.right_child[p];
        for (; u != v; u = tree.right_sib[u]) {
            children.push_back(u);
        }
        children.push_back(u);
    }

    void
    find_internal_nodes()
    {
        tsk_id_t u, v;
        stack.resize(0);
        internal_nodes.resize(0);

        // find the real root
        u = tree.virtual_root;
        u = tree.right_child[u];

        stack.push_back(u);
        while (stack.size() > 0) {
            u = stack.back();
            stack.pop_back();
            // if u is a leave node
            if (tree.right_child[u] != TSK_NULL) {
                internal_nodes.push_back(u);
            }
            /* Put nodes on the stack right-to-left, so we visit in left-to-right */
            for (v = tree.right_child[u]; v != TSK_NULL; v = tree.left_sib[v]) {
                stack.push_back(v);
            }
        }
    }
    void
    calculate_tmrca_mat(Args &args)
    {
        tsk_id_t child1, child2;
        find_internal_nodes();
        bool is_p_with_mutation = false;

        for (tsk_id_t p : internal_nodes) {
            find_children(p);
            is_p_with_mutation
                = is_mutation_in_tree_interval() && do_node_carry_mutation(p);

            for (size_t i = 1; i < children.size(); i++) {
                child1 = children[i];
                find_leaves(child1, leaves1);

                for (size_t j = 0; j < i; j++) {
                    child2 = children[j];
                    find_leaves(child2, leaves2);

                    for (tsk_id_t id1 : leaves1) {
                        for (tsk_id_t id2 : leaves2) {
                            tsk_id_t id1_copy = id1;
                            if (id1_copy < id2) {
                                swap(id1_copy, id2);
                            }
                            cur_mrca_mat.at(id1_copy, id2) = p;
                            cur_mut_mat.at(id1_copy, id2) = is_p_with_mutation;
                        }
                    }
                }
            }
        }
    }
    // output ibd and update matrix related to ibd segments
    void
    find_and_update_ibd(Args &args)
    {
        uint32_t tree_interval_start_bp = lround(tree.interval.left);
        uint32_t tree_interval_end_bp = lround(tree.interval.right);
        size_t num_pairs = ibd_start_mat.get_array_size();
        for (size_t pair = 0; pair < num_pairs; pair++) {
            // use reference to be update
            auto &ibd_mrca = ibd_mrca_mat.at(pair);
            auto &ibd_end = ibd_end_mat.at(pair);
            auto &ibd_start = ibd_start_mat.at(pair);
            auto &ibd_mut = ibd_mut_mat.at(pair);
            auto cur_mrca = cur_mrca_mat.at(pair);
            auto cur_mut = cur_mut_mat.at(pair);

            // ibd continues
            if (ibd_mrca == cur_mrca) {
                ibd_end = tree_interval_end_bp;
                // if seeing mutation, set ibd_mut to true; if not, keep the
                // same as previous. So ibd_mut will be true if it see at least
                // one mutation across marginal trees until it breaks
                if (cur_mut) {
                    ibd_mut = 1;
                }

            }
            // ibd breaks
            else {
                // if the ibd segments is long enough output ibd
                if (tree_interval_start_bp - ibd_start > args.min_length_bp) {
                    IBDSegments seg;
                    // find id1 and id2 using matrix index
                    ibd_mrca_mat.get_matrix_index(pair, seg.id1, seg.id2);
                    seg.start = ibd_start;
                    /*
                    need to use cur tree left as the seg end. Otherwise, AS-IBD coords
                    never to touch each other as we skip trees using the sampling
                    windows.
                    */
                    // seg.end = ibd_end;
                    seg.end = tree_interval_start_bp;
                    seg.mrca = ibd_mrca;
                    seg.tmrca = tables.nodes.time[ibd_mrca];
                    seg.mut = ibd_mut;
                    args.ibdseg_arr.push_back(seg);
                }
                // update ibd matrix elements
                ibd_start = tree_interval_start_bp; // should start from cur interval
                                                    // start instead of cur interval end
                ibd_end = tree_interval_end_bp;
                ibd_mrca = cur_mrca;
                ibd_mut = cur_mut; // reset mutation status
            }
        }
    }

    // when reaching the last tree, there migh be some IBD continuing to the end
    // of the chromsome we need to output the ibd
    void
    find_hanging_ibd(Args &args)
    {
        size_t num_pairs = ibd_start_mat.get_array_size();
        for (size_t pair = 0; pair < num_pairs; pair++) {
            // use reference to be update
            auto &ibd_mrca = ibd_mrca_mat.at(pair);
            auto &ibd_end = ibd_end_mat.at(pair);
            auto &ibd_start = ibd_start_mat.at(pair);
            auto &ibd_mut = ibd_mut_mat.at(pair);
            // if the ibd segments is long enough output ibd
            if (ibd_end - ibd_start > args.min_length_bp) {
                IBDSegments seg;
                // find id1 and id2 using matrix index
                ibd_mrca_mat.get_matrix_index(pair, seg.id1, seg.id2);
                seg.start = ibd_start;
                seg.end = ibd_end;
                seg.mrca = ibd_mrca;
                seg.tmrca = tables.nodes.time[ibd_mrca];
                seg.mut = ibd_mut;
                args.ibdseg_arr.push_back(seg);
            }
        }
    }

    uint32_t
    get_tree_counter() const
    {
        return tree_counter;
    }
    uint32_t
    get_sampled_tree_counter() const
    {
        return sampled_tree_counter;
    }
    uint32_t
    get_tree_interval_end_bp() const
    {
        return lround(tree.interval.right);
    }

    double
    get_seqlen() const
    {
        return tables.sequence_length;
    }
    ~Tsk()
    {
        tsk_tree_free(&tree);
        tsk_table_collection_free(&tables);
        tsk_treeseq_free(&ts);
    }
};

class FileWriter
{
    ofstream &ofs_log;
    ofstream &ofs_ibd;
    ofstream &ofs_map;

  public:
    FileWriter(ofstream &log, ofstream &ibd, ofstream &map)
        : ofs_log(log), ofs_ibd(ibd), ofs_map(map)
    {
    }

    void
    write_log(const Tsk &tsk, const Args &args, ostream &os)
    {
        os << "\ttree: " << tsk.get_tree_counter()
           << "\ttree parsed: " << tsk.get_sampled_tree_counter()
           << "\tibd no.: " << args.ibdseg_arr.size() << "\tchromsome parsed: "
           << 1.0 * tsk.get_tree_interval_end_bp() / tsk.get_seqlen() << '\n';
    }
    void
    write_log(const Tsk &tsk, const Args &args)
    {
        ofs_log << "\ttree: " << tsk.get_tree_counter()
                << "\ttree parsed: " << tsk.get_sampled_tree_counter()
                << "\tibd no.: " << args.ibdseg_arr.size() << "\tchromsome parsed: "
                << 1.0 * tsk.get_tree_interval_end_bp() / tsk.get_seqlen() << '\n';
    }
    void
    write_args_to_log(const Args &args)
    {
        // TODO
        ofs_log << "====================================\n";
        ofs_log << "chromN:               " << args.chrom << '\n';
        ofs_log << "bp_per_cm:            " << args.bp_per_cm << '\n';
        ofs_log << "true_ibd_sample_win:  " << args.true_IBD_sampling_window << '\n';
        ofs_log << "mincm:                " << args.mincm << '\n';
        ofs_log << "treeseq file:         " << args.treeseq_fn << '\n';
        if (args.nsam != -1)
            ofs_log << "nsam:                 " << args.nsam << '\n';
    }

    void
    write_map(Args &args, const Tsk &tsk)
    {
        ofs_map << args.chrom << " . " << 1.0 / args.bp_per_cm << " " << 1 << '\n';
        ofs_map << args.chrom << " . " << tsk.get_seqlen() / args.bp_per_cm << " "
                << (uint32_t) tsk.get_seqlen() << '\n';
    }

    void
    write_ibd(Args &args, bool raw = false)
    {
        if (!raw) {
            for (auto seg : args.ibdseg_arr) {
                int sample1, sample2, hap1, hap2;
                sample1 = seg.id1 / 2;
                hap1 = seg.id1 % 2;
                sample2 = seg.id2 / 2;
                hap2 = seg.id1 % 2;
                // ignore homozygous by descent (.hbd)
                if (sample1 == sample2)
                    continue;
                // ensure the order of samples
                if (sample1 > sample2) {
                    swap(sample1, sample2);
                    swap(hap1, hap2);
                }
                ofs_ibd << sample1 << '\t';
                ofs_ibd << hap1 << '\t';
                ofs_ibd << sample2 << '\t';
                ofs_ibd << hap2 << '\t';
                ofs_ibd << args.chrom << '\t';
                ofs_ibd << seg.start << '\t';
                ofs_ibd << seg.end << '\t';
                ofs_ibd << (seg.end - seg.start) * 1.0 / args.bp_per_cm << '\t';
                ofs_ibd << seg.mrca << '\t';
                ofs_ibd << seg.tmrca << '\t';
                ofs_ibd << (uint32_t) seg.mut << '\n';
            }
        } else {
            ofs_ibd << "Id1\tId2\tStart\tEnd\tAncestor\tTmrca\tHasMutation\n";
            for (auto seg : args.ibdseg_arr) {
                ofs_ibd << seg.id1 << '\t';
                ofs_ibd << seg.id2 << '\t';
                ofs_ibd << seg.start << '\t';
                ofs_ibd << seg.end << '\t';
                ofs_ibd << seg.mrca << '\t';
                ofs_ibd << seg.tmrca << '\t';
                ofs_ibd << (uint32_t) seg.mut << '\n';
            }
        }
    }
};

int
main(int argc, char *argv[])
{
    // TreeParser::test();
    auto args = Args(argc, argv);

    // create files
    ofstream ofs_ibd(to_string(args.chrom) + ".ibd");
    ofstream ofs_map(to_string(args.chrom) + ".map");
    ofstream ofs_log(to_string(args.chrom) + ".log");
    auto writer = FileWriter(ofs_log, ofs_ibd, ofs_map);
    writer.write_args_to_log(args);

    Tsk tsk(args.treeseq_fn.c_str());
    while (tsk.iter_tree()) {
        if (tsk.is_multiple_in_tree_interval(args.true_IBD_sampling_window)
            || tsk.do_tree_contains_mutation()) {
            tsk.calculate_tmrca_mat(args);
            tsk.find_and_update_ibd(args);
            if (tsk.get_sampled_tree_counter() % 20 == 1) {
                writer.write_log(tsk, args);
                writer.write_log(tsk, args, cerr);
            }
        }
    }
    tsk.find_hanging_ibd(args);
    writer.write_ibd(args, true);

    writer.write_map(args, tsk);

    return 0;
}