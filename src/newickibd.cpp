#include <algorithm>
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
#include <vector>

// structure for a tree node
struct Node {
  int index;         // position in the node array
  int label;         // corresponding to an id number of present sample
  double len_above;  // time to the parent node; or -1 if it is the root node
  int parent;        // index to parent node in the nodeArr array
  int left;     // index to left child node in the nodeArr array, -1 indicate
                // sample node
  int right;    // index to right child node in the nodeArr array, -1 indicate
                // sample node
  double time;  // time to present
};

struct IBDSegments {
  uint32_t id1;
  uint32_t id2;
  long start;
  long end;
  double tmrca;
};

bool operator<(const IBDSegments &seg1, const IBDSegments &seg2) {
  return std::tie(seg1.id1, seg2.id2, seg1.start, seg1.end) <
         std::tie(seg2.id1, seg2.id2, seg2.start, seg2.end);
}

template <class T>
class LowerTriangularMatrix {
  std::vector<T> __vec;  // convert float to int16_t by round(10xcM);
  uint32_t d;  // d = nsam if is not hap_pair matrix; otherwise d = 2 * nsam;

 public:
  // the constructor will initialize the array element to 0
  LowerTriangularMatrix(uint32_t dimension) : d(dimension) {
    size_t num_element = d;
    num_element = num_element * (num_element - 1) / 2;
    __vec.resize(num_element, 0);
  }
  LowerTriangularMatrix() = delete;

  // Lower trangular  row > col;
  //  -  -  -  -  -
  //  0  -  -  -  -
  //  1  2  -  -  -
  //  3  4  5  -  -
  //  6  7  8  9  -
  //
  size_t get_array_size() {
    size_t sz = d;
    return sz * (sz - 1) / 2;
  }

  size_t get_arr_index(size_t row, size_t col) {
    return row * (row - 1) / 2 + col;
  }

  uint32_t get_dimension() { return d; }

  void get_matrix_index(size_t arr_index, uint32_t &out_row,
                        uint32_t &out_col) {
    // r* (r - 1) / 2 + c = y;
    // 0 <= c <= r-1
    // sqrt(2*y + 0.25) + 0.5 >= r >= sqrt(2*y + 2.25) - 0.5
    //

    out_row = sqrtl(2 * arr_index + 0.25) + 0.5;
    size_t tmp = out_row;
    out_col = arr_index - tmp * (tmp - 1) / 2;
  }

  T &at(size_t row, size_t col) { return __vec[row * (row - 1) / 2 + col]; }

  T &at(size_t arr_index) { return __vec[arr_index]; }

  static void test() {
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
          std::cout << m.at(i, j) << '\t';
        else
          std::cout << "-\t";
      }
      std::cout << '\n';
    }
  }
};
class Args {
 public:
  int chrom;
  long bp_per_cm;
  long time_scale_factor;
  long true_IBD_sampling_window;
  bool hom_het;
  double delta;  // tolerance when comparing ibd tmrca with segment tmrca.
  long min_length_bp;
  double mincm;
  // inferred data
  std::string simulation_command_str;
  long seqlen;
  int nsam;
  long nsegsites;
  int tree_sample_nodes_label_start_num;  // macs starts with 0; msms starts
                                          // with 1;
  std::vector<long> positions;
  std::string haplotypes;
  std::vector<IBDSegments> ibdseg_arr;

  Args(int argc, char *argv[]) {
    std::string u;
    u = "Usage 1: \n"
        "./macs 50 10000 -t 0.0004 -r 0.0004 -h 1000 -T 2>/dev/null \\\n"
        "  | ./msformatter | newick_tree 1 1000000 40000 10000 0 2.0\n"
        "\n"
        "Usage 2: \n"
        "java -Xmx10g -jar ./msms.jar 100 1 -t 4000 -r 4000 10000000 -T \\\n"
        "  | newick_tree 1 1000000 40000 10000 0 2.0\n"
        "\n"
        "Positional parameters: \n"
        "1. chromN:                     Chromosome number (int). Use for \n"
        "                               name output file like 1.vcf, 1.ibd\n"
        "                               1.map file. It is also used as    \n"
        "                               chromosome contig name in vcf file\n"
        "2. bp_per_cm:                  1000000 for human; 15000 for p.f.\n"
        "3. time_scale_factor:          2N0 (macs); 4N0 (msms) \n"
        "4. true_IBD_sampling_window:   1: check all trees; 1000: will check \n"
        "                               trees covers 0, 1000, 2000 \n"
        " 	                            suggestions: human->10000, "
        "pf->150; \n"
        "5. hom_het:                    1: homozygous diploid; 0: heterozygous "
        "\n"
        "                               diploid when convert haplotypes into \n"
        "                               vcf.                                  "
        "\n"
        "6. mincm:                      min length of ibd segment to keep\n\n";

    if (argc != 7) {
      std::cerr << u;
      exit(-1);
    } else {
      chrom = atoi(argv[1]);
      bp_per_cm = strtol(argv[2], NULL, 10);
      time_scale_factor = (long)strtod(argv[3], NULL);
      true_IBD_sampling_window = strtol(argv[4], NULL, 10);
      if (true_IBD_sampling_window == 0) {
        fprintf(stderr, "true_IBD_sampling_window parameter error!\n");
        exit(1);
      }
      hom_het = atoi(argv[5]);
      mincm = strtod(argv[6], NULL);
    }
    delta = 0.5 / time_scale_factor;
    mincm = 2.0;
    min_length_bp = mincm * bp_per_cm;
    nsam = -1;
    seqlen = -1;
    nsegsites = -1;
    positions.reserve(100);
    haplotypes = "";
  }

  void print(std::ostream &os) {
    std::cout << "====================================\n";
    std::cout << "chromN:               " << chrom << '\n';
    std::cout << "bp_per_cm:            " << bp_per_cm << '\n';
    std::cout << "time_scale_factor:    " << time_scale_factor << '\n';
    std::cout << "true_ibd_sample_win:  " << true_IBD_sampling_window << '\n';
    std::cout << "hom_het:              " << hom_het << '\n';
    if (nsam != -1) std::cout << "nsam:                 " << nsam << '\n';
    if (seqlen != -1) std::cout << "seqlen:               " << seqlen << '\n';
    if (nsegsites != -1)
      std::cout << "nsegsites:            " << seqlen << '\n';

    if (positions.size() > 0) {
      std::cout << "positions:            " << '\n';
      std::copy(positions.begin(), positions.end(),
                std::ostream_iterator<long>(std::cout, ", "));
    }
    if (haplotypes.size() > 0) {
      std::cout << "haplotypes:           " << '\n';
      for (ssize_t i = 0; i < nsegsites; i++) {
        auto sv = std::string_view(haplotypes.c_str() + i * nsam, nsam);
        std::cout << sv << '\n';
      }
    }
  }
};

class MsLikeFileReader {
  std::istream &ifs_ms;
  std::string line;
  std::string mode;  // "cmd"
  long tree_counter;
  long sample_counter;

 public:
  MsLikeFileReader(std::istream &ms) : ifs_ms(ms) {
    tree_counter = 0;
    sample_counter = 0;
  }
  void parse_meta(Args &args) {
    mode = "cmd";

    // tokenize
    std::getline(ifs_ms, line);
    auto tokens = tokenize_string(line);

    // condition on simulators
    if (tokens[0].substr(tokens[0].size() - 2) == "ms") {
      // msms
      args.simulation_command_str = line;
      int nsam = std::stoi(tokens[1]);
      int nrep = std::stoi(tokens[2]);
      uint32_t idx_r = 0;
      long seqlen = 0;
      if (nrep != 1) {
        std::cerr << "nrep should be 1 but it is not in the mslike command\n";
        exit(1);
      }
      while (tokens[idx_r] != "-r") idx_r += 1;
      if (tokens[idx_r] == "-r" && tokens.size() > idx_r + 2) {
        seqlen = std::stol(tokens[idx_r + 2]);
      } else {
        std::cerr << "Error: can't parse number of cut sites from msLike "
                     "results\n";
        exit(1);
      }
      args.nsam = nsam;
      args.seqlen = seqlen;
      args.tree_sample_nodes_label_start_num =
          1;  // msms node label starts with 1

    } else if (tokens[0].substr(tokens[0].size() - 4) == "macs") {
      args.simulation_command_str = line;
      if (tokens.size() > 3) {
        int nsam = std::stoi(tokens[1]);
        int seqlen = std::stol(tokens[2]);
        if (nsam <= 0 || seqlen <= 0) {
          std::cerr << "nsam or seqlen not parsed\n";
          exit(1);
        }
        args.nsam = nsam;
        args.seqlen = seqlen;
        // msms node label starts with 1
        args.tree_sample_nodes_label_start_num = 0;
      } else {
        std::cerr << "cmd not recognized\n";
        exit(1);
      }

    } else {
      std::cerr << "cmd not recognized\n";
      exit(1);
    }
  }

  std::vector<std::string> tokenize_string(const std::string &str) {
    std::istringstream iss(str);
    std::string tok;
    std::vector<std::string> tokens;
    while (std::getline(iss, tok, ' ')) tokens.push_back(tok);
    return tokens;
  }

  void pass_util_tree_lines() {
    mode = "pass_util_tree_lines";

    // tokenize
    while (std::getline(ifs_ms, line) && line.substr(0, 2) != "//")
      ;
    if (line.substr(0, 2) != "//") {
      std::cerr << "expect passing over until line starts with //\n";
      exit(1);
    }
  }

  bool iter_tree_line() {
    mode = "get_tree_line";
    if (!std::getline(ifs_ms, line))
      mode = "end";
    else if (line.size() == 0 || line[0] != '[')
      mode = "";
    tree_counter++;
    return mode == "get_tree_line";
  }

  std::string &get_line_buffer() { return line; }

  void parse_segsites(Args &args, bool start_with_buffer = true) {
    mode = "segsites";
    // iterate line until the segsites line
    while (start_with_buffer || std::getline(ifs_ms, line)) {
      start_with_buffer = false;
      if (line.substr(0, 8) == "segsites") break;
    }
    if (line.size() < 8 || line.substr(0, 8) != "segsites") {
      std::cerr << "Can't find segsites lines\n";
      exit(1);
    }
    auto idx = line.find_first_of(' ');
    args.nsegsites = strtol(line.c_str() + idx + 1, NULL, 10);
  }

  // scale up positions into bp; also address position collision.
  void parse_positions(Args &args) {
    mode = "positions";
    std::getline(ifs_ms, line);
    std::istringstream iss(line);
    std::string tok;
    while (std::getline(iss, tok, ' ')) {
      if (tok.size() > 0 && (tok[0] == '0' || tok[0] == 1)) {
        auto pos01 = std::lround(std::stod(tok) * args.seqlen);
        if (args.positions.size() > 0 && pos01 <= args.positions.back())
          pos01 = args.positions.back() + 1;
        args.positions.push_back(pos01);
      }
    }
  }
  void parse_haplotypes(Args &args) {
    mode = "haplotypes";
    args.haplotypes.resize(args.nsam * args.nsegsites);
    while (std::getline(ifs_ms, line) && line != "") {
      // Note: the hanotypes are transpose to be able to convert to vcf file
      for (ssize_t i = 0; i < args.nsegsites; i++)
        args.haplotypes[i * args.nsam + sample_counter] = line[i];
      sample_counter++;
    }
  }
};

class TestData {
 public:
  std::string tree_str_nsam100;
  std::string tree_str_nsam5;
  std::string ms_output_str;
  TestData() {
    tree_str_nsam5 = "[100](((3:0.1,4:0.1):0.1,5:0.2):0.2,(7:0.3, 8:0.3):0.1);";
    ms_output_str =
        "ms 5 1 -t 40 -r 40 100000 -T  [3.2rc Build:162]\n"
        "0x2ee1e7820a517c09\n"
        "\n"
        "//\n"
        "[7529](4:0.14095,((3:0.06957,(2:0.01382,1:0.01382):0.05576):0.00522,5:"
        "0."
        "07479):0.06617);\n"
        "[382](3:0.36020,(4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):"
        "0."
        "06617):0.21924);\n"
        "[144](3:0.84656,(4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):"
        "0."
        "06617):0.70561);\n"
        "[2399](3:0.69732,(4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):"
        "0."
        "06617):0.55636);\n"
        "[1081](3:0.69732,(4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):"
        "0."
        "06617):0.55636);\n"
        "[242](3:0.69732,(4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):"
        "0."
        "06617):0.55636);\n"
        "[487]((4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):0.06617):0."
        "46698,3:0.60793);\n"
        "[542]((4:0.14095,((2:0.01382,1:0.01382):0.06097,5:0.07479):0.06617):0."
        "18547,3:0.32642);\n"
        "[3778](4:0.33899,(((2:0.01382,1:0.01382):0.06097,5:0.07479):0.25164,3:"
        "0."
        "32642):0.01257);\n"
        "[355](4:0.33899,(((2:0.01382,1:0.01382):0.06097,5:0.07479):0.25164,3:"
        "0."
        "32642):0.01257);\n"
        "[4653]((4:0.31929,((2:0.01382,1:0.01382):0.06097,5:0.07479):0.24450):"
        "0."
        "01970,3:0.33899);\n"
        "[352]((4:0.31929,((2:0.01382,1:0.01382):0.12714,5:0.14095):0.17833):0."
        "01970,3:0.33899);\n"
        "[991]((((2:0.01382,1:0.01382):0.12714,5:0.14095):0.19804,3:0.33899):0."
        "03538,4:0.37437);\n"
        "[2008](((3:0.11447,(2:0.01382,1:0.01382):0.10065):0.02649,5:0.14095):"
        "0."
        "23341,4:0.37437);\n"
        "[843](((3:0.11447,(2:0.01382,1:0.01382):0.10065):0.02649,5:0.14095):0."
        "65321,4:0.79416);\n"
        "[3300](4:0.63992,((3:0.11447,(2:0.01382,1:0.01382):0.10065):0.02649,5:"
        "0."
        "14095):0.49897);\n"
        "[61](4:0.63992,((3:0.11447,(2:0.01382,1:0.01382):0.10065):0.02649,5:0."
        "14095):0.49897);\n"
        "[491](4:0.63992,((3:0.11447,(2:0.01382,1:0.01382):0.10065):0.02649,5:"
        "0."
        "14095):0.49897);\n"
        "[784]((3:0.04214,4:0.04214):0.59778,((2:0.01382,1:0.01382):0.12714,5:"
        "0."
        "14095):0.49897);\n"
        "[1149]((3:0.04214,4:0.04214):0.59778,((2:0.01382,1:0.01382):0.12714,5:"
        "0."
        "14095):0.49897);\n"
        "[33067](((2:0.01382,1:0.01382):0.04408,(3:0.04214,4:0.04214):0.01575):"
        "0."
        "08306,5:0.14095);\n"
        "[5591]((2:0.05789,(3:0.04214,4:0.04214):0.01575):0.08306,(5:0.04902,1:"
        "0."
        "04902):0.09194);\n"
        "[722](1:1.78824,((2:0.05789,(3:0.04214,4:0.04214):0.01575):0.08306,5:"
        "0."
        "14095):1.64729);\n"
        "[1660](1:0.53261,((2:0.05789,(3:0.04214,4:0.04214):0.01575):0.08306,5:"
        "0."
        "14095):0.39165);\n"
        "[6799](((1:0.05165,2:0.05165):0.00624,(3:0.04214,4:0.04214):0.01575):"
        "0."
        "08306,5:0.14095);\n"
        "[7050]((1:0.05789,(3:0.04214,4:0.04214):0.01575):0.08306,(2:0.05404,5:"
        "0."
        "05404):0.08691);\n"
        "[2932]((1:0.05789,(3:0.04214,4:0.04214):0.01575):0.08306,(2:0.07479,5:"
        "0."
        "07479):0.06617);\n"
        "[6640]((1:0.05789,(3:0.03984,4:0.03984):0.01805):0.08306,(2:0.07479,5:"
        "0."
        "07479):0.06617);\n"
        "[3968]((2:0.07479,5:0.07479):0.04403,(1:0.05789,(3:0.03984,4:0.03984):"
        "0."
        "01805):0.06093);\n"
        "segsites: 22\n"
        "positions: 0.02426 0.04606 0.07566 0.10185 0.11728 0.13347 0.14277 "
        "0.17632 0.19185 0.26375 0.26637 0.29702 0.30515 0.39073 0.51731 "
        "0.53111 "
        "0.56699 0.73240 0.76334 0.84720 0.91126 0.99044 \n"
        "1010000100110000100001\n"
        "1010000100110000100010\n"
        "1001100000101100100000\n"
        "0110011111001100110100\n"
        "1010000100110011001000\n";
  }
};

class TreeParser {
  // a string representing a whole line of newick forat tree string
  std::string tree_str;
  long tree_span_bp;          // length in bp the tree spans
  long tree_interval_end_bp;  // current tree interval start position in bp
  int num_node_used;  // running counter for how many nodes has been processed
  int num_node;       // total number nodes allocated.
  std::vector<Node> node_arr;  // a vector of nodes
  long tree_counter;
  long sampled_tree_counter;
  // for finding leaves of a given ancestral node. Each element is the leave
  // node label;
  std::vector<int> leaves1;
  std::vector<int> leaves2;
  // each element is the running ibd start bp value of a pair of haplotypes
  LowerTriangularMatrix<long> ibd_start_mat;
  // each element is the running ibd end bp value of a pair of haplotypes
  LowerTriangularMatrix<long> ibd_end_mat;
  // each element is tmrca of a pair of haplotypes at the previous tree/ibd
  // segment time
  LowerTriangularMatrix<double> ibd_time_mat;
  // each element is tmrca of a pair of haplotypes at the current tree
  LowerTriangularMatrix<double> tmrca_mat;

 public:
  TreeParser(size_t nsam)
      : ibd_start_mat(LowerTriangularMatrix<long>(nsam)),
        ibd_end_mat(LowerTriangularMatrix<long>(nsam)),
        ibd_time_mat(LowerTriangularMatrix<double>(nsam)),
        tmrca_mat(LowerTriangularMatrix<double>(nsam)) {
    num_node = nsam * 2 - 1;
    node_arr.resize(num_node);
    leaves1.reserve(num_node);
    leaves2.reserve(num_node);
    tree_interval_end_bp = 0;
    sampled_tree_counter = 0;
    tree_counter = 0;
  }

  long get_tree_counter() const { return tree_counter; }

  long get_sampled_tree_counter() const { return sampled_tree_counter; }

  long get_tree_interval_end_bp() const { return tree_interval_end_bp; }

  static int get_nsam(std::string &tree_string) {
    // count number of many pairs of parenthesis
    int nsam = 0;
    for (auto &c : tree_string) {
      if (c == '(') nsam++;
    }
    nsam++;
    return nsam;
  }

  int parse_tree(std::string &in_tree_str, Args &args) {
    // store tree string
    tree_str = in_tree_str;

    // get the tree span
    auto s = tree_str.find_first_of('[') + 1;
    tree_span_bp = std::stol(tree_str.c_str() + s);
    tree_counter++;

    // update tree_interval_end_bp for the next tree after update ibd info
    tree_interval_end_bp += tree_span_bp;

    // condition on whether the tree interveal contains a check point for
    // tree tmrca comparison; if a mulplte of true_IBD_sampling_window is
    // not on the segment, skip parsing the tree using browning method

    long multiple = (tree_interval_end_bp / args.true_IBD_sampling_window);
    multiple *= args.true_IBD_sampling_window;

    if (multiple < tree_interval_end_bp - tree_span_bp ||
        multiple == tree_interval_end_bp) {
      return -1;
    }

    // range for the newick tree
    s = tree_str.find_first_of('(');
    auto e = tree_str.size();

    // initize num_node used
    num_node_used = 0;

    // pare parent = -1 to indicate the substree started at tree root
    const char *start = tree_str.c_str() + s;
    const char *end = tree_str.c_str() + e - 1;
    parse_subtree_recursive(start, end, -1);

    // cal node time
    calculate_node_time();

    // fill the tmrca
    calculate_tmrca_mat(args);

    sampled_tree_counter++;

    return 0;
  }

 private:
  /*
   * Recursive function to parse subtree due to a node.
   *
   * return the index of subtree root node in the node_arr
   */
  int parse_subtree_recursive(const char *start, const char *end, int parent) {
    // use a new node
    auto cur_node_idx = num_node_used;
    auto &cur_node = node_arr[cur_node_idx];

    // inc counter
    num_node_used++;

    // store parent and node array
    cur_node.parent = parent;
    cur_node.index = cur_node_idx;

    // 1. find the separator (; or :) from right to left
    const char *pSep = end;
    while (*pSep != ':' && *pSep != ';' && pSep >= start) pSep--;

    // once found the number on the right of the seperator is the branch length
    // above the current node
    if (*pSep == ';') {
      // root sample
      cur_node.len_above = 0;
    } else if (*pSep == ':')
      // non root sample
      cur_node.len_above = strtod(pSep + 1, NULL);

    // 2. find child separater or comma (,) in outmost ():
    //  The key is to ensure that there are extract one more parentheses on the
    //  right than on the left the intersted command, when looking from the
    //  right
    int diff_RLparentheses = 0;
    const char *pComma = pSep;
    while (!(*pComma == ',' && diff_RLparentheses == 1)) {
      if (*pComma == ')')
        diff_RLparentheses++;
      else if (*pComma == '(')
        diff_RLparentheses--;
      pComma--;
      if (pComma <= start) break;
    }
    // 3. if there is no comma, it means the node is not compund node and thus
    // is a current sample
    if (*pComma != ',') {
      cur_node.label = atoi(start);
      cur_node.left = -1;
      cur_node.right = -1;
    }
    // 4. otherwise it is compund node, then divide to left chid subtree and
    // right child subtree. Note: need to shift inside from as the outer
    // parenthesis/comma/colon is part of the current node.
    else {
      cur_node.left =
          parse_subtree_recursive(start + 1, pComma - 1, cur_node_idx);
      cur_node.right =
          parse_subtree_recursive(pComma + 1, pSep - 2, cur_node_idx);
      cur_node.label = -1;
    }
    return cur_node_idx;
  }

  void calculate_node_time() {
    double t;
    int i, j;
    for (i = 0; i < num_node; i++) {
      t = 0.0;
      for (j = node_arr[i].left; j != -1; j = node_arr[j].left)
        t += node_arr[j].len_above;
      node_arr[i].time = t;
    }

    // // debug
    // for (i = 0; i < num_node; i++) {
    //     if (node_arr[i].label != -1)
    //         assert(node_arr[i].time == 0);
    //     else
    //         assert(node_arr[i].time > 0);
    // }
  }

 public:
  int check_tree(const Args &arg) {
    // check labels
    ssize_t label_neg_counter = 0;
    ssize_t right_neg_counter = 0;
    ssize_t left_neg_counter = 0;
    ssize_t parent_neg_counter = 0;
    ssize_t time_zero_counter = 0;
    int num_failure = 0;

    for (int i = 0; i < num_node; i++) {
      auto n = node_arr[i];
      if (n.label == -1) label_neg_counter++;
      if (n.left == -1) left_neg_counter++;
      if (n.right == -1) right_neg_counter++;
      if (n.parent == -1) parent_neg_counter++;
      if (n.time == 0) time_zero_counter++;
    }
    if (left_neg_counter != arg.nsam) {
      num_failure++;
    }
    if (right_neg_counter != arg.nsam) {
      num_failure++;
    }
    if (label_neg_counter != arg.nsam - 1) {
      num_failure++;
    }
    if (time_zero_counter != arg.nsam) num_failure++;
    if (parent_neg_counter != 1) num_failure++;

    return num_failure;
  }
  void print_as_table() {
    std::vector<Node> copy = this->node_arr;
    // sorting
    std::sort(
        copy.begin(), copy.end(), [](const Node &a, const Node &b) -> bool {
          return a.time < b.time || ((a.time == b.time) && (a.label < b.label));
        });

    // header
    std::cout << "index" << '\t' << "time" << '\t' << "label" << '\t'
              << "parent" << '\t' << "left" << '\t' << "right" << '\t'
              << "len_above" << '\n';

    // body part
    for (auto n : copy) {
      std::cout << n.index << '\t' << n.time << '\t' << n.label << '\t'
                << n.parent << '\t' << n.left << '\t' << n.right << '\t'
                << n.len_above << '\n';
    }
  }

 private:
  // find children using array 1
  // duplicate this function to make the recursion using less stack memory
  void find_leaves_recursive_using_arr1(int node_idx) {
    const auto &node = node_arr[node_idx];

    if (node.left == -1)
      leaves1.push_back(node.label);
    else {
      find_leaves_recursive_using_arr1(node.left);
      find_leaves_recursive_using_arr1(node.right);
    }
  }

  void find_leaves_recursive_using_arr2(int node_idx) {
    const auto &node = node_arr[node_idx];
    if (node.left == -1)
      leaves2.push_back(node.label);
    else {
      find_leaves_recursive_using_arr2(node.left);
      find_leaves_recursive_using_arr2(node.right);
    }
  }

  void calculate_tmrca_mat(Args &args) {
    for (ssize_t i = 0; i < num_node; i++) {
      const auto &node = node_arr[i];

      // ignore leave node and focusing on ancestral nodes
      if (node.left == -1) continue;

      // leaves of the left child stored in leaves1
      // find_leaves_using_arr1(node.left);
      leaves1.resize(0);
      find_leaves_recursive_using_arr1(node.left);
      // leaves of the right child stored in leaves2
      leaves2.resize(0);
      find_leaves_recursive_using_arr2(node.right);

      // update tmrca matrix
      for (size_t j = 0; j < leaves1.size(); j++) {
        for (size_t k = 0; k < leaves2.size(); k++) {
          int label1 = leaves1[j];
          int label2 = leaves2[k];

          // since label1 and label2 need to be swaped label1 should be not
          // declared in outter loop
          if (label1 < label2) {
            std::swap(label1, label2);
          }
          // NOTE: msms node label start at 1 so when updating the tmrca_mat,
          // we need take 1 off; macs node label starts with 0
          if (args.tree_sample_nodes_label_start_num == 1) {
            label1--;
            label2--;
          }
          tmrca_mat.at(label1, label2) = node.time;
        }
      }
    }
  }

 public:
  // output ibd and update matrix related to ibd segments
  void find_and_update_ibd(Args &args) {
    size_t num_pairs = ibd_start_mat.get_array_size();
    for (size_t pair = 0; pair < num_pairs; pair++) {
      // use reference to be update
      auto &ibd_time = ibd_time_mat.at(pair);
      auto &ibd_end = ibd_end_mat.at(pair);
      auto &ibd_start = ibd_start_mat.at(pair);
      auto cur_tmrca = tmrca_mat.at(pair);
      // ibd continues
      if (fabs(ibd_time - cur_tmrca) < args.delta)
        ibd_end = tree_interval_end_bp;
      // ibd breaks
      else {
        // if the ibd segments is long enough output ibd
        if (ibd_end - ibd_start > args.min_length_bp) {
          IBDSegments seg;
          // find id1 and id2 using matrix index
          ibd_time_mat.get_matrix_index(pair, seg.id1, seg.id2);
          seg.start = ibd_start;
          seg.end = ibd_end;
          seg.tmrca = ibd_time;
          args.ibdseg_arr.push_back(seg);
        }
        // update ibd matrix elements
        ibd_start = ibd_end = tree_interval_end_bp;
        ibd_time = cur_tmrca;
      }
    }
  }

  // when reaching the last tree, there migh be some IBD continuing to the end
  // of the chromsome we need to output the ibd
  void find_hanging_ibd(Args &args) {
    size_t num_pairs = ibd_start_mat.get_array_size();
    for (size_t pair = 0; pair < num_pairs; pair++) {
      // use reference to be update
      auto &ibd_time = ibd_time_mat.at(pair);
      auto &ibd_end = ibd_end_mat.at(pair);
      auto &ibd_start = ibd_start_mat.at(pair);
      // if the ibd segments is long enough output ibd
      if (ibd_end - ibd_start > args.min_length_bp) {
        IBDSegments seg;
        // find id1 and id2 using matrix index
        ibd_time_mat.get_matrix_index(pair, seg.id1, seg.id2);
        seg.start = ibd_start;
        seg.end = ibd_end;
        seg.tmrca = ibd_time;
        args.ibdseg_arr.push_back(seg);
      }
    }
  }
};

class FileWriter {
  std::ofstream &ofs_log;
  std::ofstream &ofs_ibd;
  std::ofstream &ofs_map;
  std::ofstream &ofs_vcf;

 public:
  FileWriter(std::ofstream &log, std::ofstream &ibd, std::ofstream &map,
             std::ofstream &vcf)
      : ofs_log(log), ofs_ibd(ibd), ofs_map(map), ofs_vcf(vcf) {}

  void write_log(const TreeParser &parser, const Args &args) {
    ofs_log << "\ttree: " << parser.get_tree_counter()
            << "\ttree parsed: " << parser.get_sampled_tree_counter()
            << "\tibd no.: " << args.ibdseg_arr.size() << "\tchromsome parsed: "
            << 1.0 * parser.get_tree_interval_end_bp() / args.seqlen << '\n';
  }
  void write_args_to_log(const Args &args) {
    // TODO
    ofs_log << "====================================\n";
    ofs_log << "chromN:               " << args.chrom << '\n';
    ofs_log << "bp_per_cm:            " << args.bp_per_cm << '\n';
    ofs_log << "time_scale_factor:    " << args.time_scale_factor << '\n';
    ofs_log << "true_ibd_sample_win:  " << args.true_IBD_sampling_window
            << '\n';
    ofs_log << "hom_het:              " << args.hom_het << '\n';
    ofs_log << "simulation:           " << args.simulation_command_str << '\n';
    if (args.nsam != -1)
      ofs_log << "nsam:                 " << args.nsam << '\n';
    if (args.seqlen != -1)
      ofs_log << "seqlen:               " << args.seqlen << '\n';
    if (args.nsegsites != -1)
      ofs_log << "nsegsites:            " << args.seqlen << '\n';
  }

  void write_vcf(const Args &args) {
    // header
    std::string header;
    header.reserve(10240);
    header += "##fileformat=VCFv4.2\n";
    header += "##source=macs_tree_parser by B.G.\n";
    header += "##contig=<ID=" + std::to_string(args.chrom) +
              ",length=" + std::to_string(args.seqlen) + ">\n";
    header +=
        "##INFO=<ID=PR,Number=1,Type=Flag,Description=\"coverted from ms "
        "format data\">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    // add sample names (homozygous diploid versus heterzygous diploid)
    auto n_individuals = args.hom_het ? args.nsam : args.nsam / 2;
    for (long i = 0; i < n_individuals; i++) {
      header += "\ti";
      header += std::to_string(i);
    }
    header += '\n';
    ofs_vcf << header;

    // body part
    std::string line;
    std::string chrom = std::to_string(args.chrom);
    for (ssize_t site = 0; site < args.nsegsites; site++) {
      line.resize(0);
      line += chrom;
      line += "\t";
      line += std::to_string(args.positions[site]);
      line += "\t.\tA\tT\t.\tPASS\tPR\tGT";

      if (args.hom_het == 1) {
        for (ssize_t ind = 0; ind < args.nsam; ind++) {
          char nucleotide = args.haplotypes[args.nsam * site + ind];
          line += '\t';
          line += nucleotide;
          line += '|';
          line += nucleotide;
        }
      } else {
        for (ssize_t ind = 0; ind < args.nsam / 2; ind++) {
          char nucleotide1 = args.haplotypes[args.nsam * site + ind * 2];
          char nucleotide2 = args.haplotypes[args.nsam * site + ind * 2 + 1];
          line += '\t';
          line += nucleotide1;
          line += '|';
          line += nucleotide2;
        }
      }
      ofs_vcf << line << '\n';
    }
  }

  void write_map(Args &args) {
    ofs_map << args.chrom << " . " << 1.0 / args.bp_per_cm << " " << 1 << '\n';
    ofs_map << args.chrom << " . " << args.seqlen / args.bp_per_cm << " "
            << args.seqlen << '\n';
  }

  void write_ibd(Args &args) {
    for (auto seg : args.ibdseg_arr) {
      int sample1, sample2, hap1, hap2;
      sample1 = seg.id1 / 2;
      hap1 = seg.id1 % 2;
      sample2 = seg.id2 / 2;
      hap2 = seg.id1 % 2;
      // ignore homozygous by descent (.hbd)
      if (sample1 == sample2) continue;
      // ensure the order of samples
      if (sample1 > sample2) {
        std::swap(sample1, sample2);
        std::swap(hap1, hap2);
      }
      ofs_ibd << sample1 << '\t';
      ofs_ibd << hap1 << '\t';
      ofs_ibd << sample2 << '\t';
      ofs_ibd << hap2 << '\t';
      ofs_ibd << args.chrom << '\t';
      ofs_ibd << seg.start << '\t';
      ofs_ibd << seg.end << '\t';
      ofs_ibd << (seg.end - seg.start) * 1.0 / args.bp_per_cm << '\t';
      ofs_ibd << (seg.tmrca) << '\n';
    }
  }
};

int main(int argc, char *argv[]) {
  // TreeParser::test();
  auto args = Args(argc, argv);

  // create files
  std::ofstream ofs_ibd(std::to_string(args.chrom) + ".ibd");
  std::ofstream ofs_map(std::to_string(args.chrom) + ".map");
  std::ofstream ofs_log(std::to_string(args.chrom) + ".log");
  std::ofstream ofs_vcf(std::to_string(args.chrom) + ".vcf");

  auto reader = MsLikeFileReader(std::cin);
  auto writer = FileWriter(ofs_log, ofs_ibd, ofs_map, ofs_vcf);

  reader.parse_meta(args);
  writer.write_args_to_log(args);
  // parser rely nsam from meta info from the reader
  auto parser = TreeParser(args.nsam);

  reader.pass_util_tree_lines();
  while (reader.iter_tree_line()) {
    auto &tree_line = reader.get_line_buffer();
    auto res = parser.parse_tree(tree_line, args);
    // skip updating when tree is not sampled
    if (res != 0) continue;
    parser.find_and_update_ibd(args);
    writer.write_log(parser, args);
  }
  parser.find_hanging_ibd(args);
  writer.write_ibd(args);

  reader.parse_segsites(args);
  reader.parse_positions(args);
  reader.parse_haplotypes(args);

  writer.write_vcf(args);
  writer.write_map(args);

  return 0;
}
