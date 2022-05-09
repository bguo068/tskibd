import tskit
import numpy as np
import pandas as pd
import numba
import pyranges as pr


def find_true_ibd_from_ts(
    ts: tskit.TreeSequence,
    sample_window: int = 10_000,
    min_bp: int = 2 * 1000_000,
):
    ibd_cols = ["Id1", "Id2", "Start", "End", "Ancestor"]
    nsam = ts.num_samples
    seqlen = ts.sequence_length

    # initialize matrices to store running information of ibd segment and current tree
    # 1. ancestral node for current tree
    A = np.zeros(shape=(nsam, nsam), dtype=np.uint32)
    # 2. start pos for running IBD segment
    Bseg = np.zeros(shape=(nsam, nsam), dtype=np.uint32)
    # 3. end pos for running IBD segment
    Eseg = np.zeros(shape=(nsam, nsam), dtype=np.uint32)
    # 4. ancestral node for running IBD segments
    Aseg = np.zeros(shape=(nsam, nsam), dtype=np.uint32)
    # 5. store array for IBD segments, the num of rows will be expanded as needed
    ResIBD = np.zeros(shape=(nsam, len(ibd_cols)), dtype=np.uint32)  # for ibd segments
    # 6. indicator matrix for whether elements are in the lower triangle of the matrix
    is_lower_triangle = np.tril(np.ones((nsam, nsam), dtype=np.uint8), -1)

    # counters
    tree_counter = 0  # number of all trees that has been iterated
    sampled_tree_counter = 0  # number of trees that has been sampled
    ibdseg_counter = 0  # number IBD segments has been stored
    bifurcation_node_counter = 0
    multifucation_node_counter = 0  # polytomy

    # numba function to store ibd to the result matrix
    @numba.njit
    def store_ibd(Bseg, Eseg, Aseg, is_output, ibdseg_counter, ResIBD):
        # expand rows if nedded
        i = ibdseg_counter
        for r in range(Bseg.shape[0]):
            for c in range(Bseg.shape[1]):
                if is_output[r, c] > 0:
                    i += 1
                    ResIBD[i][0] = r
                    ResIBD[i][1] = c
                    ResIBD[i][2] = Bseg[r, c]
                    ResIBD[i][3] = Eseg[r, c]
                    ResIBD[i][4] = Aseg[r, c]

    # iterate tree in tree sequences
    for tree in ts.trees():

        # make sure the tree sequence is complete with coalescence
        assert tree.num_roots == 1

        tree_counter += 1
        tree_end = int(tree.interval.right)
        tree_start = int(tree.interval.left)

        # only sample trees whose interval contains >=1 multiple of sample_window
        multiple = tree_end // sample_window * sample_window
        # if not sampled, skip the following part of the loop
        if multiple == tree_end or multiple < tree_start:
            continue

        # if the tree is to be sampled
        sampled_tree_counter += 1

        # Calculate common ancestor matrix by looking at leaves of each children
        # of the ancetral nodes. Any pair of nodes with one being a leave node
        # of the right child and the other a leave of the left child share the
        # ancestral node as the most common ancestor
        internal_nodes = [n for n in tree.nodes() if tree.is_internal(n)]
        # get most recent common ancester matrix
        for n in internal_nodes:

            # find all children of the ancestor instead of just left and right
            # children. This code is try to deal to plotomy situation
            left_child = tree.left_child(n)
            right_child = tree.right_child(n)
            children = [left_child]
            while children[-1] != right_child:
                children.append(tree.right_sib(children[-1]))

            if len(children) > 2:
                multifucation_node_counter += 1
            else:
                bifurcation_node_counter += 1

            for i in range(1, len(children)):
                for j in range(i):
                    child1 = children[i]
                    child2 = children[j]
                    child1_leaves = np.fromiter(tree.leaves(child1), dtype=np.uint32)
                    child2_leaves = np.fromiter(tree.leaves(child2), dtype=np.uint32)
                    A[np.ix_(child1_leaves, child2_leaves)] = n
                    A[np.ix_(child2_leaves, child1_leaves)] = n

        # if at the two nearby sampled trees, each sample pair have different common
        # ancesters the the IBD segment is not continuing and thus there is a break
        is_break = A != Aseg

        # the IBD (due to the break) will be stored only when its length is long enough
        is_long_enough = (Eseg - Bseg) >= min_bp

        # to ingore repeated record due to the order of samples in a pair, only store ibd
        # in the lower triangule matrix
        is_output = (is_break & is_long_enough & is_lower_triangle).astype(np.uint8)

        # enlarge the result matrix when needed
        num_new_ibd = np.count_nonzero(is_output)
        if num_new_ibd > 0:
            needed = np.count_nonzero(is_output)
            if ibdseg_counter + needed + 1000 > ResIBD.shape[0]:
                ResIBD = np.vstack((ResIBD, ResIBD))

            store_ibd(Bseg, Eseg, Aseg, is_output, ibdseg_counter, ResIBD)
            ibdseg_counter += num_new_ibd

        # log info
        if sampled_tree_counter % 20 == 1:
            print(
                f"sampled: {sampled_tree_counter}\t"
                f"trees: {tree_counter}\t"
                f"Process{tree.interval.right /seqlen *100: 0.1f}%\t"
                f"Ibd: {ibdseg_counter}\t"
                f"Polytomies: {multifucation_node_counter/(multifucation_node_counter+bifurcation_node_counter)*100:0.1f}%\t"
            )

        # finishing updating the matrices
        Bseg[is_break] = tree_end
        Eseg[:] = tree_end
        Aseg[is_break] = A[is_break]

    # deal with the hanging IBD after the last tree.
    is_long_enough = (Eseg - Bseg) >= min_bp
    is_output = (is_long_enough & is_lower_triangle).astype(np.uint8)
    num_new_ibd = np.count_nonzero(is_output)
    store_ibd(Bseg, Eseg, Aseg, is_output, ibdseg_counter, ResIBD)
    ibdseg_counter += num_new_ibd

    # covert result matrix to dataframe
    df = pd.DataFrame(ResIBD[:ibdseg_counter, :])
    df.columns = ibd_cols

    # get TMCRA by look at node time of common ancestor
    df["Tmrca"] = ts.tables.nodes.time[df["Ancestor"].to_numpy()]

    return df


def make_diploid_filterred_merged_ibd(
    df: pd.DataFrame,
    bp_per_cm: int,
    chrno: int,
    diploid_mode: str = "heterozygote",
    remove_hbd: bool = True,
    remove_ibd_with_mutation=False,
    min_tmrca: float = 1.5,
    return_unmerge_ibd=True,
):

    raw_ibd_cols = ["Id1", "Id2", "Start", "End", "Ancestor", "Tmrca", "HasMutation"]
    assert list(df.columns) == raw_ibd_cols
    assert diploid_mode in ["homozygote", "heterozygote"]

    # make heterozygous dipoid
    # calculate indiviaul/hap
    df["Sample1"] = df.Id1.round().astype(int) // 2
    df["Sample2"] = df.Id2.round().astype(int) // 2
    df["Hap1"] = df.Id1.round().astype(int) % 2
    df["Hap2"] = df.Id2.round().astype(int) % 2

    # fix dtype
    df["Start"] = df.Start.round().astype(int)
    df["End"] = df.End.round().astype(int)

    # calcuate cm
    df["Cm"] = (df.End - df.Start) / bp_per_cm

    # add chromosme info
    df["Chrom"] = chrno

    # select columns
    df = df[
        [
            "Sample1",
            "Hap1",
            "Sample2",
            "Hap2",
            "Chrom",
            "Start",
            "End",
            "Cm",
            "Ancestor",
            "Tmrca",
            "HasMutation",
        ]
    ].copy()

    # reorder sample/hap pairs to ensuare the consistant order
    sel = df.Sample1 > df.Sample2
    tmp = df.loc[sel, ["Sample2", "Hap2", "Sample1", "Hap1"]].to_numpy(int)
    df.loc[sel, ["Sample1", "Hap1", "Sample2", "Hap2"]] = tmp

    # remove ibd with mutation
    if remove_ibd_with_mutation:
        df = df[df.HasMutation == 0]

    # remove hbd
    if remove_hbd:
        df = df[df.Sample1 != df.Sample2]

    # remove closed related IBD segmetns
    df = df[df.Tmrca >= min_tmrca]

    df_unmerged = df.copy()

    # merge ibd by sample1:sample2 pairs
    fake_chr = (
        df["Sample1"].astype(str).str.cat(df[["Sample2", "Chrom"]].astype(str), sep=":")
    )
    gr = pr.PyRanges(chromosomes=fake_chr, starts=df.Start, ends=df.End - 1)
    gr = gr.merge(count=True)
    # then convert back to dataframe
    df = gr.df.Chromosome.str.split(":", expand=True).copy()
    df.columns = ["Sample1", "Sample2", "Chrom"]
    df["Start"] = gr.df.Start
    df["End"] = gr.df.End
    df["Hap1"] = 0
    df["Hap2"] = 0
    df["Cm"] = (gr.df.End - gr.df.Start) / bp_per_cm
    df = df[["Sample1", "Hap1", "Sample2", "Hap2", "Chrom", "Start", "End", "Cm"]]

    return (df, df_unmerged) if return_unmerge_ibd else df


def make_plink_map(chrno: int, bp_per_cm: int, seqlen: int):
    # make a map
    df_map = pd.DataFrame(
        {
            "Chrom": [chrno, chrno],
            "SnpId": ".",
            "Cm": [1.0 / bp_per_cm, seqlen * 1.0 / bp_per_cm],
            "Bp": [1, seqlen],
        }
    )
    return df_map
