// Id1	Id2	Start	End	Ancestor	Tmrca	HasMutation
// 153	33	0	30093	4308	17	0
// 153	117	0	30093	4308	17	0
// 153	119	0	30093	3149	9	0

use csv;
use serde::{Deserialize, Serialize};
mod intervaltree;
use clap::Parser;
use std::path::PathBuf;

use itertools::{
    EitherOrBoth::{Both, Left, Right},
    Itertools,
};
use slice_group_by::{GroupBy, GroupByMut};

use intervaltree::IntervalTree;

#[derive(Debug, PartialEq, PartialOrd, Deserialize, Serialize, Clone, Copy)]
struct IbdSeg {
    #[serde(rename = "Id1")]
    id1: u32,
    #[serde(rename = "Id2")]
    id2: u32,
    #[serde(rename = "Start")]
    start: u32,
    #[serde(rename = "End")]
    end: u32,
    #[serde(rename = "Ancestor")]
    anc: u32,
    #[serde(rename = "Tmrca")]
    tmrca: f32,
    #[serde(rename = "HasMutation")]
    hasmut: u8,
}

#[derive(PartialEq, Debug)]
struct IbdSet(Vec<IbdSeg>);

impl IbdSet {
    fn from_tskibd_file(p: &str) -> Self {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(p)
            .expect(&format!("Can't read file: {p}"));
        //let header = rdr.headers().unwrap().clone();
        let mut record = csv::ByteRecord::new();

        let mut v = vec![];
        while rdr.read_byte_record(&mut record).unwrap() {
            let ibdseg: IbdSeg = record.deserialize(None).unwrap();
            v.push(ibdseg);
        }
        Self(v)
    }
    fn to_tskibd_file(&self, p: &str) {
        let mut wrt = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(p)
            .expect(&format!("Cannot write to file {}", p));
        for ibdseg in self.0.iter() {
            wrt.serialize(ibdseg).unwrap();
        }
    }
    fn normalize(&mut self) {
        for s in self.0.iter_mut() {
            if s.id2 > s.id1 {
                std::mem::swap(&mut s.id1, &mut s.id2);
            }
        }
    }
    fn sort(&mut self) {
        self.0.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }
    /// Only keep IBD segment with tmrca < 1.5
    fn filt_by_tmrca(&mut self) {
        self.0.retain(|s| s.tmrca < 1.5);
    }
    fn filt_by_overlapping(&mut self, other: &Self) {
        // ensure normalized and others is filt by tmrca
        assert!(self.0.iter().all(|a| a.id1 > a.id2));
        assert!(other.0.iter().all(|a| (a.id1 > a.id2) && (a.tmrca < 1.5)));

        // ensure sorted
        assert!(self
            .0
            .iter()
            .zip(self.0.iter().skip(1))
            .all(|(a, b)| a.le(b)));
        assert!(other
            .0
            .iter()
            .zip(other.0.iter().skip(1))
            .all(|(a, b)| a.le(b)));

        let it1 = self
            .0
            .linear_group_by_key_mut(|s| (s.id1, s.id2, s.start, s.end));
        let it2 = other.0[..].linear_group_by_key(|s| (s.id1, s.id2, s.start, s.end));

        let mut tree = IntervalTree::new(100);
        it1.merge_join_by(it2, |a, b| (a[0].id1, a[0].id2).cmp(&(b[0].id1, b[0].id2)))
            .for_each(|res| match res {
                Left(_) => {}
                Right(_) => {}
                Both(a, b) => {
                    tree.clear_and_fill_with_iter(b.iter().map(|x| (x.start..x.end, ())));
                    for seg in a {
                        // if overlap, mark this segs with (id1=0, id2=0) which will be clear out
                        // later
                        if let Some(_) = tree.query(seg.start..seg.end).next() {
                            seg.id1 = 0;
                            seg.id2 = 0;
                        }
                    }
                }
            });

        self.0.retain(|seg| (seg.id1 != 0) && (seg.id2 != 0));
    }
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// path to infer IBD file (in tskibd format)
    infer_set: PathBuf,

    /// path to infer IBD file (in tskibd format)
    true_set: PathBuf,

    /// path to output IBD file
    #[arg(short = 'o', long, required = true)]
    out: PathBuf,
}

fn main() {
    let cli = Cli::parse();

    let mut true_set = IbdSet::from_tskibd_file(cli.true_set.to_str().unwrap());
    let mut infer_set = IbdSet::from_tskibd_file(cli.infer_set.to_str().unwrap());
    true_set.normalize();
    true_set.sort();
    true_set.filt_by_tmrca();
    infer_set.normalize();
    infer_set.sort();
    infer_set.filt_by_overlapping(&true_set);
    infer_set.to_tskibd_file(cli.out.to_str().unwrap());
}

#[cfg(test)]
mod tests {
    use super::*;
    fn make_ibdset(v: Vec<(u32, u32, u32, u32, u32, f32, u8)>) -> IbdSet {
        let mut set = vec![];
        for (id1, id2, start, end, anc, tmrca, hasmut) in v {
            let seg = IbdSeg {
                id1,
                id2,
                start,
                end,
                anc,
                tmrca,
                hasmut,
            };
            set.push(seg);
        }
        IbdSet(set)
    }

    #[test]
    fn test1() {
        let mut true_set = make_ibdset(vec![
            (1, 0, 10, 100, 20, 1.0, 0),
            (1, 0, 1000, 1100, 20, 1.0, 0),
            (2, 1, 1000, 1100, 20, 20.0, 0),
            (99, 10, 20, 200, 20, 1.0, 0),
            (99, 10, 2000, 2200, 20, 30.0, 0),
            (99, 11, 2000, 2200, 20, 20.0, 0),
        ]);
        let mut infer_set = make_ibdset(vec![
            (1, 0, 90, 200, 20, 100.0, 0),
            (2, 1, 1000, 1100, 20, 100.0, 0),
            (99, 10, 20, 200, 20, 100.0, 0),
            (99, 10, 2000, 2200, 20, 100.0, 0),
        ]);
        let exp_out = make_ibdset(vec![
            (2, 1, 1000, 1100, 20, 100.0, 0),
            (99, 10, 2000, 2200, 20, 100.0, 0),
        ]);
        true_set.normalize();
        true_set.sort();
        true_set.filt_by_tmrca();
        infer_set.normalize();
        infer_set.sort();
        infer_set.filt_by_overlapping(&true_set);

        assert_eq!(infer_set, exp_out);
    }
}
