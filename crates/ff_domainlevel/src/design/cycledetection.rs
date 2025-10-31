
use nohash_hasher::IntMap;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            return false; // Cycle detected!
        }

        // Union by rank
        if self.rank[root_x] < self.rank[root_y] {
            self.parent[root_x] = root_y;
        } else if self.rank[root_x] > self.rank[root_y] {
            self.parent[root_y] = root_x;
        } else {
            self.parent[root_y] = root_x;
            self.rank[root_x] += 1;
        }
        true
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]); // Path compression
        }
        self.parent[x]
    }

    pub fn connected_components(&mut self) -> Vec<Vec<usize>> {
        let mut groups: IntMap<usize, Vec<usize>> = IntMap::default();
        for i in 0..self.parent.len() {
            let root = self.find(i); // ensures path compression too
            groups.entry(root).or_default().push(i);
        }
        groups.into_values().collect()
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_union_find() {
        let mut uf = UnionFind::new(5);

        // Initially, each node is its own parent
        for i in 0..5 {
            assert_eq!(uf.find(i), i);
        }

        // Union two elements
        assert!(uf.union(0, 1));
        assert_eq!(uf.find(0), uf.find(1));

        // Union another pair
        assert!(uf.union(2, 3));
        assert_eq!(uf.find(2), uf.find(3));

        // Connect the two sets
        assert!(uf.union(1, 2));
        assert_eq!(uf.find(0), uf.find(3));
    }

    #[test]
    fn test_cycle_detection() {
        let mut uf = UnionFind::new(4);

        assert!(uf.union(0, 1));
        assert!(uf.union(1, 2));
        assert!(uf.union(2, 3));

        // 0-1-2-3 are all connected now
        assert!(!uf.union(3, 0)); // Would form a cycle
    }

    #[test]
    fn test_union_idempotency() {
        let mut uf = UnionFind::new(2);

        assert!(uf.union(0, 1));
        // Repeated union on same components returns false (already connected)
        assert!(!uf.union(0, 1));
        assert!(!uf.union(1, 0));
        assert_eq!(uf.find(0), uf.find(1));
    }

    #[test]
    fn test_rank_greater_than_one() {
        let mut uf = UnionFind::new(8);
    
        // Step 1: union(0,1) and union(2,3) → both roots will have rank 1
        assert!(uf.union(0, 1));
        assert!(uf.union(2, 3));
    
        let r0 = uf.find(0);
        assert_eq!(uf.rank[r0], 1);
        let r2 = uf.find(2);
        assert_eq!(uf.rank[r2], 1);
    
        // Step 2: union(0,2) → ranks equal → merged root rank should become 2
        assert!(uf.union(0, 2));
        let root = uf.find(0);
        assert_eq!(uf.rank[root], 2, "Expected rank to be 2");

        assert!(!uf.union(3, 1));
    }

    #[test]
    fn test_clone_behavior() {
        let mut uf1 = UnionFind::new(3);
        uf1.union(0, 1);

        let mut uf2 = uf1.clone();
        assert_eq!(uf1.find(0), uf2.find(1));

        // Modify the clone
        uf2.union(1, 2);

        // Ensure original is unaffected
        assert_ne!(uf1.find(2), uf2.find(2));
    }

    #[test]
    fn test_connected_components() {
        let mut uf = UnionFind::new(7);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(3, 4);

        let mut components = uf.connected_components();
        println!("{:?}", components);
        for group in &mut components {
            group.sort(); // Ensure consistent order for test
        }
        components.sort(); // Sort groups for comparison

        assert_eq!(components, vec![vec![0, 1, 2], vec![3, 4], vec![5], vec![6]]);
    }

}


