/// domain.rs (the Domain and DomainRegistry)
///
/// The Domain (name, length) is a typical abstract representation for the specification of nucleic
/// acid complementarity. Because domains may occur in multiple sequences, and in many (enumerated)
/// complexes, the DomainRegistry serves as a central overview on all domains in system. It holds a
/// strong reference on each Domain, such that sequences, complexes, etc. can share the same Heap
/// allocation with runtime reference counting. It also holds the logic of which domains are
/// complementary -- which may be extended on demand.
///
/// For now, we do not support thread safety of the DomainRegistry. Even in multithreaded setups, 
/// it may well be faster to initialize separate DomainRegistries in different threads. Generally,
/// the DomainRegistry will be comparatively cheap to initiate, but then does typically not change,
/// it provides a bit faster access to complementarity checks and memory-efficient storage of
/// sequences.
///
/// Note: An attempt to use weak references to complementary domains from within the Domain has
/// turned out very complicated / impossible without RefCell. I think it is not worth it, probably
/// a bad design choice anyway. Let the DomainRegistry handle it.
///

use std::fmt;
use std::rc::Rc;
use rustc_hash::FxHashMap;

#[derive(Debug, Eq, Hash, PartialEq)]
pub struct Domain {
    pub name: String,
    pub length: usize,
}

impl fmt::Display for Domain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //write!(f, "{} [len={}]", self.name, self.length)
        write!(f, "{}", self.name)
    }
}

pub type DomainRef = Rc<Domain>;
pub type DomainRefVec = Vec<DomainRef>;

#[derive(Default)]
pub struct DomainRegistry {
    domains: FxHashMap<String, DomainRef>,
    complements: FxHashMap<DomainRef, DomainRef>, // symmetric
}

impl DomainRegistry {
    pub fn new() -> Self {
        Self {
            domains: FxHashMap::default(),        
            complements: FxHashMap::default(),
        }
    }

    pub fn intern(&mut self, name: &str, length: usize) -> DomainRef {
        // First, check if this domain already exists
        if let Some(existing) = self.domains.get(name) {
            assert_eq!(existing.length, length);
            return existing.clone();
        }

        // Compute canonical complement name
        let (name1, name2) = if name.ends_with('*') {
            (name.trim_end_matches('*').to_owned(), name.to_owned())
        } else {
            (name.to_owned(), format!("{}*", name))
        };

        assert!(!self.domains.contains_key(&name1), "Domain {} already exists", name1);
        assert!(!self.domains.contains_key(&name2), "Domain {} already exists", name2);

        let d1 = Rc::new(Domain {
            name: name1.clone(),
            length,
        });

        let d2 = Rc::new(Domain {
            name: name2.clone(),
            length,
        });

        self.domains.insert(name1, d1.clone());
        self.domains.insert(name2, d2.clone());

        self.complements.insert(d1.clone(), d2.clone());
        self.complements.insert(d2, d1.clone());

        d1
    }

    /// Retrieve the DomainRef if only the name of the Domain is known.
    pub fn get(&self, name: &str) -> Option<DomainRef> {
        self.domains.get(name).cloned()
    }

    /// Get the complementary DomainRef to a given DomainRef
    pub fn get_complement(&self, d: &DomainRef) -> DomainRef {
        self.complements.get(d).unwrap().clone()
    }

    /// Check if two Domains are complementary.
    pub fn are_complements(&self, a: &DomainRef, b: &DomainRef) -> bool {
        match self.complements.get(a) {
            Some(c) => Rc::ptr_eq(c, b),
            None => false,
        }
    }
}


/// Checks whether two sequence symbols are complementary (e.g., "a" <-> "a*")
pub fn is_complement(a: &str, b: &str) -> bool {
    a == format!("{}*", b) || b == format!("{}*", a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_complement() {
        assert!(is_complement("a", "a*"));
        assert!(is_complement("b*", "b"));
        assert!(!is_complement("a", "b*"));
        assert!(!is_complement("x", "x"));
    }

    #[test]
    fn test_symmetric_interning() {
        let mut reg = DomainRegistry::new();
        let a = reg.intern("a", 3);
        let astar = reg.intern("a*", 3);

        println!("{}", a);
        println!("{}", astar);
        println!("{}", reg.get_complement(&a));

        assert!(Rc::ptr_eq(&a, &reg.intern("a", 3)));
        assert!(Rc::ptr_eq(&astar, &reg.intern("a*", 3)));
        assert!(reg.are_complements(&a, &astar));
        assert!(reg.are_complements(&astar, &a));
        assert!(Rc::ptr_eq(&a, &reg.get_complement(&reg.get_complement(&a))));
    }

    #[test]
    #[should_panic]
    fn test_mismatched_length_panics() {
        let mut reg = DomainRegistry::new();
        reg.intern("x", 3);
        reg.intern("x", 4); // should panic
    }

}
