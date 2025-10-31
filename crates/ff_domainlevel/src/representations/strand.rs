
//!
//!NOTE: Thinking about StrandRef like implemented for Domain vs DomainRef
//!

use ahash::AHashMap;

use crate::DomainRefVec;

#[allow(dead_code)]
pub struct Strand {
    strand: DomainRefVec,
    name: String,
}

pub struct StrandRegistry {
    strands: AHashMap<String, DomainRefVec>,
    names: AHashMap<DomainRefVec, String>,
    counter: usize,
}

#[allow(dead_code)]
impl StrandRegistry {
    pub fn new() -> Self {
        Self {
            strands: AHashMap::default(),        
            names: AHashMap::default(),
            counter: 0,
        }
    }

    pub fn intern(&mut self, strand: &DomainRefVec, name: Option<String>) -> DomainRefVec {
        if let Some(existing_name) = self.names.get(strand) {
            if let Some(ref name) = name {
                assert_eq!(
                    existing_name, name,
                    "Tried to assign name '{name}' to existing strand '{existing_name}'"
                );
            }
            return strand.clone(); // already interned
        }

        let assigned_name = match name {
            Some(n) => {
                assert!(
                    !self.strands.contains_key(&n),
                    "Strand name '{}' already exists",
                    n
                );
                n
            }
            None => {
                let mut n;
                loop {
                    n = format!("s{}", self.counter);
                    self.counter += 1;
                    if !self.strands.contains_key(&n) {
                        break n;
                    }
                }
            }
        };

        self.strands.insert(assigned_name.clone(), strand.clone());
        self.names.insert(strand.clone(), assigned_name);

        strand.clone()
    }

    pub fn get_by_name(&self, name: &str) -> Option<DomainRefVec> {
        self.strands.get(name).cloned()
    }

    pub fn get_name(&self, strand: &DomainRefVec) -> Option<&String> {
        self.names.get(strand)
    }
}



