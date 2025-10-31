
use std::fmt;
use std::sync::RwLock;
use ff_structure::DotBracket;
use ff_structure::PairTable;
use ff_structure::StructureError;
use ff_structure::DotBracketVec;

//TODO: make a low-level Kernel struct that serves as I/O using 
// * DomainRefVec
// * DotBracketVec

pub struct Complex {
    kernel: String, // this is the canonical form.
    sequence: Vec<String>,
    structure: DotBracketVec, 
    pair_table: PairTable,
    name: RwLock<Option<String>>,
}

impl Complex {
    pub fn new(
        kernel: String,
        sequence: Vec<String>,
        structure: DotBracketVec,
        pair_table: PairTable,
        name: Option<&str>,
    ) -> Self {
        Self {
            kernel,
            sequence,
            structure,
            pair_table,
            name: RwLock::new(name.map(|s| s.to_string())),
        }
    }

    pub fn from_kernel(raw_kernel: &str, name: Option<&str>) -> Result<Self, StructureError> {
        let kernel = Self::cannonize_kernel(raw_kernel);
        let (sequence, structure) = parse_kernel(&kernel)?;
        let pair_table = PairTable::try_from(&structure).expect("This must be a valid kernel structure.");
        assert_eq!(kernel, get_kernel(&sequence, &structure), "Canonical kernel mismatch");

        Ok(Self::new(kernel, sequence, structure, pair_table, name))
    }
    
    pub fn set_name(&self, name: &str) -> Result<(), String> {
        let mut lock = self.name.write().map_err(|_| "RwLock poisoned")?;
        if lock.is_some() {
            return Err("Complex is already named".into());
        }
        *lock = Some(name.to_string());
        Ok(())
    }

    pub fn name(&self) -> Option<String> {
        self.name.read().ok().and_then(|r| r.clone())
    }

    pub fn kernel(&self) -> &str {
        &self.kernel
    }

    pub fn sequence(&self) -> &[String] {
        &self.sequence
    }

    pub fn structure(&self) -> &DotBracketVec {
        &self.structure
    }

    pub fn pair_table(&self) -> &PairTable {
        &self.pair_table
    }

    pub fn cannonize_kernel(raw: &str) -> String {
        raw.split_whitespace().collect::<Vec<_>>().join(" ")
    }
}

impl fmt::Debug for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = self.name.read().unwrap();
        f.debug_struct("Complex")
            .field("kernel", &self.kernel)
            .field("sequence", &self.sequence)
            .field("structure", &self.structure)
            .field("name", &*name)
            .finish()
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.name.read().unwrap().as_ref() {
            Some(name) => write!(f, "{}: {}", name, self.kernel),
            None => write!(f, "{}", self.kernel),
        }
    }
}

fn parse_kernel(kernel: &str) -> Result<(Vec<String>, DotBracketVec), StructureError> {
    let tokens: Vec<&str> = kernel.split_whitespace().collect();
    let mut sequence: Vec<String> = Vec::new();
    let mut structure: Vec<DotBracket> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    
    for (idx, &token) in tokens.iter().enumerate() {
        if token == ")" {
            match stack.pop() {
                Some(open_idx) => {
                    let opener = &sequence[open_idx];
                    let inferred = if opener.ends_with('*') {
                        opener.trim_end_matches('*').to_string()
                    } else {
                        format!("{}*", opener)
                    };
                    sequence.push(inferred);
                    structure.push(DotBracket::Close);
                }
                None => return Err(StructureError::InvalidToken(token.to_string(), "kernel".to_string(), idx )),
            }
        } else if token.ends_with('(') {
            let seq = token.strip_suffix('(').unwrap();
            if seq.is_empty() {
                return Err(StructureError::InvalidToken(token.to_string(),  "kernel".to_string(), idx));
            }
            sequence.push(seq.to_string());
            structure.push(DotBracket::Open);
            stack.push(sequence.len() - 1);
        } else if token.ends_with(')') || token.contains('(') || token.contains(')') {
            return Err(StructureError::InvalidToken(token.to_string(),  "kernel".to_string(), idx));
        } else {
            sequence.push(token.to_string());
            structure.push(DotBracket::Unpaired)
        }
    }
    if !stack.is_empty() {
        return Err(StructureError::InvalidToken("end of string".to_string(),  "kernel".to_string(), tokens.len()));
    }
    Ok((sequence, DotBracketVec(structure)))
}

pub fn get_kernel(sequence: &[String], structure: &[DotBracket]) -> String {
    assert_eq!(sequence.len(), structure.len(), "Mismatched sequence and structure lengths");
    let mut tokens = Vec::new();
    for (seq, struc) in sequence.iter().zip(structure.iter()) {
        let token = match struc {
            DotBracket::Unpaired => seq.clone(), // just the sequence
            DotBracket::Open => format!("{}(", seq), // explicit opener
            DotBracket::Close => ")".to_string(), // implicit sequence, inferred
            DotBracket::Break => unreachable!("not implemented!"), // implicit sequence, inferred
        };
        tokens.push(token);
    }
    tokens.join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_kernel_unpaired() {
        let c = Complex::from_kernel("a b c", None).unwrap();
        assert_eq!(c.kernel(), "a b c");
        assert_eq!(c.name(), None);
        assert_eq!(c.sequence(), &vec!["a", "b", "c"]);
        assert_eq!(
            c.structure(),
            &DotBracketVec(vec![DotBracket::Unpaired, DotBracket::Unpaired, DotBracket::Unpaired])
        );
    }

    #[test]
    fn test_from_kernel_with_structure() {
        let c = Complex::from_kernel("a( b )", None).unwrap();
        assert_eq!(c.kernel(), "a( b )");
        assert_eq!(c.name(), None);
        assert_eq!(c.sequence(), &vec!["a", "b", "a*"]);
        assert_eq!(
            c.structure(),
            &DotBracketVec(vec![DotBracket::Open, DotBracket::Unpaired, DotBracket::Close])
        );
    }

    #[test]
    fn test_complex_set_name() {
        let c = Complex::from_kernel("a b b", None).unwrap();
        assert_eq!(format!("{}", c), "a b b");
        let _ = c.set_name("testname");
        assert_eq!(c.name(), Some("testname".to_owned()));
        assert_eq!(format!("{}", c), "testname: a b b");
    }

    #[test]
    fn test_cannonize_kernel() {
        let messy = "   a   b    c ";
        let clean = Complex::cannonize_kernel(messy);
        assert_eq!(clean, "a b c");
    }

    #[test]
    fn test_parse_kernel_unpaired() {
        let (seq, struc) = parse_kernel("a b c").unwrap();
        assert_eq!(seq, vec!["a", "b", "c"]);
        assert_eq!(struc.0, vec![
            DotBracket::Unpaired,
            DotBracket::Unpaired,
            DotBracket::Unpaired,
        ]);
    }

    #[test]
    fn test_parse_kernel_with_pairing() {
        let (seq, struc) = parse_kernel("x( y )").unwrap();
        assert_eq!(seq, vec!["x", "y", "x*"]);
        assert_eq!(struc.0, vec![DotBracket::Open, DotBracket::Unpaired, DotBracket::Close]);
    }

    #[test]
    fn test_parse_kernel_nested() {
        let (seq, struc) = parse_kernel("a( b( c ) d )").unwrap();
        assert_eq!(seq, vec!["a", "b", "c", "b*", "d", "a*"]);
        assert_eq!(struc.0, vec![
            DotBracket::Open,
            DotBracket::Open,
            DotBracket::Unpaired,
            DotBracket::Close,
            DotBracket::Unpaired,
            DotBracket::Close,
        ]);
    }

    #[test]
    fn test_parse_kernel_errors() {
        assert!(matches!(parse_kernel("a("), Err(StructureError::InvalidToken(..))));
        assert!(matches!(parse_kernel("a)"), Err(StructureError::InvalidToken(..))));
        assert!(matches!(parse_kernel("x((y)"), Err(StructureError::InvalidToken(..))));
    }

}
