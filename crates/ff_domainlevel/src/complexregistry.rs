
use std::fmt;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, Weak};
use std::sync::atomic::{AtomicUsize, Ordering};
use once_cell::sync::Lazy;

use crate::complex::Complex;
use ff_structure::StructureError;

#[derive(Debug)]
pub enum ComplexRegistryError {
    NameAlreadyInUse(String),
    AlreadyNamed(String),
    Internal(String),
    KernelParseFailed(StructureError),
}

impl From<StructureError> for ComplexRegistryError {
    fn from(e: StructureError) -> Self {
        ComplexRegistryError::KernelParseFailed(e)
    }
}

impl fmt::Display for ComplexRegistryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ComplexRegistryError::NameAlreadyInUse(name) =>
                write!(f, "Name '{}' is already in use", name),
            ComplexRegistryError::AlreadyNamed(name) =>
                write!(f, "Complex is already named '{}', cannot rename", name),
            ComplexRegistryError::Internal(msg) =>
                write!(f, "Internal error: {}", msg),
            ComplexRegistryError::KernelParseFailed(e) =>
                write!(f, "Kernel parsing failed: {}", e),
        }
    }
}

impl std::error::Error for ComplexRegistryError {}

pub type ComplexRef = Arc<Complex>;

pub struct ComplexRegistry {
    by_kernel: HashMap<String, Weak<Complex>>,
    by_name: HashMap<String, String>, // name -> kernel
    name_counter: AtomicUsize,
}

pub static REGISTRY: Lazy<Mutex<ComplexRegistry>> = Lazy::new(|| {
    Mutex::new(ComplexRegistry {
        by_kernel: HashMap::new(),
        by_name: HashMap::new(),
        name_counter: AtomicUsize::new(1),
    })
});

impl ComplexRegistry {
    pub fn get_or_create(
        raw_kernel: &str,
        name: Option<&str>,
    ) -> Result<ComplexRef, ComplexRegistryError> {
        let kernel = Complex::cannonize_kernel(raw_kernel);

        let mut reg = REGISTRY.lock().unwrap();
        match reg.by_kernel.get(&kernel).and_then(|w| w.upgrade()) {
            Some(existing) => {
                // check naming conflict if necessary.
                if let Some(new_name) = name {
                    match existing.name() {
                        Some(current_name) if current_name != new_name => {
                            return Err(ComplexRegistryError::AlreadyNamed(current_name));
                        }
                        None => {
                            if reg.by_name.get(new_name).map(|v| v != &kernel).unwrap_or(false) {
                                return Err(ComplexRegistryError::NameAlreadyInUse(new_name.to_string()));
                            }
                            existing.set_name(new_name)
                                .map_err(|e| ComplexRegistryError::Internal(e))?;
                            reg.by_name.insert(new_name.to_string(), kernel.clone());
                        }
                        _ => {}
                    }
                }
                return Ok(existing)
            }
            None => {
                reg.by_kernel.remove(&kernel);
            }
        }

        // Name assignment (manual or automatic)
        if let Some(n) = name {
            if reg.by_name.get(n).map(|v| v != &kernel).unwrap_or(false) {
                return Err(ComplexRegistryError::NameAlreadyInUse(n.to_string()));
            }
            reg.by_name.insert(n.to_string(), kernel.clone());
        } 

        // Otherwise: parse and construct a new Complex
        let complex = Complex::from_kernel(&kernel, name)?;
        let arc = Arc::new(complex);
        reg.by_kernel.insert(kernel.clone(), Arc::downgrade(&arc));
        Ok(arc)
    }

    pub fn assign_auto_name(complex: &Arc<Complex>, prefix: &str) -> Result<String, ComplexRegistryError> {
        let mut reg = REGISTRY.lock().unwrap();

        // Already named? That's an error.
        if let Some(existing) = complex.name() {
            return Err(ComplexRegistryError::AlreadyNamed(existing));
        }

        // Generate and reserve a unique name
        loop {
            let candidate = format!("{}{}", prefix, reg.name_counter.fetch_add(1, Ordering::Relaxed));
            if reg.by_name.contains_key(&candidate) {
                continue;
            }
            complex.set_name(&candidate)
                .map_err(|e| ComplexRegistryError::Internal(e))?;
            reg.by_name.insert(candidate.clone(), complex.kernel().to_string());
            return Ok(candidate);
        }
    }

    #[cfg(test)]
    pub fn clear() {
        let mut reg = REGISTRY.lock().unwrap();
        reg.by_kernel.clear();
        reg.by_name.clear();
        reg.name_counter.store(1, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_and_retrieve_named_complex() {
        let a = ComplexRegistry::get_or_create("tcarn a b a*", Some("foo")).unwrap();
        let b = ComplexRegistry::get_or_create("tcarn a b a*", Some("foo")).unwrap();

        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(a.name(), Some("foo".to_owned()));
    }

    #[test]
    fn test_duplicate_name_error() {
        let _c1 = ComplexRegistry::get_or_create("a x a*", Some("x")).unwrap();
        let c2 = ComplexRegistry::get_or_create("a( x )", Some("x"));

        assert!(matches!(
            c2,
            Err(ComplexRegistryError::NameAlreadyInUse(name)) if name == "x"
        ));
    }

    #[test]
    fn test_already_named_error() {
        // this tests requires that a variable is assigned in the first line
        let _c1 = ComplexRegistry::get_or_create("a y a*", Some("y")).unwrap();
        let c2 = ComplexRegistry::get_or_create("a y a*", Some("other"));

        assert!(matches!(
            c2,
            Err(ComplexRegistryError::AlreadyNamed(name)) if name == "y"
        ));
    }

    #[test]
    fn test_name_assignment() {
        let c = ComplexRegistry::get_or_create("tn x y z", None).unwrap();
        let name = ComplexRegistry::assign_auto_name(&c, "tn_").unwrap();
        assert!(name.starts_with("tn_"));
    }

    #[test]
    fn test_canonicalization_equivalence() {
        let a = ComplexRegistry::get_or_create(" a   b  c ", Some("c1")).unwrap();
        let b = ComplexRegistry::get_or_create("a b c", Some("c1")).unwrap();

        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(a.kernel(), "a b c");
    }

    #[test]
    fn test_parse_failure_handling() {
        let result = ComplexRegistry::get_or_create("a( b", None);
        assert!(matches!(
            result,
            Err(ComplexRegistryError::KernelParseFailed(StructureError::InvalidToken(..)))
        ));
    }
}

