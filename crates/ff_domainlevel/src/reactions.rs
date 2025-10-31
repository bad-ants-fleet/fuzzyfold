use std::fmt;
use crate::complexregistry::ComplexRef;


#[derive(Clone)]
pub struct Reaction {
    reactant: ComplexRef,
    product: ComplexRef,
    rule: String, 
    _rate_const: f64, 
    rewrite: String, 
}

impl Reaction {
    pub fn new(
        reactant: ComplexRef,
        product: ComplexRef,
        rule: String,
        _rate_const: f64,
        rewrite: String,
    ) -> Self {
        Self {
            reactant,
            product,
            rule,
            _rate_const,
            rewrite,
        }
    }

    pub fn reactant(&self) -> &ComplexRef {
        &self.reactant
    }

    pub fn product(&self) -> &ComplexRef {
        &self.product
    }

    pub fn rule(&self) -> &str {
        &self.rule
    }

    pub fn rewrite(&self) -> &str {
        &self.rewrite
    }
}

impl fmt::Display for Reaction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}] {} → {}",
            self.rule,
            match self.reactant.name() {
                Some(name) => name,
                None => self.reactant.kernel().to_string(),
            },
            match self.product.name() {
                Some(name) => name,
                None => self.product.kernel().to_string(),
            },
        )
    }
}

impl fmt::Debug for Reaction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Reaction")
            .field("reactant", &self.reactant.kernel())
            .field("product", &self.product.kernel())
            .field("rule", &self.rule)
            .field("rate_const_", &self._rate_const)
            .field("rewrite", &self.rewrite)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::complexregistry::ComplexRegistry;

    #[test]
    fn test_basic_reaction_display() {
        let reactant = ComplexRegistry::get_or_create("a b a*", Some("R")).unwrap();
        let product = ComplexRegistry::get_or_create("a( b )", Some("P")).unwrap();

        let reaction = Reaction::new(reactant.clone(), product.clone(), "R1.1".into(), 1.0, "1a 3a* -> 1a( 3)".into());
        assert_eq!(format!("{}", reaction), "[R1.1] R → P");
    }

    #[test]
    fn test_reaction_debug_output() {
        let r = ComplexRegistry::get_or_create("x y z", Some("X")).unwrap();
        let p = ComplexRegistry::get_or_create("x( y z )", Some("Y")).unwrap();

        let rxn = Reaction::new(r.clone(), p.clone(), "foo".into(), 0.5, "x y z → x( y z )".into());
        let debug = format!("{:?}", rxn);
        assert!(debug.contains("Reaction"));
        assert!(debug.contains("x y z"));
        assert!(debug.contains("x( y z )"));
        assert!(debug.contains("foo"));
        assert!(debug.contains("0.5"));
    }

    #[test]
    fn test_reaction_fields() {
        let a = ComplexRegistry::get_or_create("a b", Some("A")).unwrap();
        let b = ComplexRegistry::get_or_create("a( b )", Some("B")).unwrap();

        let r = Reaction::new(a.clone(), b.clone(), "fold".into(), 2.0, "a b → a( b )".into());

        assert!(Arc::ptr_eq(r.reactant(), &a));
        assert!(Arc::ptr_eq(r.product(), &b));
        assert_eq!(r.rule(), "fold");
        assert_eq!(r.rewrite(), "a b → a( b )");
    }
}
