
use crate::ComplexRef;
use crate::Reaction;

/// Trait for implementing a single rewrite rule.
pub trait RewriteRule: Send + Sync {

    /// Rule identifier string.
    fn id(&self) -> &'static str;

    /// Rule pattern desciption.
    fn pattern(&self) -> &'static str;

    /// Optional category tag.
    fn category(&self) -> &'static str {
        "default"
    }

    /// Numerical cost of applying this rule
    fn cost(&self) -> f64 {
        1.0
    }

    /// Enumerate and apply all matches
    fn apply<'a>(&'a self, complex: &'a ComplexRef) 
        -> Box<dyn Iterator<Item = (ComplexRef, String)> + 'a>;

    fn get_reactions<'a>(&'a self, reactant: &'a ComplexRef) 
        -> Box<dyn Iterator<Item = Reaction> + 'a> {
        Box::new(
            self.apply(reactant)
            .map(move |(product, rewrite)| {
                //self.on_apply(reactant, &rewrite, &product);
                Reaction::new(
                    reactant.clone(),
                    product,
                    self.id().to_string(),
                    self.cost(),
                    rewrite,
                )
            })
        )
    }

    // Hook for logging/debugging after rule was applied (optional)
    //fn on_apply(&self, _c: &ComplexRef, rewrite: &str, _p: &ComplexRef) {
    //    tracing::debug!(target: "rule", "{} applied: {:?}", self.id(), rewrite);
    //}
}


