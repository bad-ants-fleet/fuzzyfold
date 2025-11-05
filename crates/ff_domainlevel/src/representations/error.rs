
/// Error type for domain-level Nussinov dynamic programming operations.
#[derive(Debug, Clone)]
pub enum RegistryError {
    /// A domain name in the sequence was not found in the registry.
    UnknownDomain(String),

    /// Other general-purpose error placeholder (for future extensions).
    Other(String),
}


