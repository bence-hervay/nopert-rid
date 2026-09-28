//! Running and recording a search.
pub mod certificate;
pub mod queue;
pub mod config;
pub mod command;

/// Error type of the caller-supplied parts: evaluators, verifiers, outputs
/// and components.
pub type BoxError = Box<dyn std::error::Error + Send + Sync>;

/// The text of a panic payload, for reporting a caught panic as an error.
pub(crate) fn panic_message(panic: &(dyn std::any::Any + Send)) -> String {
    panic
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| panic.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "non-text panic payload".into())
}
