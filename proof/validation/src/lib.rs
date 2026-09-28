//! Experiments on the `rid` crate's components: exact point catalogues, the
//! completeness, irredundancy and pilot experiments, and the checks and
//! comparisons of their transcripts. The package depends on the crate and
//! never the reverse, so nothing here changes a decision.
pub mod config;
pub mod points;
pub mod experiments;

#[cfg(test)]
pub(crate) mod testing;
