//! Figures of the rhombicosidodecahedron problem: shadows of single
//! configurations and of configuration boxes, and configuration-space slices
//! coloured by the component that eliminated each cell. The areas below
//! depend only on the areas listed before them and on the `rid` crate.
pub mod content;
pub mod drawing;
pub mod description;
pub mod classifier;

#[cfg(test)]
mod random;
