//! Exact search for, and checking of, a certificate that the
//! rhombicosidodecahedron is not Rupert. The areas below depend only on the
//! areas listed before them. The mathematics is in the accompanying article.
pub mod arithmetic;
pub mod problem;
pub mod elimination;
pub mod components;
pub mod search;

/// Identifies this implementation, its fixed data and pinned dependencies.
pub const POLICY: &str = env!("RID_POLICY");
