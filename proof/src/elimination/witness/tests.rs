use super::*;
use crate::arithmetic::exact::{frac, q, Q};
use crate::problem::configuration;
use crate::problem::geometry::{edges, screen, symmetries, EDGE_COUNT};
use std::collections::{BTreeMap, BTreeSet};

mod adversarial;

type C = QSqrt5;

fn int(n: i64) -> C {
    C::integer(n)
}

fn rat(x: Q) -> C {
    C::from_rational(x)
}

fn at(p: &Polynomial, x: &[Q; 5]) -> C {
    p.evaluate(&x.clone().map(C::from_rational))
}

/// `1/c` in Q(√5): the conjugate divided by the field norm.
fn inverse(c: &C) -> C {
    c.conjugate().scale(&(q(1) / c.norm()))
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn between(&mut self, lo: &Q, hi: &Q) -> Q {
        let denominator = 1 + self.below(1 << 20) as i64;
        let k = self.below(denominator as u64 + 1) as i64;
        lo + (hi - lo) * frac(k, denominator)
    }
    fn rational(&mut self, bound: i64) -> Q {
        self.between(&q(-bound), &q(bound))
    }
    /// A rational number or, half of the time, an irrational one of Q(√5).
    fn number(&mut self, bound: i64) -> C {
        let a = self.rational(bound);
        if self.below(2) == 0 {
            return rat(a);
        }
        let b = self.rational(bound) / q(3);
        C::new(a / q(2), b)
    }
    fn support(&mut self) -> Support {
        if self.below(3) == 0 {
            let normal = loop {
                let n = [self.number(2), self.number(2)];
                if !n.iter().all(C::is_zero) {
                    break n;
                }
            };
            let contact = self.below(VERTEX_COUNT as u64) as usize;
            Support::Direction(Direction::new(normal, contact).unwrap())
        } else {
            let [a, b] = edges()[self.below(EDGE_COUNT as u64) as usize];
            let (from, to) = if self.below(2) == 0 { (a, b) } else { (b, a) };
            Support::Edge(Edge::new(from, to).unwrap())
        }
    }
    fn gap(&mut self) -> Gap {
        let support = self.support();
        Gap::new(support, self.below(VERTEX_COUNT as u64) as usize).unwrap()
    }
}

/// The affine gap at the configuration `u = (s, t, 1)`, `r`, computed
/// independently of [`Gap::polynomial`], as the tests' reference: the plug
/// vertex is turned by quaternion conjugation with `(1, r)`, both vertices are
/// projected by `T_u`, and the support is read in the screen.
fn value_at(gap: &Gap, view: &[C; 2], rotation: &Point) -> C {
    let vertices = geometry::vertices();
    let u: Point = [view[0].clone(), view[1].clone(), C::one()];
    let [r1, r2, r3] = rotation.clone();
    let turned = geometry::rotate(&[C::one(), r1, r2, r3], &vertices[gap.plug()]);
    let scale = &C::one() + &geometry::dot(rotation, rotation);
    let hole = screen(&u, &vertices[gap.support().contact()]);
    let plug = screen(&u, &turned);
    // With u₃ = 1, n(u)·X is the 2D cross product of T_u(e) and T_u(X) for an
    // edge e, and (n₁, n₂)·T_u(X) for a planar direction.
    let read = |x: &[C; 2]| match gap.support() {
        Support::Edge(edge) => {
            let d = screen(&u, &edge.direction());
            &(&d[0] * &x[1]) - &(&d[1] * &x[0])
        }
        Support::Direction(direction) => {
            let [n1, n2] = direction.screen_normal();
            &(n1 * &x[0]) + &(n2 * &x[1])
        }
    };
    &(&scale * &read(&hole)) - &read(&plug)
}

fn constants(p: &[C]) -> Vec<Polynomial> {
    p.iter().cloned().map(Polynomial::constant).collect()
}

fn triple(p: Vec<Polynomial>) -> [Polynomial; 3] {
    p.try_into().unwrap()
}

fn value_of_constant(p: &Polynomial) -> C {
    p.evaluate(&std::array::from_fn(|_| int(0)))
}

fn affine(s: &C, t: &C) -> Point {
    [s.clone(), t.clone(), int(1)]
}

// ---- The canonical JSON form -------------------------------------------

fn to_json<T: serde::Serialize>(value: &T) -> String {
    serde_json::to_string(value).unwrap()
}

fn from_json<T: serde::de::DeserializeOwned>(text: &str) -> Result<T, serde_json::Error> {
    serde_json::from_str(text)
}

/// Serialised, parsed back and serialised again: the same value and bytes.
fn assert_round_trip(witness: &Witness, text: &str) {
    assert_eq!(to_json(witness), text);
    let parsed: Witness = from_json(text).unwrap();
    assert_eq!(&parsed, witness);
    assert_eq!(to_json(&parsed), text);
}

#[test]
fn witnesses_have_one_exact_json_form() {
    // Every oriented edge of the RID with several plug vertices.
    for &[a, b] in edges() {
        for (from, to) in [(a, b), (b, a)] {
            for plug in [0, from, to, VERTEX_COUNT - 1] {
                let gap = Gap::new(Support::Edge(Edge::new(from, to).unwrap()), plug).unwrap();
                let text = format!("{{\"edge\":[{from},{to}],\"vertex\":{plug}}}");
                assert_eq!(to_json(&gap), text);
                assert_eq!(from_json::<Gap>(&text).unwrap(), gap);
                assert_round_trip(&Witness::Gap(gap), &text);
            }
        }
    }
    // A planar direction: components a + b√5 as canonical rational strings.
    let normal = [int(1), C::new(frac(-1, 2), frac(1, 2))];
    let gap = Gap::new(Support::Direction(Direction::new(normal, 7).unwrap()), 9).unwrap();
    let text = r#"{"direction":[["1","0"],["-1/2","1/2"]],"contact":7,"vertex":9}"#;
    assert_round_trip(&Witness::Gap(gap), text);
    let mut rng = Random(77);
    for _ in 0..300 {
        let witness = Witness::Gap(rng.gap());
        let text = to_json(&witness);
        assert_round_trip(&witness, &text);
    }
    // Domain inequalities: the bare index in a record, `{"inequality":i}` as
    // a witness.
    for index in 0..CONSTRAINT_COUNT {
        let inequality = DomainInequality::new(index).unwrap();
        assert_eq!(to_json(&inequality), index.to_string());
        assert_eq!(from_json::<DomainInequality>(&index.to_string()).unwrap(), inequality);
        let text = format!("{{\"inequality\":{index}}}");
        assert_round_trip(&Witness::Domain(inequality), &text);
    }
}

/// The maximum form's JSON form is an edge gap's: every oriented edge with
/// several plug vertices round-trips byte for byte; a pair that is not an
/// edge and a vertex out of range are refused with the constructor's reason;
/// missing, unknown, repeated, `null` and mistyped fields and the other
/// witness forms are refused.
#[test]
fn the_maximum_form_has_one_exact_json_form() {
    for &[a, b] in edges() {
        for (from, to) in [(a, b), (b, a)] {
            for plug in [0, from, to, VERTEX_COUNT - 1] {
                let maximum = MaximumGap::new(Edge::new(from, to).unwrap(), plug).unwrap();
                let text = format!("{{\"edge\":[{from},{to}],\"vertex\":{plug}}}");
                assert_eq!(to_json(&maximum), text);
                assert_eq!(from_json::<MaximumGap>(&text).unwrap(), maximum);
            }
        }
    }
    let [a, b] = edges()[3];
    let refusal = |text: &str| from_json::<MaximumGap>(text).unwrap_err().to_string();
    let not_an_edge = (0..VERTEX_COUNT).find(|&c| !geometry::is_edge(a, c)).unwrap();
    assert!(refusal(&format!(r#"{{"edge":[{a},{not_an_edge}],"vertex":1}}"#)).contains(&WitnessError::NotAnEdge { from: a, to: not_an_edge }.to_string()));
    assert!(refusal(&format!(r#"{{"edge":[{a},{b}],"vertex":60}}"#)).contains(&WitnessError::VertexOutOfRange { index: 60 }.to_string()));
    assert!(refusal(&format!(r#"{{"edge":[{a},{b}],"vertex":1,"sign":1}}"#)).contains("unknown field"));
    for text in [
        format!(r#"{{"edge":[{a},{b}]}}"#),
        format!(r#"{{"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":1,"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":null}}"#),
        format!(r#"{{"edge":null,"vertex":1}}"#),
        format!(r#"{{"edge":[{a}],"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b},{a}],"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":1.0}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":-1}}"#),
        r#"{"inequality":3}"#.to_string(),
        r#"{"direction":[["1","0"],["0","0"]],"contact":0,"vertex":3}"#.to_string(),
        format!(r#"{{"edge":[{a},{b}],"contact":{a},"vertex":1}}"#),
    ] {
        assert!(from_json::<MaximumGap>(&text).is_err(), "{text}");
    }
    // Serde reads an array positionally; it is not the canonical spelling,
    // so the certificate's byte-exact round trip refuses it in a record.
    let array = format!("[[{a},{b}],1]");
    let read: MaximumGap = from_json(&array).unwrap();
    assert_ne!(to_json(&read), array);
}

#[test]
fn parsed_witnesses_are_validated_like_constructed_ones() {
    // Edges: accepted exactly when the constructor accepts the ordered pair.
    for a in 0..VERTEX_COUNT + 1 {
        for b in 0..VERTEX_COUNT + 1 {
            let parsed = from_json::<Gap>(&format!("{{\"edge\":[{a},{b}],\"vertex\":0}}"));
            assert_eq!(parsed.is_ok(), Edge::new(a, b).is_ok(), "[{a}, {b}]");
        }
    }
    let [a, b] = edges()[0];
    let refused = [
        format!(r#"{{"edge":[{a},{b}],"vertex":60}}"#),
        format!(r#"{{"edge":[{a},{b}]}}"#),
        format!(r#"{{"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":1,"extra":0}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":1,"contact":{a}}}"#),
        format!(r#"{{"edge":[{a},{b},{a}],"vertex":1}}"#),
        format!(r#"{{"edge":[{a}],"vertex":1}}"#),
        format!(r#"{{"edge":[{a},-1],"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":1.0}}"#),
        format!(r#"{{"edge":["{a}","{b}"],"vertex":1}}"#),
        format!(r#"{{"edge":[{a},{b}],"vertex":null}}"#),
        r#"{"direction":[["0","0"],["0","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[["2/4","0"],["1","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[["+1","0"],["1","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[["1","0"],["1/0","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[[1,0],[1,0]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[["1","0"],["1","0"]],"contact":60,"vertex":2}"#.into(),
        r#"{"direction":[["1","0"],["1","0"]],"vertex":2}"#.into(),
        r#"{"direction":[["1","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"direction":[["1","0","0"],["1","0"]],"contact":1,"vertex":2}"#.into(),
        r#"{"inequality":73}"#.into(),
        r#"{"inequality":-1}"#.into(),
        r#"{"inequality":"1"}"#.into(),
        format!(r#"{{"inequality":1,"edge":[{a},{b}]}}"#),
        format!(r#"{{"inequality":1,"vertex":2}}"#),
        "{}".into(),
        "[]".into(),
        "null".into(),
        "13".into(),
        // Serde reads arrays positionally; none has the fields of one form.
        format!("[[{a},{b}]]"),
        format!("[[{a},{b}],null,null,1]"),
        format!(r#"[[{a},{b}],[["1","0"],["1","0"]],1,2]"#),
        "[null,null,null,null,1]".into(),
    ];
    for text in &refused {
        assert!(from_json::<Witness>(text).is_err(), "{text}");
    }
    for text in ["73", "-1", "1.0", "\"1\"", "null"] {
        assert!(from_json::<DomainInequality>(text).is_err(), "{text}");
    }
    // The design's example record shape: [16, 17] is an edge, [16, 58] not.
    assert!(from_json::<Witness>(r#"{"edge":[16,17],"vertex":3}"#).is_ok());
    assert!(from_json::<Witness>(r#"{"edge":[16,58],"vertex":3}"#).is_err());
    // Refusals keep the constructors' reasons where serde passes them on.
    let error = from_json::<DomainInequality>("73").unwrap_err().to_string();
    assert!(error.contains("domain inequality 73 is not below 73"), "{error}");
}

#[test]
fn parsing_refusals_name_their_reason() {
    let witness = |text: &str| from_json::<Witness>(text).unwrap_err().to_string();
    let gap = |text: &str| from_json::<Gap>(text).unwrap_err().to_string();
    for (text, reason) in [
        (r#"{"edge":[16,58],"vertex":3}"#, "[16, 58] is not an edge"),
        (r#"{"edge":[16,17],"vertex":99}"#, "vertex 99 is not below 60"),
        (r#"{"edge":[60,17],"vertex":3}"#, "vertex 60 is not below 60"),
        (r#"{"inequality":73}"#, "domain inequality 73 is not below 73"),
        (r#"{"direction":[["0","0"],["0","-0"]],"contact":1,"vertex":2}"#, "non-canonical"),
        (r#"{"direction":[["0","0"],["0","0"]],"contact":1,"vertex":2}"#, "direction is zero"),
        (r#"{"direction":[["2/4","0"],["1","0"]],"contact":1,"vertex":2}"#, "\"2/4\""),
        (r#"{"direction":[["1","0"],["1","0"]],"contact":60,"vertex":2}"#, "vertex 60 is not"),
        (r#"{"direction":[["1","0"],["1","0"]],"contact":1,"vertex":60}"#, "vertex 60 is not"),
        (r#"{"edge":[16,17]}"#, r#"the fields ["edge"] is not a witness"#),
        (r#"{"inequality":1,"vertex":2}"#, r#"the fields ["vertex", "inequality"]"#),
        ("{}", "the fields [] is not a witness"),
        (r#"{"edge":[16,17],"vertex":3,"extra":0}"#, "unknown field `extra`"),
        (r#"{"edge":[16,17],"vertex":3,"vertex":3}"#, "duplicate field `vertex`"),
        (r#"{"edge":[16,17],"vertex":null}"#, "invalid type: null"),
        (r#"{"inequality":"1"}"#, "invalid type: string"),
    ] {
        let error = witness(text);
        assert!(error.contains(reason), "{text}: {error}");
        assert!(!error.contains("untagged"), "{text}: {error}");
    }
    // A gap names the gap forms, also when given a domain inequality.
    assert!(gap(r#"{"inequality":1}"#).contains(r#"["inequality"] is not a support gap"#));
    assert!(gap(r#"{"edge":[16,58],"vertex":3}"#).contains("[16, 58] is not an edge"));
    assert!(gap(r#"{"direction":[["1","0"]],"contact":1,"vertex":2}"#).contains("length 1"));
    // The same refusals through the typed conversion.
    let fields: Fields = serde_json::from_str(r#"{"edge":[16,17]}"#).unwrap();
    assert_eq!(
        Witness::try_from(fields),
        Err(WitnessError::Form { fields: vec!["edge"], expected: WITNESS_FORMS })
    );
    let fields: Fields = serde_json::from_str(r#"{"edge":[16,58],"vertex":3}"#).unwrap();
    assert_eq!(Witness::try_from(fields), Err(WitnessError::NotAnEdge { from: 16, to: 58 }));
}

// ---- Construction ------------------------------------------------------

#[test]
fn constructors_refuse_invalid_data() {
    for a in 0..VERTEX_COUNT {
        for b in 0..VERTEX_COUNT {
            match Edge::new(a, b) {
                Ok(edge) => {
                    assert!(geometry::is_edge(a, b));
                    assert_eq!((edge.from(), edge.to()), (a, b));
                    assert_eq!(Support::Edge(edge).contact(), a);
                }
                Err(error) => {
                    assert!(!geometry::is_edge(a, b));
                    assert_eq!(error, WitnessError::NotAnEdge { from: a, to: b });
                }
            }
        }
    }
    let [a, b] = edges()[0];
    assert_eq!(
        Edge::new(VERTEX_COUNT, a),
        Err(WitnessError::VertexOutOfRange { index: VERTEX_COUNT })
    );
    assert_eq!(
        Edge::new(b, usize::MAX),
        Err(WitnessError::VertexOutOfRange { index: usize::MAX })
    );
    assert_eq!(
        Direction::new([int(0), int(0)], 3),
        Err(WitnessError::ZeroDirection)
    );
    assert_eq!(
        Direction::new([int(1), int(0)], VERTEX_COUNT),
        Err(WitnessError::VertexOutOfRange { index: VERTEX_COUNT })
    );
    let direction = Direction::new([int(0), C::new(q(0), frac(1, 7))], 5).unwrap();
    assert_eq!(direction.contact(), 5);
    assert_eq!(Support::Direction(direction.clone()).contact(), 5);
    let support = Support::Direction(direction);
    assert_eq!(
        Gap::new(support.clone(), VERTEX_COUNT),
        Err(WitnessError::VertexOutOfRange { index: VERTEX_COUNT })
    );
    let gap = Gap::new(support.clone(), 59).unwrap();
    assert_eq!((gap.support(), gap.plug()), (&support, 59));
    assert_eq!(
        DomainInequality::new(CONSTRAINT_COUNT),
        Err(WitnessError::InequalityOutOfRange { index: CONSTRAINT_COUNT })
    );
    assert_eq!(DomainInequality::new(72).unwrap().index(), 72);
}

// ---- Gap polynomials against direct projection -------------------------

#[test]
fn gap_polynomials_match_direct_projection_at_random_rational_configurations() {
    let mut random = Random(41);
    let (view, rotation) = configuration::coordinates();
    for _ in 0..60 {
        let gap = random.gap();
        let polynomial = gap.polynomial(&view, &rotation).unwrap();
        for _ in 0..4 {
            let x: [Q; 5] = std::array::from_fn(|_| random.rational(1));
            let direct = value_at(
                &gap,
                &[rat(x[0].clone()), rat(x[1].clone())],
                &[rat(x[2].clone()), rat(x[3].clone()), rat(x[4].clone())],
            );
            assert_eq!(at(&polynomial, &x), direct, "{gap:?}");
        }
    }
}

#[test]
fn gap_polynomials_match_direct_projection_at_irrational_configurations_and_scalings() {
    let mut random = Random(43);
    for _ in 0..80 {
        let gap = random.gap();
        let (s, t) = (random.number(1), random.number(1));
        let r: Point = std::array::from_fn(|_| random.number(1));
        let direct = value_at(&gap, &[s.clone(), t.clone()], &r);
        let rotation = triple(constants(&r));
        let at_affine = gap
            .polynomial(
                &View::Affine {
                    s: Polynomial::constant(s.clone()),
                    t: Polynomial::constant(t.clone()),
                },
                &rotation,
            )
            .unwrap();
        assert_eq!(value_of_constant(&at_affine), direct);
        let lambda = random.number(3);
        let scaled = affine(&s, &t).map(|c| &c * &lambda);
        let at_scaled = gap
            .polynomial(&View::Homogeneous(triple(constants(&scaled))), &rotation)
            .unwrap();
        assert_eq!(value_of_constant(&at_scaled), &direct * &lambda);
    }
}

#[test]
fn gap_polynomials_accept_polynomial_views_and_rotations() {
    let mut random = Random(47);
    // A homogeneous view (x₀, x₁, x₂) and a nonlinear rotation vector.
    let x = Polynomial::variables();
    let view = View::Homogeneous([x[0].clone(), x[1].clone(), x[2].clone()]);
    let rotation = [
        x[3].mul(&x[4]).unwrap(),
        &x[4] - &Polynomial::constant(rat(frac(1, 5))),
        &x[0] + &x[3],
    ];
    for _ in 0..20 {
        let gap = random.gap();
        let polynomial = gap.polynomial(&view, &rotation).unwrap();
        for _ in 0..3 {
            let x: [Q; 5] = std::array::from_fn(|_| random.rational(1));
            if x[2] == q(0) {
                continue;
            }
            let r = [
                rat(&x[3] * &x[4]),
                rat(&x[4] - frac(1, 5)),
                rat(&x[0] + &x[3]),
            ];
            let direct = value_at(&gap, &[rat(&x[0] / &x[2]), rat(&x[1] / &x[2])], &r);
            assert_eq!(at(&polynomial, &x), direct.scale(&x[2]));
        }
    }
}

// ---- Validity against brute force in the screen ------------------------

/// The oriented edge supports the projected hole at the affine view `u`: its
/// projection is not a point and no projected vertex lies to its left.
fn brute_edge(u: &Point, from: usize, to: usize) -> bool {
    let v = geometry::vertices();
    let p = |i: usize| screen(u, &v[i]);
    let (a, b) = (p(from), p(to));
    let d = [&b[0] - &a[0], &b[1] - &a[1]];
    if d.iter().all(C::is_zero) {
        return false;
    }
    (0..VERTEX_COUNT).all(|w| {
        let x = p(w);
        let e = [&x[0] - &a[0], &x[1] - &a[1]];
        (&(&d[0] * &e[1]) - &(&d[1] * &e[0])).sign() != Ordering::Greater
    })
}

/// The contact maximises the planar direction over the projected hole.
fn brute_direction(u: &Point, normal: &[C; 2], contact: usize) -> bool {
    let v = geometry::vertices();
    let read = |i: usize| {
        let x = screen(u, &v[i]);
        &(&normal[0] * &x[0]) + &(&normal[1] * &x[1])
    };
    let h = read(contact);
    (0..VERTEX_COUNT).all(|w| h >= read(w))
}

fn oriented_edges() -> Vec<(usize, usize)> {
    edges().iter().flat_map(|[a, b]| [(*a, *b), (*b, *a)]).collect()
}

fn valid_edges_at(u: &Point) -> Vec<(usize, usize)> {
    oriented_edges()
        .into_iter()
        .filter(|(from, to)| {
            let support = Support::Edge(Edge::new(*from, *to).unwrap());
            let valid = support.check_valid(std::slice::from_ref(u));
            assert_eq!(valid.is_ok(), brute_edge(u, *from, *to), "{from}->{to} at {u:?}");
            valid.is_ok()
        })
        .collect()
}

/// At a generic view the valid oriented edges form one directed cycle.
fn assert_single_cycle(valid: &[(usize, usize)]) {
    let next: BTreeMap<usize, usize> = valid.iter().copied().collect();
    assert_eq!(next.len(), valid.len(), "a vertex with two outgoing supports");
    let start = valid[0].0;
    let (mut current, mut steps) = (start, 0);
    loop {
        current = next[&current];
        steps += 1;
        if current == start {
            break;
        }
        assert!(steps <= valid.len());
    }
    assert_eq!(steps, valid.len());
}

#[test]
fn edge_validity_agrees_with_brute_force_and_traces_the_shadow_outline() {
    let mut random = Random(53);
    for _ in 0..8 {
        let u = affine(&random.number(1), &random.number(1));
        let valid = valid_edges_at(&u);
        assert!(valid.len() >= 10, "{}", valid.len());
        assert_single_cycle(&valid);
        // The same supports are valid at positive multiples of the view.
        let lambda = rat(random.between(&frac(1, 100), &q(100)));
        let scaled = u.clone().map(|c| &c * &lambda);
        for (from, to) in &valid {
            let support = Support::Edge(Edge::new(*from, *to).unwrap());
            assert_eq!(support.check_valid(&[scaled.clone()]), Ok(()));
            // Reversing a supporting edge always fails at a generic view.
            let reversed = Support::Edge(Edge::new(*to, *from).unwrap());
            let Err(WitnessError::Invalid { view: 0, vertex }) = reversed.check_valid(&[u.clone()])
            else {
                panic!("reversed edge accepted");
            };
            let n = reversed.normal_at(&u);
            let v = geometry::vertices();
            let h = geometry::dot(&n, &v[*to]);
            assert!(h < geometry::dot(&n, &v[vertex]));
        }
    }
}

#[test]
fn edge_validity_agrees_with_brute_force_at_boundary_views() {
    let phi = C::new(frac(1, 2), frac(1, 2));
    let inverse_phi = &phi - &int(1);
    let inverse_phi_squared = &int(2) - &phi;
    let half = rat(frac(1, 2));
    let views = [
        affine(&int(0), &int(0)),
        affine(&inverse_phi, &int(0)),
        affine(&int(0), &inverse_phi_squared),
        affine(&(&inverse_phi * &half), &int(0)),
        affine(&int(0), &(&inverse_phi_squared * &half)),
        affine(&(&inverse_phi * &half), &(&inverse_phi_squared * &half)),
        affine(&rat(frac(2, 3)), &rat(frac(2, 5))),
    ];
    for u in &views {
        let valid = valid_edges_at(u);
        assert!(!valid.is_empty());
        // Every projected outline vertex starts some supporting edge.
        let starts: BTreeSet<usize> = valid.iter().map(|e| e.0).collect();
        let ends: BTreeSet<usize> = valid.iter().map(|e| e.1).collect();
        assert_eq!(starts, ends);
    }
}

#[test]
fn direction_validity_agrees_with_brute_force_including_ties() {
    let mut random = Random(59);
    let axis_views = [affine(&int(0), &int(0)), affine(&rat(frac(1, 3)), &int(0))];
    for round in 0..40 {
        let u = if round < 2 {
            axis_views[round].clone()
        } else {
            affine(&random.number(1), &random.number(1))
        };
        let normal = match round % 4 {
            0 => [int(1), int(0)],
            1 => [int(0), int(-1)],
            _ => [random.number(2), random.number(2)],
        };
        let mut valid = 0;
        for contact in 0..VERTEX_COUNT {
            let support = Support::Direction(Direction::new(normal.clone(), contact).unwrap());
            let result = support.check_valid(&[u.clone()]);
            assert_eq!(result.is_ok(), brute_direction(&u, &normal, contact));
            if let Err(error) = result {
                let WitnessError::Invalid { view: 0, vertex } = error else {
                    panic!("unexpected {error:?}");
                };
                // The reported vertex lies strictly beyond the contact.
                let v = geometry::vertices();
                let read = |i: usize| {
                    let x = screen(&u, &v[i]);
                    &(&normal[0] * &x[0]) + &(&normal[1] * &x[1])
                };
                assert!(read(contact) < read(vertex));
            } else {
                valid += 1;
            }
        }
        assert!(valid >= 1);
    }
    // Along the z axis the direction (1, 0) is supported by a whole face:
    // the vertices with the largest x coordinate tie.
    let u = affine(&int(0), &int(0));
    let ties = (0..VERTEX_COUNT)
        .filter(|&h| {
            Support::Direction(Direction::new([int(1), int(0)], h).unwrap())
                .check_valid(&[u.clone()])
                .is_ok()
        })
        .count();
    assert!(ties >= 2, "{ties}");
}

#[test]
fn validity_on_view_rectangles_extends_to_their_interiors() {
    let mut random = Random(61);
    let mut accepted = 0;
    let mut refused = 0;
    for _ in 0..6 {
        let s = random.between(&q(0), &frac(1, 2));
        let t = random.between(&q(0), &frac(1, 3));
        let width = random.between(&frac(1, 10_000), &frac(1, 20));
        let corners: Vec<Point> = [(0, 0), (1, 0), (0, 1), (1, 1)]
            .iter()
            .map(|(i, j)| {
                affine(
                    &rat(&s + &width * q(*i)),
                    &rat(&t + &width * q(*j)),
                )
            })
            .collect();
        let centre = affine(&rat(&s + &width / q(2)), &rat(&t + &width / q(2)));
        for (from, to) in valid_edges_at(&centre) {
            let support = Support::Edge(Edge::new(from, to).unwrap());
            match support.check_valid(&corners) {
                Ok(()) => {
                    accepted += 1;
                    for _ in 0..4 {
                        let a = random.between(&q(0), &q(1));
                        let b = random.between(&q(0), &q(1));
                        let u = affine(&rat(&s + &width * &a), &rat(&t + &width * &b));
                        assert!(brute_edge(&u, from, to));
                    }
                    // A cone of positive multiples of the corners is as valid.
                    let cone: Vec<Point> = corners
                        .iter()
                        .map(|u| {
                            let lambda = rat(random.between(&frac(1, 10), &q(10)));
                            u.clone().map(|c| &c * &lambda)
                        })
                        .collect();
                    assert_eq!(support.check_valid(&cone), Ok(()));
                }
                Err(WitnessError::Invalid { view, vertex }) => {
                    refused += 1;
                    assert!(!brute_edge(&corners[view], from, to));
                    let n = support.normal_at(&corners[view]);
                    let v = geometry::vertices();
                    assert!(geometry::dot(&n, &v[from]) < geometry::dot(&n, &v[vertex]));
                    assert!((0..view).all(|k| brute_edge(&corners[k], from, to)));
                }
                Err(error) => panic!("unexpected {error:?}"),
            }
        }
    }
    assert!(accepted > 0 && refused > 0, "{accepted} {refused}");
}

#[test]
fn degenerate_view_sets_are_refused() {
    let [a, b] = edges()[0];
    let support = Support::Edge(Edge::new(a, b).unwrap());
    let good = affine(&rat(frac(1, 10)), &rat(frac(1, 10)));
    assert_eq!(support.check_valid(&[]), Err(WitnessError::NoViews));
    for (bad, position) in [
        ([int(0), int(0), int(0)], 1),
        ([rat(frac(1, 10)), int(0), int(-1)], 0),
        ([int(0), int(1), C::new(q(2), q(-1))], 1),
    ] {
        let views = if position == 0 {
            vec![bad.clone(), good.clone()]
        } else {
            vec![good.clone(), bad.clone()]
        };
        assert_eq!(
            support.check_valid(&views),
            Err(WitnessError::ViewNotPositive { view: position })
        );
        let domain = Witness::Domain(DomainInequality::new(0).unwrap());
        assert_eq!(
            domain.check_valid(&views),
            Err(WitnessError::ViewNotPositive { view: position })
        );
    }
}

#[test]
fn views_along_an_edge_have_a_zero_normal_and_are_refused() {
    let v = geometry::vertices();
    let mut seen = 0;
    for (from, to) in oriented_edges() {
        let e = geometry::difference(&v[to], &v[from]);
        if e[2].sign() != Ordering::Greater {
            continue;
        }
        seen += 1;
        let support = Support::Edge(Edge::new(from, to).unwrap());
        let along = e.clone();
        let other = affine(&rat(frac(1, 7)), &rat(frac(1, 9)));
        assert!(support.normal_at(&along).iter().all(C::is_zero));
        // The check stops at the first failing view.
        assert_eq!(
            support.check_valid(&[along.clone(), other.clone()]),
            Err(WitnessError::ZeroNormal { view: 0 })
        );
        assert_ne!(
            support.check_valid(&[other.clone(), along.clone()]),
            Ok(())
        );
        assert!(!brute_edge(&along.map(|c| &c * &inverse(&e[2])), from, to));
    }
    assert!(seen > 0);
}

// ---- Necessity at configurations with identical shadows ----------------

#[test]
fn valid_gaps_are_nonnegative_where_the_plug_is_a_symmetric_copy_of_the_hole() {
    let mut random = Random(67);
    // r = 0 and r = g⃗/g₀ for group elements with g₀ ≠ 0 give plug = hole.
    let mut rotations: Vec<Point> = vec![[int(0), int(0), int(0)]];
    for g in symmetries().iter().filter(|g| !g[0].is_zero()).take(12) {
        let w = inverse(&g[0]);
        rotations.push([&g[1] * &w, &g[2] * &w, &g[3] * &w]);
    }
    for _ in 0..3 {
        let (s, t) = (random.number(1), random.number(1));
        let u = affine(&s, &t);
        let valid = valid_edges_at(&u);
        for r in &rotations {
            for (from, to) in valid.iter().take(4) {
                for plug in 0..VERTEX_COUNT {
                    let support = Support::Edge(Edge::new(*from, *to).unwrap());
                    let gap = Gap::new(support, plug).unwrap();
                    let value = value_at(&gap, &[s.clone(), t.clone()], r);
                    assert_ne!(value.sign(), Ordering::Less);
                }
            }
        }
    }
}

// ---- The maximum form -----------------------------------------------------

/// Member `w` of the maximum form at `u = (s, t, 1)`, `r`, by the tests'
/// direct reading: `(1 + |r|²) n(u)·w − n(u)·R̂(r)p` with the plug vertex
/// turned by quaternion conjugation and both points read in the screen by
/// `T_u` (as in [`value_at`]).
fn member_at(maximum: &MaximumGap, w: usize, view: &[C; 2], rotation: &Point) -> C {
    let vertices = geometry::vertices();
    let u: Point = [view[0].clone(), view[1].clone(), C::one()];
    let [r1, r2, r3] = rotation.clone();
    let turned = geometry::rotate(&[C::one(), r1, r2, r3], &vertices[maximum.plug()]);
    let scale = &C::one() + &geometry::dot(rotation, rotation);
    let d = screen(&u, &maximum.edge().direction());
    let read = |x: &[C; 2]| &(&d[0] * &x[1]) - &(&d[1] * &x[0]);
    &(&scale * &read(&screen(&u, &vertices[w]))) - &read(&screen(&u, &turned))
}

fn random_maximum(random: &mut Random) -> MaximumGap {
    let [a, b] = edges()[random.below(EDGE_COUNT as u64) as usize];
    let (from, to) = if random.below(2) == 0 { (a, b) } else { (b, a) };
    MaximumGap::new(Edge::new(from, to).unwrap(), random.below(VERTEX_COUNT as u64) as usize).unwrap()
}

#[test]
fn maximum_members_match_direct_projection() {
    let mut random = Random(91);
    let (view, rotation) = configuration::coordinates();
    for _ in 0..12 {
        let maximum = random_maximum(&mut random);
        let members = maximum.members(&view, &rotation).unwrap();
        assert_eq!(members.len(), VERTEX_COUNT);
        // The member of the edge's first vertex is the edge's gap.
        let gap = Gap::new(Support::Edge(maximum.edge().clone()), maximum.plug()).unwrap();
        assert_eq!(members[maximum.edge().from()], gap.polynomial(&view, &rotation).unwrap());
        for _ in 0..3 {
            let x: [C; 5] = std::array::from_fn(|_| random.number(1));
            let (st, r) = ([x[0].clone(), x[1].clone()], [x[2].clone(), x[3].clone(), x[4].clone()]);
            for (w, member) in members.iter().enumerate() {
                assert_eq!(member.evaluate(&x), member_at(&maximum, w, &st, &r), "member {w}");
            }
        }
    }
    let [a, b] = edges()[0];
    assert_eq!(MaximumGap::new(Edge::new(a, b).unwrap(), 60), Err(WitnessError::VertexOutOfRange { index: 60 }));
}

/// With no validity hypothesis: at `r = 0` and at the group's Cayley vectors
/// the plug's vertex set is the hole's, so for every oriented edge (valid at
/// the view or not) and every plug vertex the maximum over the members is
/// nonnegative, and over the plug vertices its minimum is exactly 0. (The
/// members share `n(u)`, so the maximum is `(1 + |r|²) max_w n·w − n·R̂p`,
/// read in the screen as in [`member_at`].)
#[test]
fn the_maximum_form_is_nonnegative_where_the_plug_is_a_symmetric_copy_of_the_hole() {
    let mut random = Random(92);
    let mut rotations: Vec<Point> = vec![[int(0), int(0), int(0)]];
    for g in symmetries().iter().filter(|g| !g[0].is_zero()).take(6) {
        let w = inverse(&g[0]);
        rotations.push([&g[1] * &w, &g[2] * &w, &g[3] * &w]);
    }
    let vertices = geometry::vertices();
    for _ in 0..2 {
        let u: Point = [random.number(1), random.number(1), C::one()];
        let hole: Vec<[C; 2]> = vertices.iter().map(|w| screen(&u, w)).collect();
        for r in &rotations {
            let [r1, r2, r3] = r.clone();
            let turned: Vec<[C; 2]> =
                vertices.iter().map(|p| screen(&u, &geometry::rotate(&[C::one(), r1.clone(), r2.clone(), r3.clone()], p))).collect();
            let scale = &C::one() + &geometry::dot(r, r);
            for (from, to) in oriented_edges() {
                let d = screen(&u, &Edge::new(from, to).unwrap().direction());
                let read = |x: &[C; 2]| &(&d[0] * &x[1]) - &(&d[1] * &x[0]);
                let support = &scale * &hole.iter().map(read).max().unwrap();
                let values: Vec<C> = turned.iter().map(|x| &support - &read(x)).collect();
                assert!(values.iter().all(|v| v.sign() != Ordering::Less), "{from} → {to}");
                assert!(values.iter().min().unwrap().is_zero(), "{from} → {to}");
            }
        }
    }
}

#[test]
fn domain_witnesses_are_minus_the_domain_polynomials() {
    let (view, rotation) = configuration::coordinates();
    let x = Polynomial::variables();
    let homogeneous = View::Homogeneous([x[0].clone(), x[1].clone(), x[2].clone()]);
    for index in 0..CONSTRAINT_COUNT {
        let witness = Witness::Domain(DomainInequality::new(index).unwrap());
        for v in [&view, &homogeneous] {
            let expected = -&domain::polynomial(index, v, &rotation).unwrap();
            assert_eq!(witness.polynomial(v, &rotation), Ok(expected));
        }
        assert_eq!(witness.check_valid(&[affine(&int(1), &int(1))]), Ok(()));
        assert_eq!(witness.check_valid(&[]), Err(WitnessError::NoViews));
    }
    let mut random = Random(71);
    for _ in 0..10 {
        let gap = random.gap();
        let witness = Witness::Gap(gap.clone());
        assert_eq!(witness.polynomial(&view, &rotation), gap.polynomial(&view, &rotation));
        let u = affine(&random.number(1), &random.number(1));
        assert_eq!(witness.check_valid(&[u.clone()]), gap.support().check_valid(&[u]));
    }
}

#[test]
fn direction_normals_are_perpendicular_to_the_view_and_read_the_screen() {
    let mut random = Random(73);
    for _ in 0..40 {
        let support = random.support();
        let u: Point = [random.number(2), random.number(2), random.number(2)];
        let n = support.normal_at(&u);
        assert!(geometry::dot(&n, &u).is_zero());
        let x: Point = std::array::from_fn(|_| random.number(3));
        if let Support::Direction(direction) = &support {
            let [n1, n2] = direction.screen_normal();
            let image = screen(&u, &x);
            assert_eq!(geometry::dot(&n, &x), &(n1 * &image[0]) + &(n2 * &image[1]));
        }
        // The polynomial normal agrees with the exact one.
        let polynomial = support.normal(&triple(constants(&u)));
        for j in 0..3 {
            assert_eq!(value_of_constant(&polynomial[j]), n[j]);
        }
    }
}

/// Heavy campaign: `cargo test --release --lib -- --ignored
/// elimination::witness::tests::gap_campaign --test-threads 2`.
#[test]
#[ignore]
fn gap_campaign() {
    let mut random = Random(1_000_033);
    let (view, rotation) = configuration::coordinates();
    for _ in 0..3000 {
        let gap = random.gap();
        let polynomial = gap.polynomial(&view, &rotation).unwrap();
        let x: [Q; 5] = std::array::from_fn(|_| random.rational(1));
        let direct = value_at(
            &gap,
            &[rat(x[0].clone()), rat(x[1].clone())],
            &[rat(x[2].clone()), rat(x[3].clone()), rat(x[4].clone())],
        );
        assert_eq!(at(&polynomial, &x), direct);
    }
    for _ in 0..60 {
        let u = affine(&random.number(1), &random.number(1));
        let valid = valid_edges_at(&u);
        assert_single_cycle(&valid);
    }
}

/// Costs of the per-box operations, printed:
/// `cargo test --release --lib -- --ignored elimination::witness::tests::timing
/// --nocapture --test-threads 1`.
#[test]
#[ignore]
fn timing() {
    use crate::problem::configuration::ConfigurationBox;
    use std::time::Instant;
    let mut random = Random(97);
    let (view, rotation) = configuration::coordinates();
    let gaps: Vec<Gap> = (0..400).map(|_| random.gap()).collect();
    let start = Instant::now();
    for gap in &gaps {
        gap.polynomial(&view, &rotation).unwrap();
    }
    println!("gap polynomial (identity zoom): {:?} each", start.elapsed() / 400);
    let boxes: Vec<ConfigurationBox> = (0..10)
        .map(|k| {
            let path: String = (0..20 + 3 * k)
                .map(|_| if random.below(2) == 0 { '0' } else { '1' })
                .collect();
            ConfigurationBox::from_path(&path).unwrap()
        })
        .collect();
    let start = Instant::now();
    let mut accepted = 0;
    for b in &boxes {
        let corners = b.view_corners();
        for (from, to) in oriented_edges() {
            let support = Support::Edge(Edge::new(from, to).unwrap());
            accepted += usize::from(support.check_valid(&corners).is_ok());
        }
    }
    let checks = (boxes.len() * 2 * EDGE_COUNT) as u32;
    println!(
        "edge validity at 4 corner views: {:?} per support ({accepted} of {checks} valid)",
        start.elapsed() / checks
    );
}
