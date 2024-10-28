use std::{fmt, usize};

use itertools::Itertools;

#[derive(Clone, Copy)]
pub enum BoundaryConditions {
    Periodic,
    SymmetryBreaking,
}

impl fmt::Display for BoundaryConditions {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                BoundaryConditions::Periodic => "Periodic boundary conditions",
                BoundaryConditions::SymmetryBreaking => "Symmetry breaking",
            }
        )
    }
}

#[derive(Clone, Copy)]
pub struct Lattice {
    pub d: usize,
    pub n: usize,
    pub vertices: usize,
    pub edges: usize,
    pub boundary_conditions: BoundaryConditions,
}

impl Lattice {
    pub fn new(d: usize, n: usize, boundary_conditions: BoundaryConditions) -> Lattice {
        let vertices = n.pow(d as u32);
        let edges_cnt = match boundary_conditions {
            BoundaryConditions::Periodic => 2 * vertices,
            BoundaryConditions::SymmetryBreaking => 2 * vertices + 2 * n,
        };
        Lattice {
            d,
            n,
            vertices,
            edges: edges_cnt,
            boundary_conditions: boundary_conditions,
        }
    }
    pub fn idx_to_coords(&self, idx: usize) -> Vec<usize> {
        (0..self.d)
            .scan(idx, |state, _| {
                let (q, r) = (*state / self.n, *state % self.n);
                *state = q;
                Some(r)
            })
            .collect()
    }
    pub fn coords_to_idx(&self, coords: Vec<usize>) -> usize {
        coords
            .iter()
            .enumerate()
            .map(|(i, x)| x * self.n.pow(i as u32))
            .sum()
    }
    pub fn neighbours(&self, idx: usize) -> Vec<usize> {
        let coords = self.idx_to_coords(idx);
        match self.boundary_conditions {
            BoundaryConditions::Periodic => {
                let wrap = |x: i32| {
                    (if x < 0 {
                        x + self.n as i32
                    } else if x >= self.n as i32 {
                        x - self.n as i32
                    } else {
                        x
                    }) as usize
                };
                (0..self.d)
                    .flat_map(|pos| -> Vec<Vec<i32>> {
                        vec![
                            coords
                                .iter()
                                .enumerate()
                                .map(|(i, v)| if i == pos { *v as i32 + 1 } else { *v as i32 })
                                .collect(),
                            coords
                                .iter()
                                .enumerate()
                                .map(|(i, v)| if i == pos { *v as i32 - 1 } else { *v as i32 })
                                .collect(),
                        ]
                    })
                    .map(|neighbour| {
                        self.coords_to_idx(neighbour.iter().map(|x| wrap(*x)).collect())
                    })
                    .collect()
            }
            BoundaryConditions::SymmetryBreaking => (0..self.d)
                .flat_map(|pos| -> Vec<Vec<i32>> {
                    vec![
                        coords
                            .iter()
                            .enumerate()
                            .map(|(i, v)| if i == pos { *v as i32 + 1 } else { *v as i32 })
                            .collect(),
                        coords
                            .iter()
                            .enumerate()
                            .map(|(i, v)| if i == pos { *v as i32 - 1 } else { *v as i32 })
                            .collect(),
                    ]
                })
                .filter(|candidate| candidate.iter().all(|x| 0 <= *x && *x < self.n as i32))
                .map(|neighbour| {
                    self.coords_to_idx(neighbour.iter().map(|x| *x as usize).collect())
                })
                .collect(),
        }
    }
    pub fn neighbouring(&self, idx: usize, spins: &Vec<bool>) -> i32 {
        self.neighbours(idx)
            .into_iter()
            .map(|idx| if spins[idx] { 1 } else { -1 })
            .pad_using(2 * self.d, |_| 1)
            .sum()
    }
}
