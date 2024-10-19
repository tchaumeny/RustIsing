use image::{ImageBuffer, Rgb};
use plotters::prelude::*;
use slugify::slugify;
use std::path::Path;

use crate::lattice::Lattice;


#[derive(Clone)]
pub struct Series {
    pub name: String,
    pub color: RGBColor,
    pub data: Vec<(i32, i32)>,
}

pub fn plot_histogram(series: Series) {
    let slug = slugify!(&series.name);
    let file_name = format!("{slug}.png");
    let root_area = BitMapBackend::new(&file_name, (800, 600)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let x_min = series.data.iter().min_by_key(|(x, _)| x).unwrap().0 - 1;
    let x_max = series.data.iter().max_by_key(|(x, _)| x).unwrap().0 + 2;
    let w_max = series.data.iter().max_by_key(|(_, w)| w).unwrap().1 + 2;

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(series.name, ("sans-serif", 20))
        .build_cartesian_2d(x_min..x_max, 0..w_max)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();

    ctx.draw_series(series.data.iter().map(|(x, w)| {
        let mut bar = Rectangle::new([(*x - 1, 0), (*x + 1, *w)], series.color.filled());
        bar.set_margin(0, 0, 5, 5);
        bar
    }))
    .unwrap();
    println!("Generated file {}", &file_name);
}


pub fn plot_lattice(name: String, lattice: Lattice, spins: &Vec<bool>) {
    if lattice.d != 2 {
        eprintln!("Cannot plot lattice in dimension != 2.");
        return;
    }
    let slug = slugify!(&name);
    let file_name = format!("{slug}.png");
    let target_width = 500_u32;
    let unit = u32::max(target_width / lattice.n as u32, 1);
    let real_width = unit * lattice.n as u32;
    let mut img = ImageBuffer::new(real_width, real_width);
    for (idx, spin) in spins.into_iter().enumerate() {
        let coords = lattice.idx_to_coords(idx);
        let color: Rgb<u8> = if *spin {
            Rgb([0, 105, 92])
        } else {
            Rgb([128, 203, 196])
        };
        for (x, y) in itertools::iproduct!(0..unit, 0..unit) {
            img.put_pixel(unit * coords[0] as u32 + x, unit * coords[1] as u32 + y, color);
        }
    }

    let path = Path::new(&file_name);
    img.save(path).unwrap();

    println!("Lattice image saved to {:?}", path);
}
