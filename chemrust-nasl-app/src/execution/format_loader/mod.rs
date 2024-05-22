use std::{error::Error, fs::read_to_string, path::Path};

use castep_cell_io::{CellParser, LatticeParam};
use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::{
    atom::Atoms,
    geom::coordinates::CoordData,
    lattice::{
        cell_param::{CellConstants, LatticeVectors},
        LatticeCell,
    },
};
use nalgebra::{Matrix3, Point3};

pub fn load_cell_file<P: AsRef<Path>>(cell_path: P) -> Result<LatticeCell, Box<dyn Error>> {
    let content = read_to_string(cell_path)?;
    let model = CellParser::from(content.as_str()).parse()?;
    let lattice_param_block = model.lattice();
    let mut indices = Vec::new();
    let mut symbols: Vec<ElementSymbol> = Vec::new();
    let mut coordinates = Vec::new();
    let mut labels = Vec::new();
    let lattice_vec = match lattice_param_block.parameter() {
        LatticeParam::LatticeCart(lat_cart) => {
            let data = [
                lat_cart.a().to_vec(),
                lat_cart.b().to_vec(),
                lat_cart.c().to_vec(),
            ];
            LatticeVectors::new(Matrix3::from_vec(data.concat()))
        }
        LatticeParam::LatticeABC(lat_abc) => CellConstants::new(
            lat_abc.a(),
            lat_abc.b(),
            lat_abc.c(),
            lat_abc.alpha().value(),
            lat_abc.beta().value(),
            lat_abc.gamma().value(),
        )
        .into(),
    };
    model
        .ionic_positions()
        .positions()
        .iter()
        .enumerate()
        .for_each(|(id, pos)| {
            indices.push(id);
            symbols.push(pos.symbol());
            coordinates.push(CoordData::Fractional(Point3::from(pos.coordinate())));
            labels.push(None);
        });
    let atoms = Atoms::new(indices, symbols, coordinates, labels);
    Ok(LatticeCell::new(lattice_vec, atoms))
}
