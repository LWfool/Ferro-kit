//! ferro-workflow — 计算任务输入文件生成
//!
//! 支持 Gaussian、CP2K、Quantum ESPRESSO 输入文件构建。

pub mod job_builder;
pub mod templates;
pub mod cp2k;
pub mod cp2k_basis_db;
pub mod qe;

pub use job_builder::*;
pub use cp2k::{
    Cp2kJobBuilder,
    Cp2kTask, Cp2kFunctional, Cp2kBasis, Cp2kDispersion,
    Cp2kScf, Cp2kPbc, Cp2kCubePrint, Cp2kAtomCharge,
    Cp2kThermostat, Cp2kMdParams,
};
pub use cp2k_basis_db::{DbEntry, DbFunc, DbKind};
pub use qe::{
    QeJobBuilder, QeTask, QeFunctional, QeSmearing, QeMdParams,
};

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Frame};
    use nalgebra::Vector3;

    #[test]
    fn test_builders_constructible() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        assert!(Cp2kJobBuilder::new(f.clone()).build().is_ok());
        assert!(QeJobBuilder::new(f).build().is_ok());
    }
}
