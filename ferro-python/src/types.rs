//! `PyTrajectory` — `#[pyclass]` 包装 `ferro_core::Trajectory`。
//!
//! 纯 Rust 库对 Python 零感知；所有转换逻辑集中在本 crate。

use ferro_core::Trajectory;
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

/// 轨迹对象（单帧文件也用此类型）。
#[pyclass(name = "Trajectory")]
#[derive(Clone)]
pub struct PyTrajectory {
    pub(crate) inner: Trajectory,
}

impl PyTrajectory {
    pub(crate) fn new(inner: Trajectory) -> Self {
        Self { inner }
    }

    fn frame_ref(&self, index: usize) -> PyResult<&ferro_core::Frame> {
        self.inner
            .frame(index)
            .ok_or_else(|| PyIndexError::new_err(format!("frame index {index} out of range")))
    }
}

#[pymethods]
impl PyTrajectory {
    /// 帧数。
    fn n_frames(&self) -> usize {
        self.inner.n_frames()
    }

    /// 原子数（各帧一致时返回 `Some`，空轨迹返回 `None`）。
    fn n_atoms(&self) -> Option<usize> {
        self.inner.n_atoms()
    }

    /// 首帧出现的元素（按首现顺序去重）。
    fn elements(&self) -> Vec<String> {
        self.inner
            .first()
            .map(|f| f.unique_elements())
            .unwrap_or_default()
    }

    /// 指定帧的原子坐标 \[Å\]，列表 of (x, y, z)。
    fn positions(&self, frame: usize) -> PyResult<Vec<(f64, f64, f64)>> {
        let f = self.frame_ref(frame)?;
        Ok(f.atoms
            .iter()
            .map(|a| (a.position.x, a.position.y, a.position.z))
            .collect())
    }

    /// 指定帧的元素符号列表（与 `positions` 顺序一致）。
    fn symbols(&self, frame: usize) -> PyResult<Vec<String>> {
        let f = self.frame_ref(frame)?;
        Ok(f.atoms.iter().map(|a| a.element.clone()).collect())
    }

    /// 指定帧的晶胞矩阵（行向量 a/b/c，单位 Å）；非周期性返回 `None`。
    fn cell(&self, frame: usize) -> PyResult<Option<Vec<Vec<f64>>>> {
        let f = self.frame_ref(frame)?;
        Ok(f.cell.as_ref().map(|c| {
            (0..3)
                .map(|i| (0..3).map(|j| c.matrix[(i, j)]).collect())
                .collect()
        }))
    }

    /// 指定帧的总电荷。
    fn charge(&self, frame: usize) -> PyResult<i32> {
        Ok(self.frame_ref(frame)?.charge)
    }

    /// 指定帧的自旋多重度 2S+1。
    fn multiplicity(&self, frame: usize) -> PyResult<u32> {
        Ok(self.frame_ref(frame)?.multiplicity)
    }

    /// 仅保留尾部 `n` 帧的新轨迹。
    fn tail(&self, n: usize) -> PyTrajectory {
        PyTrajectory::new(self.inner.tail(n))
    }

    fn __len__(&self) -> usize {
        self.inner.n_frames()
    }

    fn __repr__(&self) -> String {
        format!(
            "Trajectory(frames={}, atoms={})",
            self.inner.n_frames(),
            self.inner
                .n_atoms()
                .map(|n| n.to_string())
                .unwrap_or_else(|| "?".into()),
        )
    }
}
