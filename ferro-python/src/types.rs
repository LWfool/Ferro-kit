//! `PyTrajectory` — `#[pyclass]` 包装 `ferro_core::Trajectory`。
//!
//! 纯 Rust 库对 Python 零感知；所有转换逻辑集中在本 crate。

use ferro_core::Trajectory;
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::pyerr;

/// 轨迹对象（单帧文件也用此类型）。
///
/// `skip_from_py_object`：所有绑定函数都按 `&PyTrajectory` 借用，没有按值提取的场景。
/// pyo3 0.29 起 `#[pyclass]` + `Clone` 的隐式 `FromPyObject` 需显式选择，
/// 而隐式实现会在每次提取时克隆整条轨迹 —— 对大轨迹是不该发生的深拷贝。
#[pyclass(name = "Trajectory", skip_from_py_object)]
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

    /// 指定帧的总能量 \[eV\]；结构文件没有能量时为 `None`。
    fn energy(&self, frame: usize) -> PyResult<Option<f64>> {
        Ok(self.frame_ref(frame)?.energy)
    }

    /// 指定帧每个原子的力 \[eV/Å\]，顺序与 `positions` 一致；没有时为 `None`。
    fn forces(&self, frame: usize) -> PyResult<Option<Vec<(f64, f64, f64)>>> {
        Ok(self.frame_ref(frame)?.forces.as_ref()
            .map(|fs| fs.iter().map(|f| (f.x, f.y, f.z)).collect()))
    }

    /// 指定帧的应力张量 \[eV/Å³\]，3×3 行优先，**正 = 压缩**（CP2K / VASP / QE 的符号）；
    /// 没有时为 `None`。
    fn stress(&self, frame: usize) -> PyResult<Option<Vec<Vec<f64>>>> {
        Ok(self.frame_ref(frame)?.stress.as_ref().map(|s| {
            (0..3).map(|i| (0..3).map(|j| s[(i, j)]).collect()).collect()
        }))
    }

    /// 指定帧三个方向的周期性 (x, y, z)。
    fn pbc(&self, frame: usize) -> PyResult<(bool, bool, bool)> {
        let [x, y, z] = self.frame_ref(frame)?.pbc;
        Ok((x, y, z))
    }

    /// 仅保留尾部 `n` 帧的新轨迹。
    fn tail(&self, n: usize) -> PyTrajectory {
        PyTrajectory::new(self.inner.tail(n))
    }

    /// 选帧，语义同 `ferro convert` 的 `--start/--end/--stride/--number`：0 起、两端闭区间，
    /// 默认全部帧；`number` 在区间内均匀取（两端都取），与 `stride` 互斥。返回原帧号。
    ///
    /// 参数校验与 `ferro-cli/src/cmd/convert.rs` 的 `run` 各写一份（两处、报错措辞
    /// 一边是 CLI 旗标一边是 Python 参数名，按 R6 不合）；取帧本身都走 ferro-core
    #[pyo3(signature = (start = None, end = None, stride = None, number = None))]
    fn frame_indices(
        &self,
        start: Option<usize>,
        end: Option<usize>,
        stride: Option<usize>,
        number: Option<usize>,
    ) -> PyResult<Vec<usize>> {
        let s = start.unwrap_or(0);
        if let Some(e) = end {
            if s > e {
                return Err(pyerr(format!("start {s} is past end {e} (both 0-based and inclusive)")));
            }
        }
        if stride.is_some() && number.is_some() {
            return Err(pyerr("give stride or number, not both"));
        }
        if stride == Some(0) || number == Some(0) {
            return Err(pyerr("stride and number must be at least 1"));
        }
        let idx = match number {
            Some(n) => self.inner.spread_indices(s, end, n),
            None => self.inner.select_indices(s, end, stride.unwrap_or(1)),
        };
        if idx.is_empty() {
            return Err(pyerr(format!(
                "no frames selected: the trajectory holds {} frame(s), start {s} selects none",
                self.inner.n_frames()
            )));
        }
        Ok(idx)
    }

    /// 按原帧号取出若干帧的新轨迹（越界报错，不静默跳过）。
    fn subset(&self, indices: Vec<usize>) -> PyResult<PyTrajectory> {
        for &i in &indices {
            self.frame_ref(i)?;
        }
        Ok(PyTrajectory::new(self.inner.subset(&indices)))
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
