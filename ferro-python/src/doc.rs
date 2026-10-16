//! `ferro-py doc` 的手册：`docs/src/python/` 经 `include_str!` 编进扩展。
//!
//! wheel 在超算上离线安装，旁边没有源码树，也没有 Rust 的 `ferro` 命令，所以手册
//! 只能随扩展走。渲染用 `ferro-render`，与 `ferro doc` 同一个渲染器。
//!
//! 页面表与 `ferro-cli/src/doc.rs` 的 `PAGES` 是两张表、各收各的页（用户定：Rust 的
//! `ferro doc` 不收 ferro-py 的页，免得一份手册两处收录而走样）。主题解析、列表、分页器
//! 在 `cli.py`，与 `ferro-cli/src/doc.rs` 的同名逻辑两处相似，按 R6 不合并：一边 Rust、
//! 一边 Python，合并只能经 FFI 绕一圈。

use pyo3::prelude::*;

/// (主题, 源文件, 全文)，SUMMARY.md 的顺序 —— `ferro-py doc` 不带参数时照此列出。
/// 主题跟着子命令走（`ferro-py doc job vasp`），帮助页末行就是下一条要敲的命令。
const PAGES: &[(&str, &str, &str)] = &[
    ("api", "python/api.md", include_str!("../../docs/src/python/api.md")),
    ("config", "python/config.md", include_str!("../../docs/src/python/config.md")),
    ("job", "python/job.md", include_str!("../../docs/src/python/job.md")),
    ("job vasp", "python/job-vasp.md", include_str!("../../docs/src/python/job-vasp.md")),
    ("job qe", "python/job-qe.md", include_str!("../../docs/src/python/job-qe.md")),
    ("job cp2k", "python/job-cp2k.md", include_str!("../../docs/src/python/job-cp2k.md")),
    ("conv", "python/conv.md", include_str!("../../docs/src/python/conv.md")),
];

/// 全部手册页：[(主题, 源文件, 全文)]。
#[pyfunction]
#[pyo3(name = "_doc_pages")]
fn doc_pages() -> Vec<(&'static str, &'static str, &'static str)> {
    PAGES.to_vec()
}

/// 手册 markdown → 终端文本（ANSI 与框线字符；`NO_COLOR` 关颜色）。
/// `width` 缺省取 stdout 所在终端的列数，再退回 `$COLUMNS`、80。
#[pyfunction]
#[pyo3(name = "_render_markdown", signature = (text, width = None))]
fn render_markdown(text: &str, width: Option<usize>) -> String {
    let width = width.unwrap_or_else(ferro_render::terminal_width);
    ferro_render::render(text, width, ferro_render::Style::for_terminal())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(doc_pages, m)?)?;
    m.add_function(wrap_pyfunction!(render_markdown, m)?)?;
    Ok(())
}
