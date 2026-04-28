//! 任务模板管理
//! 
//! 存储和管理常用的计算任务模板

/// 获取 Gaussian 优化任务模板
pub fn gaussian_opt_template() -> &'static str {
    "%chk=opt.chk
# B3LYP/6-31G* Opt

Geometry Optimization

0 1
"
}

/// 获取 ORCA 单点能模板
pub fn orca_sp_template() -> &'static str {
    "! B3LYP def2-TZVP

* xyz 0 1
* end
"
}

// 可以添加更多模板...
