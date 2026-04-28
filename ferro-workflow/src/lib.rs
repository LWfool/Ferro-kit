//! ChemTools Workflow - 计算任务创建和管理模块
//! 
//! 用于生成各种量子化学软件的输入文件，如 Gaussian, ORCA, Gromacs 等

pub mod job_builder;
pub mod templates;

pub use job_builder::*;

#[cfg(test)]
mod tests {
    #[test]
    fn test_workflow_module() {
        assert!(true);
    }
}
