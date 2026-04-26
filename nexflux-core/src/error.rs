//! Error type definitions.

use thiserror::Error;

/// ChemTools 核心错误类型
#[derive(Error, Debug)]
pub enum ChemError {
    /// 文件不存在或无法访问
    #[error("File not found: {0}")]
    FileNotFound(String),
    
    /// 文件格式错误
    #[error("Invalid file format: {0}")]
    InvalidFormat(String),
    
    /// 解析错误
    #[error("Parse error: {0}")]
    ParseError(String),
    
    /// 索引越界
    #[error("Index out of bounds: {0}")]
    IndexOutOfBounds(usize),
    
    /// 数据验证错误
    #[error("Validation error: {0}")]
    ValidationError(String),
    
    /// IO 错误
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    
    /// 其他错误
    #[error("Error: {0}")]
    Other(String),
}

/// Result 类型别名
pub type Result<T> = std::result::Result<T, ChemError>;
