//! 把构建信息导出给 `ferro --version`（`FERRO_BUILD_*` 编译期环境变量）。
//!
//! 不引入依赖：日期由 `SystemTime` 自行换算 UTC，git / rustc 走子进程，失败一律给 `unknown`。

use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");

    // 跟踪全部 workspace 成员的源码：它们任何一处改动本来就会让 ferro-cli 重编，
    // 跟踪它们不增加编译次数，却让日期与 commit 不至于停在上一次跑脚本的时刻。
    // 不存在的路径不登记 —— cargo 会把缺失的 rerun-if-changed 视为每次都脏
    let mut watched = vec!["Cargo.toml".to_string(), "Cargo.lock".to_string(), "ferro-cli/build.rs".to_string()];
    for krate in ["ferro-core", "ferro-io", "ferro-structure", "ferro-analysis", "ferro-workflow", "ferro-cli"] {
        watched.push(format!("{krate}/src"));
        watched.push(format!("{krate}/Cargo.toml"));
    }
    // HEAD 管切分支，refs 管提交，index 管 add（dirty 标记随之变化）
    watched.extend([".git/HEAD", ".git/index", ".git/refs"].map(String::from));
    for p in &watched {
        let path = root.join(p);
        if path.exists() {
            println!("cargo:rerun-if-changed={}", path.display());
        }
    }
    println!("cargo:rerun-if-env-changed=SOURCE_DATE_EPOCH");

    let commit = match run("git", &["rev-parse", "--short", "HEAD"], &root) {
        Some(hash) => {
            let dirty = run("git", &["status", "--porcelain", "--untracked-files=no"], &root)
                .is_some_and(|s| !s.is_empty());
            if dirty { format!("{hash} (dirty)") } else { hash }
        }
        None => "unknown".into(),
    };

    // SOURCE_DATE_EPOCH 优先：可复现构建的通行约定
    let secs = std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.parse::<i64>().ok())
        .unwrap_or_else(|| SystemTime::now().duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs() as i64));

    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    let rustc_version = run(&rustc, &["--version"], &root).unwrap_or_else(|| "unknown".into());

    // cargo 给构建脚本的 CARGO_FEATURE_<NAME>：大写、`-` 换成 `_`，原名已不可逆，按小写写回
    let mut features: Vec<String> = std::env::vars()
        .filter_map(|(k, _)| k.strip_prefix("CARGO_FEATURE_").map(|f| f.to_lowercase()))
        .collect();
    features.sort();
    let features = if features.is_empty() { "none".into() } else { features.join(", ") };

    let env = |k: &str| std::env::var(k).unwrap_or_else(|_| "unknown".into());
    println!("cargo:rustc-env=FERRO_BUILD_COMMIT={commit}");
    println!("cargo:rustc-env=FERRO_BUILD_DATE={}", utc_string(secs));
    println!("cargo:rustc-env=FERRO_BUILD_PROFILE={} (opt-level {})", env("PROFILE"), env("OPT_LEVEL"));
    println!("cargo:rustc-env=FERRO_BUILD_TARGET={}", env("TARGET"));
    println!("cargo:rustc-env=FERRO_BUILD_RUSTC={rustc_version}");
    println!("cargo:rustc-env=FERRO_BUILD_FEATURES={features}");
}

fn run(cmd: &str, args: &[&str], dir: &Path) -> Option<String> {
    let out = Command::new(cmd).args(args).current_dir(dir).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Unix 秒 → `YYYY-MM-DD HH:MM UTC`。日期部分是 Howard Hinnant 的 civil_from_days
fn utc_string(secs: i64) -> String {
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02} {:02}:{:02} UTC", rem / 3600, rem % 3600 / 60)
}
