//! `ferro` — the single entry point.
//!
//! One binary, subcommands grouped by what they produce. The old `fe-*` binaries are
//! gone: they were separate front-ends over code that a future script/REPL mode has to
//! link in anyway (a stateful `read` → `gr` → `sq` session cannot shell out to another
//! process without re-reading the trajectory each time), so keeping them would have
//! meant maintaining two copies of every command surface.

use anyhow::Result;
use clap::{Parser, Subcommand};

use ferro::cmd;
use ferro::doc;
use ferro::help;

#[derive(Parser)]
#[command(
    name = "ferro",
    version,
    about = "Computational chemistry toolkit for periodic systems",
    disable_help_subcommand = true,
    arg_required_else_help = false,
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand)]
enum Command {
    /// Trajectory analyses — one stacked CSV per run (gr, sq, msd, angle, vacf, ...)
    Traj {
        #[command(subcommand)]
        cmd: Option<cmd::traj::TrajCmd>,
    },
    /// Spatial maps — one .cube grid file per input (density, sdf, ...)
    Map {
        #[command(subcommand)]
        cmd: Option<cmd::map::MapCmd>,
    },
    /// Glass network topology — Qn speciation, ligand types, coordination numbers
    Net(cmd::net::NetCmd),
    /// Bader charge partitioning from CHGCAR or cube
    Bader(cmd::bader::BaderCmd),
    /// Convert between structure/trajectory formats
    Convert(cmd::convert::ConvertCmd),
    /// Print structure / trajectory information
    Info(cmd::info::InfoCmd),
    /// Generate quantum-chemistry input files
    Job(Box<cmd::job::JobCmd>),
    /// Build machine-learning training sets from AIMD output
    Dataset {
        #[command(subcommand)]
        cmd: Option<cmd::dataset::DatasetCmd>,
    },
    /// Read the user manual (`ferro doc` lists the topics)
    Doc {
        /// Topic, named after the subcommand: `dataset filter`, `traj gr`, `net`
        #[arg(value_name = "TOPIC", num_args = 0..)]
        topic: Vec<String>,
    },
}

fn main() -> Result<()> {
    // `--P-O=2.3` 这类 cutoff 的元素对写在**参数名**里，clap 无法建模，
    // 所以在解析前先从 argv 剥离，剩下的交给 clap（只有 `net` 用得到）
    let argv: Vec<String> = std::env::args().collect();
    let (pair_args, clap_argv) = cmd::net::split_pair_args(&argv);

    let cli = Cli::parse_from(clap_argv);

    let Some(command) = cli.command else {
        help::print_overview();
        return Ok(());
    };

    let failures = match &command {
        Command::Traj { cmd: None } => {
            help::print_traj_overview();
            0
        }
        Command::Traj { cmd: Some(c) } => {
            if cmd::traj::wants_help(c) {
                cmd::traj::print_help(c);
                0
            } else {
                cmd::traj::run(c)?
            }
        }
        Command::Map { cmd: None } => {
            help::print_map_overview();
            0
        }
        Command::Map { cmd: Some(c) } => {
            if cmd::map::wants_help(c) {
                cmd::map::print_help(c);
                0
            } else {
                cmd::map::run(c)?
            }
        }
        Command::Net(c) => {
            if cmd::net::wants_help(c) {
                cmd::net::print_help();
                0
            } else {
                cmd::net::run(c, &pair_args)?
            }
        }
        // 裸命令（无 -i）打富文本帮助，与 traj / map / net 一致；`-h` 仍归 clap 的参数表
        Command::Bader(c) => {
            if cmd::bader::wants_help(c) { help::print_bader() } else { cmd::bader::run(c)? }
            0
        }
        Command::Convert(c) => {
            if cmd::convert::wants_help(c) { help::print_convert() } else { cmd::convert::run(c)? }
            0
        }
        Command::Info(c) => {
            if cmd::info::wants_help(c) { help::print_info() } else { cmd::info::run(c)? }
            0
        }
        Command::Job(c)     => { cmd::job::run(c)?; 0 }
        Command::Dataset { cmd: None } => {
            help::print_dataset_overview();
            0
        }
        Command::Doc { topic } => {
            doc::run(topic)?;
            0
        }
        Command::Dataset { cmd: Some(c) } => {
            if cmd::dataset::wants_help(c) {
                cmd::dataset::print_help(c);
                0
            } else {
                cmd::dataset::run(c)?
            }
        }
    };

    if failures > 0 {
        // 退出码非零：否则 shell 里 `ferro traj gr ... && next-step` 会把批内失败当成功。
        // 不承诺失败原因在哪个文件里：`[inputs]` 块只有分析命令的 csv 才有，
        // 而 dataset 的产物是目录。原因每个 SKIP 行都已经说过了，这句的
        // 唯一职责是解释退出码 —— 指错地方比不指路更糟
        eprintln!("{failures} input(s) failed; see the messages above");
        std::process::exit(1);
    }
    Ok(())
}
