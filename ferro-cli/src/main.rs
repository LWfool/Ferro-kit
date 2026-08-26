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

/// The clap definitions and the rich help pages state the same facts twice.
///
/// `--shuffle` was added to `dataset filter` and its page's parameter table was
/// not updated; the drift went unnoticed until a user looked for the flag. These
/// tests turn that from "someone notices eventually" into a red `cargo test`.
///
/// The pages are extracted from this crate's own source rather than captured
/// from stdout: `help::print_*` writes to the terminal, and refactoring 25
/// functions to return their text would be a larger change than the check is
/// worth.
#[cfg(test)]
mod help_sync {
    use clap::CommandFactory;

    const HELP_SRC: &str = include_str!("help.rs");
    /// `ferro net` keeps its page next to the argv-stripping it needs, not in
    /// `help.rs`. Left out of the scan it would be the one command the check
    /// never sees — and it has 10 options.
    const NET_SRC: &str = include_str!("cmd/net.rs");

    /// Which rich page belongs to which subcommand path.
    ///
    /// A wrong function name here fails loudly (the page is not found), so the
    /// table cannot silently rot the way the two hand-written lists it guards do.
    const PAGES: &[(&[&str], &str)] = &[
        (&["traj", "gr"], "print_gr"),
        (&["traj", "sq"], "print_sq"),
        (&["traj", "msd"], "print_msd"),
        (&["traj", "angle"], "print_angle"),
        (&["traj", "vacf"], "print_vacf"),
        (&["traj", "rotcorr"], "print_rotcorr"),
        (&["traj", "vanhove"], "print_vanhove"),
        (&["map", "density"], "print_cube_density"),
        (&["map", "velocity"], "print_cube_velocity"),
        (&["map", "force"], "print_cube_force"),
        (&["map", "radius"], "print_cube_radius"),
        (&["map", "sdf"], "print_cube_sdf"),
        (&["map", "chg-sdf"], "print_cube_chg_sdf"),
        (&["convert"], "print_convert"),
        (&["info"], "print_info"),
        (&["bader"], "print_bader"),
        (&["net"], "HELP_EXTRA"),
        (&["dataset", "collect"], "print_dataset_collect"),
        (&["dataset", "filter"], "print_dataset_filter"),
        (&["dataset", "merge"], "print_dataset_merge"),
    ];

    /// Options clap accepts on a command that its page deliberately omits.
    ///
    /// `SelectArgs` is shared by `gr` and `angle`, so `gr` accepts `-c`/`-z` even
    /// though a radial distribution has no third atom. Removing the *use* of a
    /// shared argument group does not remove the argument, and documenting one
    /// that means nothing here would be worse than the silence.
    const UNDOCUMENTED: &[(&[&str], &str)] =
        &[(&["traj", "gr"], "atom-c"), (&["traj", "gr"], "label-z")];

    /// The body of `fn <name>(` up to the closing brace in column 0, or — for
    /// `ferro net` — the `HELP_EXTRA` string constant in `cmd/net.rs`.
    fn page(name: &str) -> &'static str {
        if name == "HELP_EXTRA" {
            let at = NET_SRC.find("const HELP_EXTRA").expect("cmd/net.rs has no HELP_EXTRA");
            let rest = &NET_SRC[at..];
            let end = rest.find("\";").map(|i| i + 2).unwrap_or(rest.len());
            return &rest[..end];
        }
        let needle = format!("fn {name}() {{");
        let at = HELP_SRC
            .find(&needle)
            .unwrap_or_else(|| panic!("help.rs has no `{name}`; the PAGES table is stale"));
        let rest = &HELP_SRC[at..];
        let end = rest.find("\n}\n").unwrap_or(rest.len());
        &rest[..end]
    }

    fn subcommand(path: &[&str]) -> clap::Command {
        let mut cmd = super::Cli::command();
        for step in path {
            let next = cmd
                .get_subcommands()
                .find(|c| c.get_name() == *step)
                .unwrap_or_else(|| panic!("no subcommand `{}`", path.join(" ")))
                .clone();
            cmd = next;
        }
        cmd
    }

    #[test]
    fn every_option_appears_on_its_page() {
        for (path, fname) in PAGES {
            let text = page(fname);
            for arg in subcommand(path).get_arguments() {
                let Some(long) = arg.get_long() else { continue };
                if long == "help" || long == "version" {
                    continue;
                }
                if UNDOCUMENTED.iter().any(|(p, o)| p == path && *o == long) {
                    continue;
                }
                // 页里写短名（`-a ELEM`、`-o SUFFIX`）也算写了 —— 要断言的是
                // 这个参数被交代过,不是它以哪种拼法出现
                let named = text.contains(&format!("--{long}"))
                    || arg.get_short().is_some_and(|c| text.contains(&format!("-{c} ")));
                assert!(
                    named,
                    "`ferro {}` takes --{long} but {fname} never mentions it; \
                     if that is deliberate, add it to UNDOCUMENTED with the reason",
                    path.join(" ")
                );
            }
        }
    }

    #[test]
    fn every_option_named_on_a_page_exists() {
        // 反向的一半：删掉参数忘了改帮助，帮助会教人敲一条会报错的命令 ——
        // 比漏写更坏，而正向断言抓不到
        // 允许交叉引用别的命令的参数(`ferro net --export-traj` 就出现在 gr 页里),
        // 只要它在 ferro 的某个命令上真的存在。要抓的是**哪儿都不存在**的那种
        let mut known: Vec<String> = Vec::new();
        collect_longs(&super::Cli::command(), &mut known);

        for (_, fname) in PAGES {
            for line in page(fname).lines() {
                // 「there is no --from / --to flag」这类**否定陈述**正是在做对的事:
                // 明说某个拼法不存在。断言它存在会把这句话判成错的
                if line.contains("no --") {
                    continue;
                }
                for name in longs_in(line) {
                    if name == "help" || name == "version" {
                        continue;
                    }
                    assert!(
                        known.contains(&name),
                        "{fname} names --{name}, which no ferro command takes"
                    );
                }
            }
        }
    }

    /// Every `--name` written in one line of text.
    ///
    /// Scans for the `--` marker rather than splitting on whitespace, so a
    /// slash-separated list (`--start/--end/--stride`) yields all three.
    fn longs_in(line: &str) -> Vec<String> {
        let bytes: Vec<char> = line.chars().collect();
        let mut out = Vec::new();
        let mut i = 0;
        while i + 2 < bytes.len() {
            if bytes[i] == '-' && bytes[i + 1] == '-' {
                let mut j = i + 2;
                while j < bytes.len() && (bytes[j].is_ascii_alphanumeric() || bytes[j] == '-') {
                    j += 1;
                }
                let name: String =
                    bytes[i + 2..j].iter().collect::<String>().trim_matches('-').to_string();
                // `--P-O=2.4` 是 net 的配对参数：元素对写在**参数名**里,clap 建模
                // 不了,所以 main 在解析前就从 argv 剥离了。clap 永远不知道它,
                // 断言它存在必然失败。判据取那个 `=` —— 别处的参数从不这么写
                let assigned = j < bytes.len() && bytes[j] == '=';
                if !name.is_empty() && !assigned {
                    out.push(name);
                }
                i = j.max(i + 2);
            } else {
                i += 1;
            }
        }
        out
    }

    fn collect_longs(cmd: &clap::Command, out: &mut Vec<String>) {
        out.extend(cmd.get_arguments().filter_map(|a| a.get_long()).map(String::from));
        for sub in cmd.get_subcommands() {
            collect_longs(sub, out);
        }
    }
}
