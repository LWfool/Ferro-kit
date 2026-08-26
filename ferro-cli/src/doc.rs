//! `ferro doc` — the user manual, compiled into the binary.
//!
//! The manual lives in `docs/src/` as mdBook sources and is embedded here with
//! `include_str!`. Reading it from disk would work only when the source tree is
//! at hand, which is exactly the situation a `cargo install`ed binary is not in.
//! At 208 KB it is a rounding error next to the basis-set tables already inside.
//!
//! Topics are named after the subcommand tree (`ferro doc dataset filter`), so
//! the pointer a help page prints is the next command to type rather than a path
//! to go looking for. Pages that describe no single command keep a flat name.
//!
//! The markdown is printed as written. These sources are meant to be read by
//! people — headings are `##`, code is fenced — so rendering them buys little,
//! and every renderer worth having is a dependency.

use std::io::{IsTerminal, Write};
use std::process::{Command, Stdio};

use anyhow::Result;

/// One manual page: how it is addressed, and what it contains.
struct Page {
    /// What the user types after `ferro doc`.
    topic: &'static str,
    /// The mdBook source path, shown when listing.
    source: &'static str,
    text: &'static str,
}

/// Every page of `docs/src/`, in `SUMMARY.md` order.
///
/// The order is the manual's own, not alphabetical: `ferro doc` with no topic
/// prints this list, and a reader scanning it should meet the pages in the order
/// the book introduces them.
const PAGES: &[Page] = &[
    Page {
        topic: "introduction",
        source: "introduction.md",
        text: include_str!("../../docs/src/introduction.md"),
    },
    Page {
        topic: "installation",
        source: "installation.md",
        text: include_str!("../../docs/src/installation.md"),
    },
    Page {
        topic: "data-model",
        source: "data-model.md",
        text: include_str!("../../docs/src/data-model.md"),
    },
    Page {
        topic: "traj gr",
        source: "analysis/gr.md",
        text: include_str!("../../docs/src/analysis/gr.md"),
    },
    Page {
        topic: "traj sq",
        source: "analysis/sq.md",
        text: include_str!("../../docs/src/analysis/sq.md"),
    },
    Page {
        topic: "traj msd",
        source: "analysis/msd.md",
        text: include_str!("../../docs/src/analysis/msd.md"),
    },
    Page {
        topic: "traj angle",
        source: "analysis/angle.md",
        text: include_str!("../../docs/src/analysis/angle.md"),
    },
    Page {
        topic: "traj vanhove",
        source: "analysis/vanhove.md",
        text: include_str!("../../docs/src/analysis/vanhove.md"),
    },
    Page {
        topic: "traj vacf",
        source: "analysis/vacf.md",
        text: include_str!("../../docs/src/analysis/vacf.md"),
    },
    Page {
        topic: "traj rotcorr",
        source: "analysis/rotcorr.md",
        text: include_str!("../../docs/src/analysis/rotcorr.md"),
    },
    Page {
        topic: "map density",
        source: "analysis/cube-density.md",
        text: include_str!("../../docs/src/analysis/cube-density.md"),
    },
    Page {
        topic: "map radius",
        source: "analysis/cube-radius.md",
        text: include_str!("../../docs/src/analysis/cube-radius.md"),
    },
    Page {
        topic: "map sdf",
        source: "analysis/cube-sdf.md",
        text: include_str!("../../docs/src/analysis/cube-sdf.md"),
    },
    Page {
        topic: "map chg-sdf",
        source: "analysis/chg-sdf.md",
        text: include_str!("../../docs/src/analysis/chg-sdf.md"),
    },
    Page {
        topic: "cube-jump",
        source: "analysis/cube-jump.md",
        text: include_str!("../../docs/src/analysis/cube-jump.md"),
    },
    Page {
        topic: "net",
        source: "analysis/network.md",
        text: include_str!("../../docs/src/analysis/network.md"),
    },
    Page {
        topic: "dataset collect",
        source: "dataset/collect.md",
        text: include_str!("../../docs/src/dataset/collect.md"),
    },
    Page {
        topic: "dataset filter",
        source: "dataset/filter.md",
        text: include_str!("../../docs/src/dataset/filter.md"),
    },
    Page {
        topic: "dataset merge",
        source: "dataset/merge.md",
        text: include_str!("../../docs/src/dataset/merge.md"),
    },
    Page {
        topic: "job",
        source: "workflow/job-builders.md",
        text: include_str!("../../docs/src/workflow/job-builders.md"),
    },
    Page {
        topic: "spin",
        source: "workflow/spin.md",
        text: include_str!("../../docs/src/workflow/spin.md"),
    },
    Page {
        topic: "python",
        source: "python.md",
        text: include_str!("../../docs/src/python.md"),
    },
    Page {
        topic: "cli-reference",
        source: "cli-reference.md",
        text: include_str!("../../docs/src/cli-reference.md"),
    },
];

/// Topics that are one command's page but reachable under a shorter name too.
///
/// `map velocity` and `map force` share the density page because the three modes
/// differ only in what they accumulate into the same grid.
const ALIASES: &[(&str, &str)] = &[
    ("map velocity", "map density"),
    ("map force", "map density"),
    ("network", "net"),
    ("gr", "traj gr"),
    ("sq", "traj sq"),
    ("msd", "traj msd"),
    ("angle", "traj angle"),
    ("collect", "dataset collect"),
    ("filter", "dataset filter"),
    ("merge", "dataset merge"),
];

/// Prints one topic, or the list of topics when none is given.
pub fn run(topic: &[String]) -> Result<()> {
    if topic.is_empty() {
        print!("{}", topic_list());
        return Ok(());
    }
    let asked = topic.join(" ").to_lowercase();
    let resolved = ALIASES
        .iter()
        .find(|(from, _)| *from == asked)
        .map(|(_, to)| *to)
        .unwrap_or(&asked);

    match PAGES.iter().find(|p| p.topic == resolved) {
        Some(page) => page_out(page.text),
        None => {
            // 猜错的主题名不该只说「没有」——列表就在手边,直接给出来
            eprintln!("no manual page for `{asked}`\n");
            eprint!("{}", topic_list());
            std::process::exit(1);
        }
    }
    Ok(())
}

fn topic_list() -> String {
    let width = PAGES.iter().map(|p| p.topic.len()).max().unwrap_or(0);
    let mut s = String::from("ferro doc — the user manual\n\nTopics:\n");
    for p in PAGES {
        s.push_str(&format!("  {:<width$}  {}\n", p.topic, p.source));
    }
    s.push_str("\nRead one with:\n  ferro doc dataset filter\n");
    s
}

/// Writes the page out, through `$PAGER` when there is a terminal to page for.
///
/// `git` and `man` behave this way: interactive output is paged, redirected
/// output is not, so `ferro doc net > net.md` stays a plain file. A missing or
/// unusable pager falls back to printing rather than failing — the page is what
/// was asked for; how it scrolls is not.
fn page_out(text: &str) {
    if !std::io::stdout().is_terminal() {
        print!("{text}");
        return;
    }
    let pager = std::env::var("PAGER").unwrap_or_else(|_| "less -R".to_string());
    let mut parts = pager.split_whitespace();
    let Some(program) = parts.next() else {
        print!("{text}");
        return;
    };

    let spawned = Command::new(program)
        .args(parts)
        .stdin(Stdio::piped())
        .spawn();
    let Ok(mut child) = spawned else {
        print!("{text}");
        return;
    };
    if let Some(stdin) = child.stdin.as_mut() {
        // 用户在读到末尾前按 q,写入会收到 EPIPE —— 那是正常退出,不是错误
        let _ = stdin.write_all(text.as_bytes());
    }
    let _ = child.wait();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_topic_is_unique() {
        let mut seen: Vec<&str> = PAGES.iter().map(|p| p.topic).collect();
        seen.sort_unstable();
        let n = seen.len();
        seen.dedup();
        assert_eq!(seen.len(), n, "two pages share a topic name");
    }

    #[test]
    fn every_alias_points_at_a_real_topic() {
        for (from, to) in ALIASES {
            assert!(
                PAGES.iter().any(|p| p.topic == *to),
                "alias `{from}` points at `{to}`, which is not a topic"
            );
            assert!(
                !PAGES.iter().any(|p| p.topic == *from),
                "alias `{from}` shadows a real topic"
            );
        }
    }

    #[test]
    fn no_page_is_empty() {
        for p in PAGES {
            assert!(p.text.len() > 200, "{} looks truncated", p.source);
        }
    }
}
