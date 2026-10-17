#!/bin/bash
# QE pw.x 输出 fixture 的生成脚本（2026-10-17）：QE 单点 / AIMD reader 的测试样例全部由它在超算上算出，
# 每个算例都完整跑完（不像 smoke_hpc.sh 只看输入解析）。小体系几分钟内跑完。
#
# 用法:  bash qe_fixtures_hpc.sh STRUCTURE WORKDIR
#   STRUCTURE  小的周期结构（10 个原子以内，如冒烟用的 6 原子 AlPO4），ferro 能读的任何格式
#   WORKDIR    放全部算例的目录（不存在会新建）
# 先决条件: ferro-py 在 PATH 上，配置里有 [qe] pseudo_dir（同 smoke_hpc.sh）
# 环境变量:  PW_CMD  默认 "mpirun -np 4 pw.x"
#
# 算例（目录名 = 用途）:
#   sp          单点，D3，verbosity high：能量、总力（后面跟着各项力分量，不能读成总力）、应力
#   sp_spin     同上，--multiplicity 3（nspin = 2）：自旋输出插在中间时解析不乱
#   sp_noconv   electron_maxstep = 3：SCF 不收敛
#   md          NVT（SVR）5 步，1000 K：逐步的能量、力、应力、坐标
#   vcmd        vc-md 5 步（beeman + pr，press / wmass 取默认，无热浴）：逐步变化的晶胞
# 结束后打包 WORKDIR/qe_fixtures.tar.gz（只含 pw.in、pw.out、run.log），拷回即可。
set -u

STRUCTURE=${1:?usage: bash qe_fixtures_hpc.sh STRUCTURE WORKDIR}
WORK=${2:?usage: bash qe_fixtures_hpc.sh STRUCTURE WORKDIR}
PW_CMD=${PW_CMD:-"mpirun -np 4 pw.x"}

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
STRUCTURE=$(cd "$(dirname "$STRUCTURE")" && pwd)/$(basename "$STRUCTURE")
STEM=$(basename "$STRUCTURE"); STEM=${STEM%.*}

# 生成一个算例：$1 目录名  $2 任务  其余参数原样交给 ferro-py job；成功时打印 pw.in 所在目录
gen() {
    local name=$1 task=$2; shift 2
    rm -rf "$WORK/$name"
    # stderr 单独存：stdout 有缓冲、退出时才落盘，合在一起时末行是 Created: 而不是报错
    if ! ferro-py job -s qe -t "$task" -i "$STRUCTURE" -o "$WORK/$name" --mkdir "$@" \
            > "$WORK/$name.gen.log" 2> "$WORK/$name.gen.err"; then
        echo "$name: GENERATION FAILED  $(tail -1 "$WORK/$name.gen.err")" >&2
        return 1
    fi
    echo "$WORK/$name/$STEM"
}

# 按正则改 pw.in；改不到就报错（模板变了要知道）
edit() {   # $1 pw.in  $2 python 代码，对变量 s 做替换
    python - "$1" "$2" <<'EOF'
import re, sys
path, code = sys.argv[1], sys.argv[2]
s = open(path).read()
before = s
exec(code)
if s == before:
    sys.exit(f"edit made no change in {path}: {code}")
open(path, "w").write(s)
EOF
}

run() {    # $1 目录
    ( cd "$1" && bash -c "$PW_CMD -in pw.in > pw.out" > run.log 2>&1 )
    local rc=$?
    if grep -q "JOB DONE" "$1/pw.out" 2>/dev/null; then echo "done (rc=$rc)"; else echo "NOT DONE (rc=$rc), see $1/run.log"; fi
}

dirs=()

d=$(gen sp sp) && { echo "sp: $(run "$d")"; dirs+=("$d"); }

d=$(gen sp_spin sp --multiplicity 3) && { echo "sp_spin: $(run "$d")"; dirs+=("$d"); }

if d=$(gen sp_noconv sp); then
    edit "$d/pw.in" 's = re.sub(r"(?i)electron_maxstep\s*=\s*\d+", "electron_maxstep = 3", s)' \
        && { echo "sp_noconv: $(run "$d") (expected: convergence NOT achieved)"; dirs+=("$d"); }
fi

d=$(gen md aimd --steps 5 --temperature 1000) && { echo "md: $(run "$d")"; dirs+=("$d"); }

if d=$(gen vcmd aimd --steps 5); then
    # vc-md：INPUT_PW 只支持 beeman；热浴不支持，去掉 ion_temperature / tempw / nraise；加 &CELL
    edit "$d/pw.in" '
s = re.sub(r"(?i)calculation\s*=\s*.md.", "calculation = '"'"'vc-md'"'"'", s)
s = re.sub(r"(?i)ion_dynamics\s*=\s*.verlet.", "ion_dynamics = '"'"'beeman'"'"'", s)
s = re.sub(r"(?im)^\s*(ion_temperature|tempw|nraise)\s*=.*\n", "", s)
s = re.sub(r"(?is)(&IONS.*?\n\s*/\s*\n)", r"\1&CELL\n  cell_dynamics = '"'"'pr'"'"'\n/\n", s, count=1)
' && { echo "vcmd: $(run "$d")"; dirs+=("$d"); }
fi

[ ${#dirs[@]} -eq 0 ] && { echo "nothing to pack: no case was generated" >&2; exit 1; }
cd "$WORK" && tar czf qe_fixtures.tar.gz $(for d in "${dirs[@]}"; do
    r=${d#"$WORK"/}; for f in pw.in pw.out run.log; do [ -f "$d/$f" ] && echo "$r/$f"; done
done)
echo "packed: $WORK/qe_fixtures.tar.gz"
