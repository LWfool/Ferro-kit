#!/bin/bash
# ferro-py 超算冒烟测试（ferro.jobs 第 20 步）：VASP / QE / CP2K × sp / opt / aimd 九个算例，
# 每个只跑到「输入解析通过、第一个 SCF 开始」就算通过；再生成收敛测试，可选跑完并 analyze。
#
# 用法:  bash smoke_hpc.sh STRUCTURE WORKDIR
#   STRUCTURE  一个小的周期结构（几十个原子以内，跑得快），ferro 能读的任何格式
#   WORKDIR    放全部算例的目录（不存在会新建）
#
# 先决条件: ferro-py 在 PATH 上；FERRO_CONFIG（或 ~/.config/ferro/config.toml）里配好了
# [vasp] potcar_dir、[qe] pseudo_dir、[cp2k] data_dir。
#
# 按超算修改的环境变量:
#   VASP_CMD       默认 "mpirun -np 4 vasp_std"（Γ 点网格可换 vasp_gam）
#   PW_CMD         默认 "mpirun -np 4 pw.x"
#   CP2K_CMD       默认 "mpirun -np 4 cp2k.psmp"
#   SMOKE_SECONDS  每个算例最多跑几秒，默认 300（看到 SCF 开始即可，超时被杀也算通过）
#   SKIP           跳过的程序，逗号分隔，如 SKIP=qe
#   RUN_CONV=1     收敛测试的单点完整跑完并 analyze（默认 0：只生成）
set -u

STRUCTURE=${1:?usage: bash smoke_hpc.sh STRUCTURE WORKDIR}
WORK=${2:?usage: bash smoke_hpc.sh STRUCTURE WORKDIR}
VASP_CMD=${VASP_CMD:-"mpirun -np 4 vasp_std"}
PW_CMD=${PW_CMD:-"mpirun -np 4 pw.x"}
CP2K_CMD=${CP2K_CMD:-"mpirun -np 4 cp2k.psmp"}
SMOKE_SECONDS=${SMOKE_SECONDS:-300}
SKIP=${SKIP:-}
RUN_CONV=${RUN_CONV:-0}

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
STRUCTURE=$(cd "$(dirname "$STRUCTURE")" && pwd)/$(basename "$STRUCTURE")
STEM=$(basename "$STRUCTURE"); STEM=${STEM%.*}
REPORT=$WORK/smoke_report.txt
: > "$REPORT"

say() { echo "$*" | tee -a "$REPORT"; }
skipped() { [[ ",$SKIP," == *",$1,"* ]]; }

say "ferro-py smoke test  $(date '+%F %T')  host $(hostname)"
say "python: $(python -c 'import sys, ferro; print(sys.version.split()[0], "ferro", ferro.__version__)' 2>&1)"
say "config: ${FERRO_CONFIG:-$HOME/.config/ferro/config.toml}"
say "structure: $STRUCTURE"
say ""

# 跑一个算例：$1 程序  $2 任务目录  $3 命令；按输出里的标记判定
run_case() {
    local prog=$1 dir=$2 cmd=$3 log marker bad
    case $prog in
        vasp) log=$dir/OUTCAR;  marker='POTLOK:';                        bad='ERROR|error reading|Error EDDDAV' ;;
        qe)   log=$dir/pw.out;  marker='Self-consistent Calculation';    bad='Error in routine' ;;
        cp2k) log=$dir/cp2k.out; marker='SCF WAVEFUNCTION OPTIMIZATION'; bad='ABORT' ;;
    esac
    ( cd "$dir" && timeout "$SMOKE_SECONDS" bash -c "$cmd" > run.log 2>&1 )
    local rc=$?
    if [ -f "$log" ] && grep -q "$marker" "$log" && ! grep -Eq "$bad" "$log" "$dir/run.log" 2>/dev/null; then
        if [ $rc -eq 124 ]; then echo "PASS (stopped by the ${SMOKE_SECONDS} s limit)"; else echo "PASS (rc=$rc)"; fi
    else
        echo "FAIL (rc=$rc) see $dir/run.log and $log"
    fi
}

command_for() {
    case $1 in
        vasp) echo "$VASP_CMD" ;;
        qe)   echo "$PW_CMD -in pw.in > pw.out" ;;
        cp2k) echo "$CP2K_CMD -i cp2k.inp > cp2k.out" ;;
    esac
}

say "== 9 jobs"
for prog in vasp qe cp2k; do
    skipped "$prog" && { say "$prog: skipped"; continue; }
    for task in sp opt aimd; do
        extra=""
        [ "$task" = aimd ] && extra="--steps 3"
        out=$WORK/${prog}_${task}
        rm -rf "$out"
        if ! ferro-py job -s "$prog" -t "$task" -i "$STRUCTURE" -o "$out" --mkdir $extra \
                > "$WORK/${prog}_${task}.gen.log" 2>&1; then
            say "$prog $task: GENERATION FAILED  $(tail -1 "$WORK/${prog}_${task}.gen.log")"
            continue
        fi
        say "$prog $task: $(run_case "$prog" "$out/$STEM" "$(command_for "$prog")")"
    done
done

say ""
say "== convergence tests"
conv_case() {   # $1 程序  $2 --scan 参数
    local prog=$1 scan=$2 name=${2%%=*}
    skipped "$prog" && return
    local out=$WORK/conv_${prog}_${name}
    rm -rf "$out"
    if ! ferro-py conv -s "$prog" --scan "$scan" -i "$STRUCTURE" -o "$out" --mkdir > "$out.gen.log" 2>&1; then
        say "$prog conv $scan: GENERATION FAILED  $(tail -1 "$out.gen.log")"
        return
    fi
    if [ "$RUN_CONV" != 1 ]; then
        say "$prog conv $scan: generated ($(( $(ls -d "$out"/*/ | wc -l) )) values)"
        return
    fi
    for d in "$out"/*/"$STEM"; do
        ( cd "$d" && bash -c "$(command_for "$prog")" > run.log 2>&1 )
    done
    say "$prog conv $scan:"
    ferro-py conv analyze "$out" 2>&1 | tee -a "$REPORT"
}
conv_case vasp "fine-grid=2,3"
conv_case vasp "encut=400:600:100"
conv_case cp2k "cutoff=300:500:100"
conv_case cp2k "finer-grid"

say ""
say "Report: $REPORT"
