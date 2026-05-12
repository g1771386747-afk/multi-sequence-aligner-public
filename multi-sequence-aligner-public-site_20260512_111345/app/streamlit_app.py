from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime
from pathlib import Path

import streamlit as st


SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "multi_sequence_align.py"
APP_DATA_ROOT = Path(os.environ.get("MULTISEQ_APP_DATA", Path(tempfile.gettempdir()) / "multi_sequence_aligner"))
DEFAULT_NUMBER_MAP = "1=11,2=51,3=1,4=66,5=13,6=15,7=64,8=41,9=58"


st.set_page_config(
    page_title="多序列比对分析平台",
    page_icon="DNA",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.3rem; padding-bottom: 2.5rem; max-width: 1440px;}
    [data-testid="stSidebar"] {border-right: 1px solid #e6e9ef;}
    h1 {font-size: 2rem; margin-bottom: .25rem; letter-spacing: 0;}
    h2, h3 {letter-spacing: 0;}
    h2 {font-size: 1.2rem;}
    h3 {font-size: 1rem;}
    div[data-testid="stMetric"] {
        background: #f8fafc;
        border: 1px solid #e6e9ef;
        border-radius: 8px;
        padding: .78rem .88rem;
    }
    div[data-testid="stMetric"] label {font-size: .85rem; color: #526071;}
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {font-size: 1.35rem;}
    .status-ok {color: #087443; font-weight: 600;}
    .status-warn {color: #9a5b00; font-weight: 600;}
    .status-bad {color: #b42318; font-weight: 600;}
    .muted {color: #5f6b7a; font-size: .92rem;}
    </style>
    """,
    unsafe_allow_html=True,
)


def session_id() -> str:
    if "session_id" not in st.session_state:
        st.session_state.session_id = uuid.uuid4().hex
    return st.session_state.session_id


def session_root() -> Path:
    root = APP_DATA_ROOT / "sessions" / session_id()
    root.mkdir(parents=True, exist_ok=True)
    return root


def safe_upload_name(name: str) -> str:
    base = Path(name).name
    base = re.sub(r'[\\/:*?"<>|]+', "_", base)
    return base.strip(" ._") or "uploaded_file"


def uploaded_bytes(uploaded) -> bytes:
    return bytes(uploaded.getbuffer())


def file_size_label(size: int) -> str:
    if size >= 1024 * 1024:
        return f"{size / 1024 / 1024:.1f} MB"
    if size >= 1024:
        return f"{size / 1024:.1f} KB"
    return f"{size} B"


def write_uploaded_file(uploaded, folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / safe_upload_name(uploaded.name)
    if path.exists():
        digest = hashlib.sha1(uploaded_bytes(uploaded)).hexdigest()[:8]
        path = folder / f"{path.stem}_{digest}{path.suffix}"
    path.write_bytes(uploaded_bytes(uploaded))
    return path


def save_inputs(primer_upload, fasta_uploads, sanger_uploads) -> tuple[Path, list[Path], Path | None, Path]:
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_root = session_root() / "runs" / run_id
    input_root = run_root / "inputs"
    primer_path = write_uploaded_file(primer_upload, input_root / "primer")
    fasta_paths = [write_uploaded_file(item, input_root / "genomes") for item in fasta_uploads]
    sanger_dir: Path | None = None
    if sanger_uploads:
        sanger_dir = input_root / "sanger"
        for item in sanger_uploads:
            write_uploaded_file(item, sanger_dir)
    return primer_path, fasta_paths, sanger_dir, run_root


def cache_primer_for_preview(uploaded) -> Path | None:
    if uploaded is None:
        return None
    data = uploaded_bytes(uploaded)
    digest = hashlib.sha256(data).hexdigest()[:16]
    folder = session_root() / "preview"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{digest}_{safe_upload_name(uploaded.name)}"
    if not path.exists():
        path.write_bytes(data)
    return path


def read_primer_preview(path: Path | None) -> tuple[list[dict[str, str]], str | None]:
    if path is None or not path.exists():
        return [], "请上传引物表 Excel"
    try:
        sys.path.insert(0, str(SCRIPT.parent))
        from multi_sequence_align import read_primers  # type: ignore

        primers = read_primers(path)
        rows = [
            {
                "序号": str(p.order),
                "引物对": p.pair_name,
                "F名称": p.f_name,
                "F序列": p.f_seq,
                "R名称": p.r_name,
                "R序列": p.r_seq,
            }
            for p in primers[:100]
        ]
        return rows, None
    except Exception as exc:  # pragma: no cover - displayed in UI
        return [], f"无法读取引物表：{exc}"


def build_command(
    *,
    primer_table: Path,
    fasta_paths: list[Path],
    sanger_dir: Path | None,
    out_dir: Path,
    number_map: str,
    do_stage1: bool,
    do_align: bool,
    do_sanger: bool,
    do_word: bool,
    zip_pdfs: bool,
    flank: int,
    max_product: int,
    align_samples: str,
    colors: dict[str, str],
) -> list[str]:
    if do_sanger and do_align:
        mode = "full"
    elif do_sanger:
        mode = "sanger-compare"
    elif do_align:
        mode = "full"
    elif do_stage1:
        mode = "stage1"
    else:
        return []

    cmd = [
        sys.executable,
        str(SCRIPT),
        mode,
        "--primer-table",
        str(primer_table),
        "--out-dir",
        str(out_dir),
        "--flank",
        str(flank),
        "--max-product",
        str(max_product),
    ]
    for fasta in fasta_paths:
        cmd += ["--fasta", str(fasta)]

    if do_align or mode == "full":
        cmd += [
            "--align-samples",
            align_samples.strip() or "all",
            "--color-a",
            colors["A"],
            "--color-t",
            colors["T"],
            "--color-c",
            colors["C"],
            "--color-g",
            colors["G"],
            "--color-gap",
            colors["gap"],
        ]
    if do_sanger:
        cmd += ["--sanger-dir", str(sanger_dir or ""), "--number-map", number_map]
        if do_word:
            cmd += ["--word-report"]
    if zip_pdfs:
        cmd += ["--zip-pdfs"]
    return cmd


def run_command(args: list[str]) -> tuple[int, str]:
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )
    return proc.returncode, proc.stdout


def make_result_zip(result_dir: Path) -> Path:
    zip_base = result_dir.with_name(result_dir.name + "_完整结果")
    return Path(shutil.make_archive(str(zip_base), "zip", result_dir))


def clear_session_files():
    root = session_root()
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)


def cleanup_old_sessions(max_age_hours: int = 24):
    sessions_root = APP_DATA_ROOT / "sessions"
    if not sessions_root.exists():
        return
    cutoff = time.time() - max_age_hours * 3600
    for item in sessions_root.iterdir():
        try:
            if item.is_dir() and item.stat().st_mtime < cutoff:
                shutil.rmtree(item, ignore_errors=True)
        except OSError:
            continue


cleanup_old_sessions()
st.title("多序列比对分析平台")
st.caption("上传引物表和样本 FASTA，在线生成理论扩增片段、多序列比对结果和可下载报告。")

with st.sidebar:
    st.header("上传数据")
    primer_upload = st.file_uploader("引物表 Excel", type=["xlsx"])
    fasta_uploads = st.file_uploader(
        "样本 FASTA 文件",
        type=["fa", "fas", "fasta", "fna"],
        accept_multiple_files=True,
    )

    st.header("分析内容")
    do_stage1 = st.checkbox("理论扩增片段 PDF", value=True)
    do_align = st.checkbox("多样本序列比对 PDF", value=True)
    do_sanger = st.checkbox("胶回收测序 vs 理论序列 PDF", value=False)
    do_word = st.checkbox("生成 Word 汇总报告", value=False, disabled=not do_sanger)
    zip_pdfs = st.checkbox("生成 PDF 压缩包", value=True)

    sanger_uploads = st.file_uploader(
        "胶回收测序序列文件",
        type=["fa", "fas", "fasta", "seq", "txt"],
        accept_multiple_files=True,
        disabled=not do_sanger,
    )

    st.header("参数")
    align_samples = st.text_input("参与比对的样本", "all", help="填写 all，或用英文逗号分隔 FASTA 文件名对应的样本名。")
    number_map = st.text_input("测序编号对应样本编号", DEFAULT_NUMBER_MAP, disabled=not do_sanger)
    flank = st.number_input("两侧延伸长度 bp", min_value=0, max_value=5000, value=500, step=50)
    max_product = st.number_input("最大理论扩增长度 bp", min_value=100, max_value=50000, value=5000, step=100)

    st.header("碱基颜色")
    color_left, color_right = st.columns(2)
    with color_left:
        color_a = st.color_picker("A", "#1aa260")
        color_c = st.color_picker("C", "#2568c9")
        color_gap = st.color_picker("gap", "#d9d9d9")
    with color_right:
        color_t = st.color_picker("T", "#d83b36")
        color_g = st.color_picker("G", "#111111")

    run = st.button("开始分析", type="primary", use_container_width=True)
    clear = st.button("清除本次上传与结果", use_container_width=True)

fasta_uploads = fasta_uploads or []
sanger_uploads = sanger_uploads or []

if clear:
    clear_session_files()
    st.session_state.pop("last_zip", None)
    st.session_state.pop("last_result_dir", None)
    st.rerun()

primer_preview_path = cache_primer_for_preview(primer_upload)
primer_rows, primer_error = read_primer_preview(primer_preview_path)
fasta_preview = [
    {"样本": Path(item.name).stem, "文件名": item.name, "大小": file_size_label(item.size)}
    for item in fasta_uploads
]
sanger_preview = [
    {"文件名": item.name, "大小": file_size_label(item.size)}
    for item in sanger_uploads
]

metric_cols = st.columns(4)
metric_cols[0].metric("引物对", len(primer_rows) if not primer_error else 0)
metric_cols[1].metric("FASTA 文件", len(fasta_preview))
metric_cols[2].metric("测序文件", len(sanger_preview))
metric_cols[3].metric("运行模式", "完整分析" if do_sanger and do_align else "测序比对" if do_sanger else "多序列比对" if do_align else "理论扩增")

check_tab, primer_tab, fasta_tab, sanger_tab, result_tab = st.tabs(
    ["输入检查", "引物预览", "FASTA 预览", "测序文件", "运行与下载"]
)

with check_tab:
    st.subheader("输入检查")
    checks = [
        ("引物表", primer_upload is not None, "必须上传 Excel 引物表。"),
        ("FASTA 文件", bool(fasta_uploads), "至少上传一个样本 FASTA。"),
        ("测序文件", bool(sanger_uploads), "仅在勾选测序比对时需要上传。"),
    ]
    for name, ok, note in checks:
        optional = name == "测序文件" and not do_sanger
        status = "正常" if ok else ("可选缺失" if optional else "缺失")
        css = "status-ok" if ok else ("status-warn" if optional else "status-bad")
        st.markdown(f"**{name}**：<span class='{css}'>{status}</span>", unsafe_allow_html=True)
        st.caption(note)
    st.markdown(
        "<div class='muted'>每次运行的数据会保存在服务器临时目录中，结果通过压缩包下载。</div>",
        unsafe_allow_html=True,
    )

with primer_tab:
    st.subheader("引物表预览")
    if primer_error:
        st.info(primer_error)
    else:
        st.dataframe(primer_rows, use_container_width=True, hide_index=True)
        if len(primer_rows) >= 100:
            st.caption("这里只显示前 100 条，引物会按表格原始顺序参与分析。")

with fasta_tab:
    st.subheader("样本 FASTA 预览")
    if fasta_preview:
        st.dataframe(fasta_preview, use_container_width=True, hide_index=True)
    else:
        st.info("请上传样本 FASTA 文件。")

with sanger_tab:
    st.subheader("胶回收测序文件预览")
    if sanger_preview:
        st.dataframe(sanger_preview, use_container_width=True, hide_index=True)
    else:
        st.info("未上传测序文件；未勾选测序比对时可以忽略。")

with result_tab:
    st.subheader("运行与下载")
    if st.session_state.get("last_zip"):
        zip_path = Path(st.session_state.last_zip)
        if zip_path.exists():
            st.download_button(
                "下载上一次完整结果压缩包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )

    if not run:
        st.info("确认输入无误后，点击左侧“开始分析”。")

    if run:
        if not (do_stage1 or do_align or do_sanger):
            st.error("请至少选择一个分析内容。")
            st.stop()
        if primer_upload is None:
            st.error("请上传引物表 Excel。")
            st.stop()
        if not fasta_uploads:
            st.error("请至少上传一个样本 FASTA 文件。")
            st.stop()
        if do_sanger and not sanger_uploads:
            st.error("已选择测序比对，请上传胶回收测序序列文件。")
            st.stop()

        primer_path, fasta_paths, sanger_dir, run_root = save_inputs(primer_upload, fasta_uploads, sanger_uploads)
        result_dir = run_root / "result"
        colors = {"A": color_a, "T": color_t, "C": color_c, "G": color_g, "gap": color_gap}
        command = build_command(
            primer_table=primer_path,
            fasta_paths=fasta_paths,
            sanger_dir=sanger_dir,
            out_dir=result_dir,
            number_map=number_map,
            do_stage1=do_stage1,
            do_align=do_align,
            do_sanger=do_sanger,
            do_word=do_word,
            zip_pdfs=zip_pdfs,
            flank=int(flank),
            max_product=int(max_product),
            align_samples=align_samples,
            colors=colors,
        )

        progress = st.progress(10, text="正在准备分析...")
        with st.spinner("正在分析，请等待。大型数据可能需要几分钟。"):
            progress.progress(35, text="正在运行分析程序...")
            code, output = run_command(command)
            progress.progress(100, text="分析结束")

        st.text_area("运行日志", output, height=260)
        if code == 0:
            zip_path = make_result_zip(result_dir)
            st.session_state.last_result_dir = str(result_dir)
            st.session_state.last_zip = str(zip_path)
            st.success("分析完成")
            st.download_button(
                "下载完整结果压缩包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )
        else:
            st.error("分析失败，请根据日志检查输入文件或参数。")
