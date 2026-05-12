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
import streamlit.components.v1 as components


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
    :root {
        --ink: #182230;
        --muted: #667085;
        --line: #e6e9ef;
        --panel: #f8fafc;
        --panel-strong: #eef7f8;
        --accent: #1f7a8c;
        --accent-2: #7a9a01;
        --good: #087443;
        --warn: #9a5b00;
        --bad: #b42318;
    }
    .stApp {
        background: #f6f8fb;
        color: var(--ink);
    }
    .block-container {
        padding-top: 1rem;
        padding-bottom: 2.8rem;
        max-width: 1480px;
    }
    [data-testid="stSidebar"] {
        background: #ffffff;
        border-right: 1px solid var(--line);
    }
    [data-testid="stSidebarContent"] {
        background: #ffffff;
    }
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3 {
        letter-spacing: 0;
    }
    h1, h2, h3 {
        color: var(--ink);
        letter-spacing: 0;
    }
    h1 {font-size: 2rem; margin-bottom: .25rem;}
    h2 {font-size: 1.18rem; margin-top: 1.1rem;}
    h3 {font-size: 1rem;}
    .app-header {
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: 1rem 1.1rem;
        background: #ffffff;
        margin-bottom: .9rem;
    }
    .app-title {
        font-size: 1.72rem;
        line-height: 1.15;
        font-weight: 760;
        color: var(--ink);
        margin: 0 0 .28rem 0;
    }
    .app-copy {
        max-width: 880px;
        color: #475467;
        font-size: .96rem;
        line-height: 1.55;
        margin: 0;
    }
    .flow-row {
        display: flex;
        flex-wrap: wrap;
        gap: .5rem;
        margin-top: .75rem;
    }
    .flow-chip {
        border: 1px solid var(--line);
        border-radius: 999px;
        background: #f8fafc;
        padding: .38rem .62rem;
        color: #344054;
        font-size: .86rem;
    }
    .flow-chip strong {
        color: var(--accent);
        font-weight: 760;
        margin-right: .25rem;
    }
    .section-kicker {
        font-size: .78rem;
        font-weight: 760;
        color: var(--accent);
        letter-spacing: 0;
        margin: .15rem 0 .45rem;
    }
    .status-card {
        border: 1px solid var(--line);
        border-radius: 8px;
        background: #fff;
        padding: .85rem .92rem;
        min-height: 98px;
    }
    .status-card .label {
        color: var(--muted);
        font-size: .85rem;
        margin-bottom: .36rem;
    }
    .status-card .value {
        color: var(--ink);
        font-weight: 760;
        font-size: 1.28rem;
        margin-bottom: .18rem;
    }
    .status-card .note {
        color: var(--muted);
        font-size: .84rem;
        line-height: 1.35;
    }
    .status-ok {color: var(--good); font-weight: 700;}
    .status-warn {color: var(--warn); font-weight: 700;}
    .status-bad {color: var(--bad); font-weight: 700;}
    .muted {color: var(--muted); font-size: .92rem; line-height: 1.55;}
    .download-panel {
        border: 1px solid #b9dfc8;
        background: #f2fbf5;
        border-radius: 8px;
        padding: 1rem;
        margin: .6rem 0 .85rem;
    }
    .download-panel strong {
        color: var(--good);
        display: block;
        margin-bottom: .25rem;
    }
    .preview-panel {
        border: 1px solid var(--line);
        border-radius: 8px;
        background: #fff;
        padding: .8rem .9rem;
        margin: .7rem 0 1rem;
    }
    .asset-list {
        border: 1px solid var(--line);
        border-radius: 8px;
        background: #fbfcfe;
        padding: .7rem .8rem;
        margin: .5rem 0;
    }
    .asset-list strong {
        display: block;
        color: var(--ink);
        margin-bottom: .25rem;
    }
    div[data-testid="stMetric"] {
        background: #fff;
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: .78rem .88rem;
    }
    div[data-testid="stMetric"] label {font-size: .85rem; color: var(--muted);}
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {font-size: 1.34rem; color: var(--ink);}
    .stTabs [data-baseweb="tab-list"] {
        gap: .25rem;
        border-bottom: 1px solid var(--line);
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px 8px 0 0;
        padding: .55rem .8rem;
    }
    [data-testid="stFileUploader"] section {
        background: #f8fafc;
        border: 1px dashed #cfd6e3;
        border-radius: 8px;
        color: var(--ink);
    }
    [data-testid="stFileUploader"] small,
    [data-testid="stFileUploader"] p {
        color: #667085 !important;
    }
    [data-testid="stBaseButton-secondary"],
    [data-testid="stBaseButton-primary"] {
        border-radius: 8px;
    }
    @media (max-width: 900px) {
        .app-title {font-size: 1.42rem;}
    }
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


def result_assets(result_dir: Path) -> tuple[list[Path], list[Path]]:
    if not result_dir.exists():
        return [], []
    html_files = sorted(
        result_dir.glob("**/*.html"),
        key=lambda p: (preview_rank(p), str(p.relative_to(result_dir)).lower()),
    )
    pdf_files = sorted(
        result_dir.glob("**/*.pdf"),
        key=lambda p: (preview_rank(p), str(p.relative_to(result_dir)).lower()),
    )
    return html_files, pdf_files


def preview_rank(path: Path) -> int:
    text = str(path).lower()
    if "02_multi_sequence_alignment" in text or ".alignment" in text:
        return 0
    if "03_sanger_vs_theory" in text or "sanger_vs_theory" in text:
        return 1
    if "01_sample_amplicon_reports" in text or "amplicons" in text:
        return 2
    return 9


def asset_label(path: Path, root: Path) -> str:
    rel = path.relative_to(root)
    label = str(rel).replace("\\", " / ")
    label = label.replace("01_sample_amplicon_reports / html /", "理论扩增 / ")
    label = label.replace("01_sample_amplicon_reports / pdf /", "理论扩增 PDF / ")
    label = label.replace("02_multi_sequence_alignment / html /", "多序列比对 / ")
    label = label.replace("02_multi_sequence_alignment / pdf /", "多序列比对 PDF / ")
    label = label.replace("03_sanger_vs_theory / html /", "测序比对 / ")
    label = label.replace("03_sanger_vs_theory / pdf /", "测序比对 PDF / ")
    return label


def matching_pdf_for_html(html_path: Path, pdf_files: list[Path]) -> Path | None:
    for pdf in pdf_files:
        if pdf.stem == html_path.stem:
            return pdf
    return None


def show_result_preview(result_dir: Path):
    html_files, pdf_files = result_assets(result_dir)
    if not html_files and not pdf_files:
        st.warning("没有找到可预览或可下载的结果文件。")
        return

    st.markdown('<div class="section-kicker">Preview</div>', unsafe_allow_html=True)
    st.subheader("结果预览")
    if html_files:
        selected_html = st.selectbox(
            "选择要查看的结果",
            html_files,
            format_func=lambda p: asset_label(p, result_dir),
        )
        st.markdown(
            f"<div class='preview-panel'><strong>{asset_label(selected_html, result_dir)}</strong><div class='muted'>下方为网页预览，适合快速查看比对结果；正式引用或保存建议下载 PDF。</div></div>",
            unsafe_allow_html=True,
        )
        components.html(selected_html.read_text(encoding="utf-8", errors="replace"), height=760, scrolling=True)
        matched_pdf = matching_pdf_for_html(selected_html, pdf_files)
        if matched_pdf and matched_pdf.exists():
            st.download_button(
                "下载当前预览的 PDF",
                data=matched_pdf.read_bytes(),
                file_name=matched_pdf.name,
                mime="application/pdf",
                use_container_width=True,
                key=f"download-current-{hashlib.sha1(str(matched_pdf).encode()).hexdigest()}",
            )
    else:
        st.info("本次结果没有生成 HTML 预览，但可以下载 PDF。")

    if pdf_files:
        st.markdown('<div class="section-kicker">PDF files</div>', unsafe_allow_html=True)
        st.subheader("单独下载 PDF")
        for index, pdf in enumerate(pdf_files):
            cols = st.columns([3, 1])
            cols[0].markdown(
                f"<div class='asset-list'><strong>{asset_label(pdf, result_dir)}</strong><span class='muted'>{file_size_label(pdf.stat().st_size)}</span></div>",
                unsafe_allow_html=True,
            )
            cols[1].download_button(
                "下载 PDF",
                data=pdf.read_bytes(),
                file_name=pdf.name,
                mime="application/pdf",
                use_container_width=True,
                key=f"download-pdf-{index}-{hashlib.sha1(str(pdf).encode()).hexdigest()}",
            )


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


def status_card(label: str, value: str, note: str):
    st.markdown(
        f"""
        <div class="status-card">
            <div class="label">{label}</div>
            <div class="value">{value}</div>
            <div class="note">{note}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def status_text(ok: bool, optional: bool = False) -> tuple[str, str]:
    if ok:
        return "正常", "status-ok"
    if optional:
        return "可选", "status-warn"
    return "缺失", "status-bad"


cleanup_old_sessions()

st.markdown(
    """
    <div class="app-header">
        <div class="app-title">多序列比对分析平台</div>
        <p class="app-copy">
            上传引物表与样本 FASTA，在线生成理论扩增片段、多样本序列比对、变异位点表和可下载结果包。
        </p>
        <div class="flow-row">
            <span class="flow-chip"><strong>1</strong>上传 Excel 与 FASTA</span>
            <span class="flow-chip"><strong>2</strong>选择分析模块</span>
            <span class="flow-chip"><strong>3</strong>查看比对预览</span>
            <span class="flow-chip"><strong>4</strong>下载 PDF 或完整结果</span>
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.sidebar:
    st.markdown('<div class="section-kicker">数据</div>', unsafe_allow_html=True)
    st.subheader("上传数据")
    primer_upload = st.file_uploader("引物表 Excel", type=["xlsx"])
    fasta_uploads = st.file_uploader(
        "样本 FASTA 文件",
        type=["fa", "fas", "fasta", "fna"],
        accept_multiple_files=True,
    )

    st.markdown('<div class="section-kicker">分析</div>', unsafe_allow_html=True)
    st.subheader("分析内容")
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

    st.markdown('<div class="section-kicker">参数</div>', unsafe_allow_html=True)
    st.subheader("参数设置")
    align_samples = st.text_input("参与比对的样本", "all", help="填写 all，或用英文逗号分隔 FASTA 文件名对应的样本名。")
    number_map = st.text_input("测序编号对应样本编号", DEFAULT_NUMBER_MAP, disabled=not do_sanger)
    flank = st.number_input("两侧延伸长度 bp", min_value=0, max_value=5000, value=500, step=50)
    max_product = st.number_input("最大理论扩增长度 bp", min_value=100, max_value=50000, value=5000, step=100)

    st.markdown('<div class="section-kicker">颜色</div>', unsafe_allow_html=True)
    st.subheader("碱基颜色")
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

mode_label = "完整分析" if do_sanger and do_align else "测序比对" if do_sanger else "多序列比对" if do_align else "理论扩增"
total_upload_mb = sum(item.size for item in fasta_uploads + sanger_uploads)
if primer_upload is not None:
    total_upload_mb += primer_upload.size

top_cols = st.columns(4)
with top_cols[0]:
    status_card("引物表", "已上传" if primer_upload else "等待上传", "Excel 模板将用于读取引物对。")
with top_cols[1]:
    status_card("样本 FASTA", str(len(fasta_preview)), "每个文件作为一个样本参与分析。")
with top_cols[2]:
    status_card("测序文件", str(len(sanger_preview)), "仅在测序比对模块中使用。")
with top_cols[3]:
    status_card("数据体量", file_size_label(total_upload_mb), f"当前运行模式：{mode_label}")

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
        label, css = status_text(ok, optional)
        st.markdown(f"**{name}**：<span class='{css}'>{label}</span>", unsafe_allow_html=True)
        st.caption(note)
    st.markdown(
        "<div class='muted'>每次运行的数据会保存在服务器临时目录中，结果通过压缩包下载；超过 24 小时的临时数据会自动清理。</div>",
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
    st.subheader("运行、预览与下载")
    if st.session_state.get("last_zip"):
        zip_path = Path(st.session_state.last_zip)
        if zip_path.exists():
            st.markdown(
                "<div class='download-panel'><strong>上一次分析结果已就绪</strong><span class='muted'>可以继续下载，也可以重新上传文件开始新的分析。</span></div>",
                unsafe_allow_html=True,
            )
            st.download_button(
                "下载上一次完整结果压缩包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )
            last_result_dir = Path(st.session_state.get("last_result_dir", ""))
            if last_result_dir.exists():
                show_result_preview(last_result_dir)

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

        progress = st.progress(8, text="正在准备分析...")
        with st.spinner("正在分析，请等待。大型数据可能需要几分钟。"):
            progress.progress(35, text="正在运行分析程序...")
            code, output = run_command(command)
            progress.progress(100, text="分析结束")

        with st.expander("查看运行日志", expanded=code != 0):
            st.text_area("运行日志", output, height=260)
        if code == 0:
            zip_path = make_result_zip(result_dir)
            st.session_state.last_result_dir = str(result_dir)
            st.session_state.last_zip = str(zip_path)
            st.success("分析完成")
            show_result_preview(result_dir)
            st.markdown(
                "<div class='download-panel'><strong>结果已生成</strong><span class='muted'>压缩包包含 PDF、CSV 表格、序列比对文件和可选 Word 报告。</span></div>",
                unsafe_allow_html=True,
            )
            st.download_button(
                "下载完整结果压缩包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )
        else:
            st.error("分析失败，请根据日志检查输入文件或参数。")
