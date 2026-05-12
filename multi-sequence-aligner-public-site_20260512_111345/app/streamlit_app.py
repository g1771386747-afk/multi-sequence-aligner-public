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
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    :root {
        --bg: #f5f7fa;
        --panel: #ffffff;
        --ink: #182230;
        --muted: #667085;
        --line: #dfe5ee;
        --soft: #f8fafc;
        --accent: #146c7c;
        --accent-soft: #eaf6f8;
        --green: #087443;
        --amber: #9a5b00;
        --red: #b42318;
    }
    .stApp {
        background: var(--bg);
        color: var(--ink);
    }
    header[data-testid="stHeader"] {
        background: rgba(245, 247, 250, .92);
    }
    [data-testid="stSidebar"] {
        display: none;
    }
    .block-container {
        padding-top: 1.25rem;
        padding-bottom: 3rem;
        max-width: 1320px;
    }
    h1, h2, h3, p {
        letter-spacing: 0;
    }
    .topbar {
        display: flex;
        justify-content: space-between;
        gap: 1rem;
        align-items: flex-start;
        border: 1px solid var(--line);
        border-radius: 8px;
        background: var(--panel);
        padding: 1.25rem 1.35rem;
        margin-bottom: 1rem;
    }
    .brand-title {
        font-size: 1.78rem;
        line-height: 1.15;
        font-weight: 780;
        color: var(--ink);
        margin: 0 0 .36rem;
    }
    .brand-copy {
        color: #475467;
        font-size: .98rem;
        line-height: 1.56;
        margin: 0;
        max-width: 850px;
    }
    .privacy-note {
        white-space: nowrap;
        border: 1px solid #c7e5db;
        background: #f2fbf7;
        color: #087443;
        border-radius: 999px;
        padding: .38rem .68rem;
        font-size: .84rem;
        font-weight: 700;
    }
    .workflow {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: .7rem;
        margin-bottom: 1rem;
    }
    .workflow-step {
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: .75rem .85rem;
    }
    .workflow-step span {
        display: inline-flex;
        width: 1.45rem;
        height: 1.45rem;
        align-items: center;
        justify-content: center;
        background: var(--accent-soft);
        color: var(--accent);
        border-radius: 999px;
        font-weight: 800;
        font-size: .82rem;
        margin-right: .42rem;
    }
    .workflow-step strong {
        color: var(--ink);
        font-size: .95rem;
    }
    .section-card {
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .section-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 1rem;
        margin-bottom: .75rem;
    }
    .section-title {
        font-size: 1.08rem;
        font-weight: 780;
        color: var(--ink);
        margin: 0;
    }
    .section-subtitle {
        color: var(--muted);
        font-size: .88rem;
        margin: .12rem 0 0;
    }
    .status-pill {
        border-radius: 999px;
        padding: .28rem .58rem;
        font-size: .8rem;
        font-weight: 760;
        border: 1px solid var(--line);
        background: var(--soft);
        color: var(--muted);
    }
    .summary-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: .7rem;
        margin-bottom: 1rem;
    }
    .summary-card {
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: .78rem .85rem;
    }
    .summary-card .label {
        color: var(--muted);
        font-size: .82rem;
        margin-bottom: .28rem;
    }
    .summary-card .value {
        color: var(--ink);
        font-size: 1.34rem;
        font-weight: 780;
        line-height: 1.15;
    }
    .summary-card .note {
        color: var(--muted);
        font-size: .82rem;
        margin-top: .24rem;
    }
    .ok { color: var(--green); font-weight: 760; }
    .warn { color: var(--amber); font-weight: 760; }
    .bad { color: var(--red); font-weight: 760; }
    .muted {
        color: var(--muted);
        font-size: .9rem;
        line-height: 1.55;
    }
    .result-ready {
        border: 1px solid #b9dfc8;
        background: #f2fbf5;
        border-radius: 8px;
        padding: .9rem 1rem;
        margin: .6rem 0 .85rem;
    }
    .result-ready strong {
        color: var(--green);
        display: block;
        margin-bottom: .22rem;
    }
    .asset-row {
        border: 1px solid var(--line);
        background: #fbfcfe;
        border-radius: 8px;
        padding: .62rem .72rem;
        margin-bottom: .45rem;
    }
    .asset-row strong {
        color: var(--ink);
        display: block;
        font-size: .9rem;
    }
    div[data-testid="stFileUploader"] section {
        background: #fbfcfe;
        border: 1px dashed #bcc8d8;
        border-radius: 8px;
        min-height: 116px;
        color: var(--ink);
    }
    div[data-testid="stFileUploader"] button {
        border-radius: 8px;
    }
    div[data-testid="stFileUploader"] small,
    div[data-testid="stFileUploader"] p {
        color: #667085 !important;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: .2rem;
        border-bottom: 1px solid var(--line);
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px 8px 0 0;
        padding: .5rem .75rem;
    }
    .stButton > button,
    .stDownloadButton > button {
        border-radius: 8px;
    }
    @media (max-width: 980px) {
        .topbar {flex-direction: column;}
        .privacy-note {white-space: normal;}
        .workflow, .summary-grid {grid-template-columns: repeat(2, minmax(0, 1fr));}
    }
    @media (max-width: 640px) {
        .workflow, .summary-grid {grid-template-columns: 1fr;}
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
    except Exception as exc:  # pragma: no cover
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


def preview_rank(path: Path) -> int:
    text = str(path).lower()
    if "02_multi_sequence_alignment" in text or ".alignment" in text:
        return 0
    if "03_sanger_vs_theory" in text or "sanger_vs_theory" in text:
        return 1
    if "01_sample_amplicon_reports" in text or "amplicons" in text:
        return 2
    return 9


def result_assets(result_dir: Path) -> tuple[list[Path], list[Path]]:
    if not result_dir.exists():
        return [], []
    html_files = sorted(result_dir.glob("**/*.html"), key=lambda p: (preview_rank(p), str(p).lower()))
    pdf_files = sorted(result_dir.glob("**/*.pdf"), key=lambda p: (preview_rank(p), str(p).lower()))
    return html_files, pdf_files


def asset_label(path: Path, root: Path) -> str:
    label = str(path.relative_to(root)).replace("\\", " / ")
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

    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("结果预览")
    if html_files:
        selected_html = st.selectbox("选择要查看的结果", html_files, format_func=lambda p: asset_label(p, result_dir))
        components.html(selected_html.read_text(encoding="utf-8", errors="replace"), height=720, scrolling=True)
        matched_pdf = matching_pdf_for_html(selected_html, pdf_files)
        if matched_pdf and matched_pdf.exists():
            st.download_button(
                "下载当前预览 PDF",
                data=matched_pdf.read_bytes(),
                file_name=matched_pdf.name,
                mime="application/pdf",
                use_container_width=True,
                key=f"download-current-{hashlib.sha1(str(matched_pdf).encode()).hexdigest()}",
            )
    else:
        st.info("本次结果没有生成 HTML 预览，但可以下载 PDF。")

    if pdf_files:
        st.markdown("##### 单独下载 PDF")
        for index, pdf in enumerate(pdf_files):
            cols = st.columns([3.2, 1])
            cols[0].markdown(
                f"<div class='asset-row'><strong>{asset_label(pdf, result_dir)}</strong><span class='muted'>{file_size_label(pdf.stat().st_size)}</span></div>",
                unsafe_allow_html=True,
            )
            cols[1].download_button(
                "下载",
                data=pdf.read_bytes(),
                file_name=pdf.name,
                mime="application/pdf",
                use_container_width=True,
                key=f"download-pdf-{index}-{hashlib.sha1(str(pdf).encode()).hexdigest()}",
            )
    st.markdown("</div>", unsafe_allow_html=True)


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


def file_table(uploaded_files) -> list[dict[str, str]]:
    return [{"文件名": item.name, "样本": Path(item.name).stem, "大小": file_size_label(item.size)} for item in uploaded_files]


cleanup_old_sessions()

st.markdown(
    """
    <div class="topbar">
        <div>
            <div class="brand-title">多序列比对分析平台</div>
            <p class="brand-copy">面向叶绿体/基因组片段分析的在线工具。上传引物表和样本 FASTA，即可生成理论扩增片段、多序列比对图、变异位点表，并按需下载 PDF。</p>
        </div>
        <div class="privacy-note">公共上传工具</div>
    </div>
    <div class="workflow">
        <div class="workflow-step"><span>1</span><strong>上传数据</strong></div>
        <div class="workflow-step"><span>2</span><strong>选择模块</strong></div>
        <div class="workflow-step"><span>3</span><strong>预览结果</strong></div>
        <div class="workflow-step"><span>4</span><strong>下载报告</strong></div>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.container():
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.markdown(
        """
        <div class="section-head">
            <div>
                <div class="section-title">1. 上传数据</div>
                <div class="section-subtitle">引物表为 Excel；FASTA 可一次上传多个样本。测序文件仅在启用测序比对时需要。</div>
            </div>
            <div class="status-pill">支持 xlsx / fa / fasta / seq / txt</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    up1, up2, up3 = st.columns([1, 1.35, 1])
    with up1:
        primer_upload = st.file_uploader("引物表 Excel", type=["xlsx"])
    with up2:
        fasta_uploads = st.file_uploader(
            "样本 FASTA 文件",
            type=["fa", "fas", "fasta", "fna"],
            accept_multiple_files=True,
        )
    with up3:
        sanger_uploads = st.file_uploader(
            "胶回收测序序列文件",
            type=["fa", "fas", "fasta", "seq", "txt"],
            accept_multiple_files=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

fasta_uploads = fasta_uploads or []
sanger_uploads = sanger_uploads or []

primer_preview_path = cache_primer_for_preview(primer_upload)
primer_rows, primer_error = read_primer_preview(primer_preview_path)
fasta_preview = file_table(fasta_uploads)
sanger_preview = [{"文件名": item.name, "大小": file_size_label(item.size)} for item in sanger_uploads]
total_size = sum(item.size for item in fasta_uploads + sanger_uploads) + (primer_upload.size if primer_upload else 0)

mode_cols = st.columns([1.3, 1])
with mode_cols[0]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.markdown(
        """
        <div class="section-head">
            <div>
                <div class="section-title">2. 分析模块</div>
                <div class="section-subtitle">默认生成理论扩增片段和多样本序列比对；测序验证可按需开启。</div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    m1, m2, m3 = st.columns(3)
    with m1:
        do_stage1 = st.checkbox("理论扩增片段 PDF", value=True)
    with m2:
        do_align = st.checkbox("多样本序列比对 PDF", value=True)
    with m3:
        do_sanger = st.checkbox("胶回收测序验证", value=False)
    m4, m5 = st.columns(2)
    with m4:
        do_word = st.checkbox("生成 Word 汇总报告", value=False, disabled=not do_sanger)
    with m5:
        zip_pdfs = st.checkbox("生成 PDF 压缩包", value=True)
    st.markdown("</div>", unsafe_allow_html=True)

with mode_cols[1]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.markdown(
        """
        <div class="section-head">
            <div>
                <div class="section-title">3. 参数设置</div>
                <div class="section-subtitle">常规分析保持默认值即可。</div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    align_samples = st.text_input("参与比对的样本", "all")
    number_map = st.text_input("测序编号对应样本编号", DEFAULT_NUMBER_MAP, disabled=not do_sanger)
    p1, p2 = st.columns(2)
    with p1:
        flank = st.number_input("两侧延伸 bp", min_value=0, max_value=5000, value=500, step=50)
    with p2:
        max_product = st.number_input("最大扩增长度 bp", min_value=100, max_value=50000, value=5000, step=100)
    st.markdown("</div>", unsafe_allow_html=True)

color_expander = st.expander("碱基颜色设置", expanded=False)
with color_expander:
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        color_a = st.color_picker("A", "#1aa260")
    with c2:
        color_t = st.color_picker("T", "#d83b36")
    with c3:
        color_c = st.color_picker("C", "#2568c9")
    with c4:
        color_g = st.color_picker("G", "#111111")
    with c5:
        color_gap = st.color_picker("gap", "#d9d9d9")

mode_label = "完整分析" if do_sanger and do_align else "测序验证" if do_sanger else "多序列比对" if do_align else "理论扩增"
st.markdown(
    f"""
    <div class="summary-grid">
        <div class="summary-card"><div class="label">引物表</div><div class="value">{'已上传' if primer_upload else '缺失'}</div><div class="note">解析到 {0 if primer_error else len(primer_rows)} 对引物</div></div>
        <div class="summary-card"><div class="label">样本 FASTA</div><div class="value">{len(fasta_preview)}</div><div class="note">每个 FASTA 作为一个样本</div></div>
        <div class="summary-card"><div class="label">测序文件</div><div class="value">{len(sanger_preview)}</div><div class="note">测序验证模块使用</div></div>
        <div class="summary-card"><div class="label">数据体量</div><div class="value">{file_size_label(total_size)}</div><div class="note">当前模式：{mode_label}</div></div>
    </div>
    """,
    unsafe_allow_html=True,
)

tab_check, tab_primer, tab_fasta, tab_sanger, tab_run = st.tabs(["检查", "引物表", "FASTA", "测序文件", "运行与结果"])

with tab_check:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("输入检查")
    checks = [
        ("引物表", primer_upload is not None, "必须上传 Excel 引物表。"),
        ("FASTA 文件", bool(fasta_uploads), "至少上传一个样本 FASTA。"),
        ("测序文件", bool(sanger_uploads), "启用测序验证时需要上传。"),
    ]
    for name, ok, note in checks:
        optional = name == "测序文件" and not do_sanger
        css = "ok" if ok else ("warn" if optional else "bad")
        text = "正常" if ok else ("可选" if optional else "缺失")
        st.markdown(f"**{name}**：<span class='{css}'>{text}</span>", unsafe_allow_html=True)
        st.caption(note)
    st.markdown("<div class='muted'>服务器临时数据会自动清理。建议用户在分析结束后及时下载结果包。</div>", unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)

with tab_primer:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("引物表预览")
    if primer_error:
        st.info(primer_error)
    else:
        st.dataframe(primer_rows, use_container_width=True, hide_index=True)
    st.markdown("</div>", unsafe_allow_html=True)

with tab_fasta:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("FASTA 文件")
    if fasta_preview:
        st.dataframe(fasta_preview, use_container_width=True, hide_index=True)
    else:
        st.info("请上传样本 FASTA 文件。")
    st.markdown("</div>", unsafe_allow_html=True)

with tab_sanger:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("测序文件")
    if sanger_preview:
        st.dataframe(sanger_preview, use_container_width=True, hide_index=True)
    else:
        st.info("未上传测序文件；未启用测序验证时可以忽略。")
    st.markdown("</div>", unsafe_allow_html=True)

with tab_run:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("运行与结果")
    r1, r2 = st.columns([1, 1])
    with r1:
        run = st.button("开始分析", type="primary", use_container_width=True)
    with r2:
        clear = st.button("清除本次上传与结果", use_container_width=True)
    st.markdown("</div>", unsafe_allow_html=True)

    if clear:
        clear_session_files()
        st.session_state.pop("last_zip", None)
        st.session_state.pop("last_result_dir", None)
        st.rerun()

    if st.session_state.get("last_zip"):
        zip_path = Path(st.session_state.last_zip)
        if zip_path.exists():
            st.markdown(
                "<div class='result-ready'><strong>上一次分析结果已就绪</strong><span class='muted'>可以继续下载，也可以重新上传文件开始新的分析。</span></div>",
                unsafe_allow_html=True,
            )
            st.download_button(
                "下载上一次完整结果包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )
            last_result_dir = Path(st.session_state.get("last_result_dir", ""))
            if last_result_dir.exists():
                show_result_preview(last_result_dir)

    if run:
        if not (do_stage1 or do_align or do_sanger):
            st.error("请至少选择一个分析模块。")
            st.stop()
        if primer_upload is None:
            st.error("请上传引物表 Excel。")
            st.stop()
        if not fasta_uploads:
            st.error("请至少上传一个样本 FASTA 文件。")
            st.stop()
        if do_sanger and not sanger_uploads:
            st.error("已启用测序验证，请上传胶回收测序序列文件。")
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
                "<div class='result-ready'><strong>完整结果包已生成</strong><span class='muted'>压缩包包含 PDF、CSV 表格、序列比对文件和可选 Word 报告。</span></div>",
                unsafe_allow_html=True,
            )
            st.download_button(
                "下载完整结果包",
                data=zip_path.read_bytes(),
                file_name=zip_path.name,
                mime="application/zip",
                use_container_width=True,
            )
        else:
            st.error("分析失败，请根据日志检查输入文件或参数。")
