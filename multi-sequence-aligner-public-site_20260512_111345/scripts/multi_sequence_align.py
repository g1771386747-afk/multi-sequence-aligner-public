#!/usr/bin/env python
"""Multi-sequence PCR amplicon workflow.

This script intentionally uses only the Python standard library plus Microsoft
Edge for PDF printing when available.
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
import shutil
import subprocess
import sys
import textwrap
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET


SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_GENOME_DIR = Path.home() / ".codex" / "skills" / "chloroplast-amplicon-aligner" / "data" / "01gbf"
DEFAULT_EDGE = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
DESKTOP = Path.home() / "Desktop"


BASE_COLORS = {
    "A": "#1aa260",
    "T": "#d83b36",
    "C": "#2568c9",
    "G": "#111111",
    "-": "#d9d9d9",
    "N": "#eeeeee",
}


@dataclass
class PrimerPair:
    order: int
    pair_name: str
    f_name: str
    f_seq: str
    r_name: str
    r_seq: str


@dataclass
class FastaRecord:
    sample: str
    path: Path
    header: str
    seq: str


@dataclass
class AmpliconHit:
    sample: str
    pair_name: str
    pair_order: int
    status: str
    strand: str = ""
    product_start: int = 0
    product_end: int = 0
    product_len: int = 0
    window_start: int = 0
    window_end: int = 0
    f_start: int = 0
    f_end: int = 0
    r_start: int = 0
    r_end: int = 0
    product_seq: str = ""
    window_seq: str = ""
    message: str = ""


def safe_name(text: str) -> str:
    text = re.sub(r"[\\/:*?\"<>|]+", "_", str(text))
    text = re.sub(r"\s+", "_", text).strip("_")
    return text or "unnamed"


def run_root(user_out: str | None = None) -> Path:
    if user_out:
        root = Path(user_out)
        root.mkdir(parents=True, exist_ok=True)
        return root
    root = DESKTOP / f"多序列比对_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
    root.mkdir(parents=True, exist_ok=False)
    return root


def colnum(cell_ref: str) -> int:
    m = re.match(r"([A-Z]+)", cell_ref)
    if not m:
        return 0
    n = 0
    for ch in m.group(1):
        n = n * 26 + ord(ch) - 64
    return n


def read_xlsx_sheet1(path: Path) -> dict[int, dict[int, str]]:
    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    rows: dict[int, dict[int, str]] = {}
    with zipfile.ZipFile(path) as z:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", ns):
                shared.append("".join(t.text or "" for t in si.findall(".//a:t", ns)))
        workbook = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        sheet = workbook.find(".//a:sheet", ns)
        if sheet is None:
            return rows
        rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        target = relmap[rid]
        sheet_path = "xl/" + target.lstrip("/") if not target.startswith("xl/") else target
        root = ET.fromstring(z.read(sheet_path))
        for row in root.findall(".//a:sheetData/a:row", ns):
            vals: dict[int, str] = {}
            for c in row.findall("a:c", ns):
                v = c.find("a:v", ns)
                val = "" if v is None else (v.text or "")
                if c.attrib.get("t") == "s" and val:
                    val = shared[int(val)]
                elif c.attrib.get("t") == "inlineStr":
                    val = "".join(t.text or "" for t in c.findall(".//a:t", ns))
                vals[colnum(c.attrib.get("r", ""))] = str(val).strip()
            rows[int(row.attrib["r"])] = vals
    return rows


def read_primers(path: Path) -> list[PrimerPair]:
    rows = read_xlsx_sheet1(path)
    primers: list[PrimerPair] = []
    for r in sorted(rows):
        if r == 1:
            continue
        vals = rows[r]
        pair_name = vals.get(1, "").strip()
        f_name = vals.get(2, "").strip()
        f_seq = clean_seq(vals.get(3, ""))
        r_name = vals.get(4, "").strip()
        r_seq = clean_seq(vals.get(5, ""))
        if pair_name and f_seq and r_seq:
            if not f_name:
                f_name = f"{pair_name}-F"
            if not r_name:
                r_name = f"{pair_name}-R"
            primers.append(PrimerPair(len(primers) + 1, pair_name, f_name, f_seq, r_name, r_seq))
    return primers


def clean_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def read_fasta(path: Path) -> FastaRecord:
    header = path.stem
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if not parts:
                header = line[1:].strip() or path.stem
        else:
            parts.append(clean_seq(line))
    return FastaRecord(path.stem, path, header, "".join(parts))


def read_fastas(paths: list[str] | None, genome_dir: str | None) -> list[FastaRecord]:
    files: list[Path] = []
    if paths:
        for item in paths:
            p = Path(item)
            if p.is_dir():
                files.extend(sorted(p.glob("*.fa*")))
            else:
                files.append(p)
    else:
        gdir = Path(genome_dir) if genome_dir else DEFAULT_GENOME_DIR
        files.extend(sorted(gdir.glob("*.fa*"), key=lambda p: natural_key(p.stem)))
    return [read_fasta(p) for p in files if p.exists() and p.is_file()]


def natural_key(text: str):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(text))]


_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1].upper()


def all_indices(seq: str, needle: str) -> list[int]:
    out: list[int] = []
    start = 0
    while True:
        i = seq.find(needle, start)
        if i < 0:
            return out
        out.append(i)
        start = i + 1


def circular_slice(seq: str, start1: int, length: int) -> str:
    n = len(seq)
    if n == 0 or length <= 0:
        return ""
    start0 = (start1 - 1) % n
    doubled = seq + seq
    if start0 + length <= len(doubled):
        return doubled[start0 : start0 + length]
    reps = (length // n) + 2
    return (seq * reps)[start0 : start0 + length]


def coord_add(pos1: int, offset: int, n: int) -> int:
    return ((pos1 - 1 + offset) % n) + 1


def find_amplicons(record: FastaRecord, primer: PrimerPair, flank: int, max_product: int) -> list[AmpliconHit]:
    seq = record.seq
    n = len(seq)
    if not seq:
        return [AmpliconHit(record.sample, primer.pair_name, primer.order, "NOT_FOUND", message="empty sequence")]
    hits: list[AmpliconHit] = []

    def add_hit(strand: str, f0: int, r0: int, product_len: int, f_len: int, r_len: int, product_seq: str):
        product_start = f0 + 1
        product_end = coord_add(product_start, product_len - 1, n)
        window_start = coord_add(product_start, -flank, n)
        window_len = product_len + flank * 2
        window_end = coord_add(window_start, window_len - 1, n)
        window_seq = circular_slice(seq, window_start, window_len)
        if strand == "+":
            f_start = product_start
            f_end = coord_add(f_start, f_len - 1, n)
            r_start = coord_add(product_start, product_len - r_len, n)
            r_end = product_end
        else:
            # Coordinates are on the original genome; displayed product sequence is F-to-R orientation.
            r_start = product_start
            r_end = coord_add(r_start, r_len - 1, n)
            f_start = coord_add(product_start, product_len - f_len, n)
            f_end = product_end
        hits.append(
            AmpliconHit(
                sample=record.sample,
                pair_name=primer.pair_name,
                pair_order=primer.order,
                status="FOUND",
                strand=strand,
                product_start=product_start,
                product_end=product_end,
                product_len=product_len,
                window_start=window_start,
                window_end=window_end,
                f_start=f_start,
                f_end=f_end,
                r_start=r_start,
                r_end=r_end,
                product_seq=product_seq,
                window_seq=window_seq,
            )
        )

    seq2 = seq + seq
    f_seq = primer.f_seq
    r_seq = primer.r_seq
    rc_r = revcomp(r_seq)
    rc_f = revcomp(f_seq)

    for f0 in all_indices(seq, f_seq):
        search_start = f0 + len(f_seq)
        search_end = min(f0 + max_product, f0 + n)
        cursor = search_start
        while True:
            r0 = seq2.find(rc_r, cursor, search_end)
            if r0 < 0:
                break
            product_len = r0 + len(rc_r) - f0
            product_seq = circular_slice(seq, f0 + 1, product_len)
            add_hit("+", f0, r0 % n, product_len, len(f_seq), len(r_seq), product_seq)
            cursor = r0 + 1

    for r0 in all_indices(seq, r_seq):
        search_start = r0 + len(r_seq)
        search_end = min(r0 + max_product, r0 + n)
        cursor = search_start
        while True:
            f0 = seq2.find(rc_f, cursor, search_end)
            if f0 < 0:
                break
            product_len = f0 + len(rc_f) - r0
            raw_product = circular_slice(seq, r0 + 1, product_len)
            product_seq = revcomp(raw_product)
            add_hit("-", r0, f0 % n, product_len, len(f_seq), len(r_seq), product_seq)
            cursor = f0 + 1

    if not hits:
        return [AmpliconHit(record.sample, primer.pair_name, primer.order, "NOT_FOUND", message="no paired primer hit")]
    if len(hits) > 1:
        for h in hits:
            h.status = "WARN_MULTIPLE_HITS"
            h.message = f"{len(hits)} candidate products"
    return hits


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def xml_escape(text) -> str:
    return html.escape(str(text), quote=True)


def make_docx(path: Path, title: str, sections: list[tuple[str, list[str], list[list[str]] | None]]):
    """Create a simple Word-compatible .docx without external dependencies.

    sections contains (heading, paragraphs, table_rows). table_rows includes the
    header row as the first row.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    def p(text: str, style: str | None = None) -> str:
        style_xml = f'<w:pStyle w:val="{style}"/>' if style else ""
        return f"<w:p><w:pPr>{style_xml}</w:pPr><w:r><w:t>{xml_escape(text)}</w:t></w:r></w:p>"

    def table(rows: list[list[str]]) -> str:
        if not rows:
            return ""
        trs = []
        for row in rows:
            cells = []
            for cell in row:
                cells.append(f"<w:tc><w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/></w:tcPr>{p(str(cell))}</w:tc>")
            trs.append("<w:tr>" + "".join(cells) + "</w:tr>")
        return "<w:tbl><w:tblPr><w:tblW w:w=\"0\" w:type=\"auto\"/><w:tblBorders><w:top w:val=\"single\" w:sz=\"4\"/><w:left w:val=\"single\" w:sz=\"4\"/><w:bottom w:val=\"single\" w:sz=\"4\"/><w:right w:val=\"single\" w:sz=\"4\"/><w:insideH w:val=\"single\" w:sz=\"4\"/><w:insideV w:val=\"single\" w:sz=\"4\"/></w:tblBorders></w:tblPr>" + "".join(trs) + "</w:tbl>"

    body = [p(title, "Title")]
    for heading, paragraphs, rows in sections:
        body.append(p(heading, "Heading1"))
        for para in paragraphs:
            body.append(p(para))
        if rows:
            body.append(table(rows))
    document = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        "<w:body>"
        + "".join(body)
        + '<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"/></w:sectPr>'
        "</w:body></w:document>"
    )
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        '<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:rPr><w:b/><w:sz w:val="24"/></w:rPr></w:style>'
        "</w:styles>"
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>')
        z.writestr("_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')
        z.writestr("word/_rels/document.xml.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>')
        z.writestr("word/document.xml", document)
        z.writestr("word/styles.xml", styles)


def html_doc(title: str, body: str, landscape: bool = False) -> str:
    page = "A4 landscape" if landscape else "A4"
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>
body{{font-family:Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#222}}
h1{{font-size:24px}} h2{{font-size:18px;margin-top:22px}} h3{{font-size:15px;margin:14px 0 6px}}
table{{border-collapse:collapse;font-size:12px;margin:8px 0;width:100%}}
th,td{{border:1px solid #aaa;padding:4px 6px;vertical-align:top}} th{{background:#eee}}
.page{{page-break-after:always}} .warn{{color:#b42318;font-weight:700}} .ok{{color:#067647;font-weight:700}}
.seq{{font-family:Consolas,monospace;font-size:11px;line-height:1.45;white-space:pre-wrap;word-break:break-all}}
.legend span{{display:inline-block;margin-right:14px}} .chip{{display:inline-block;width:13px;height:13px;border:1px solid #777;vertical-align:-2px;margin-right:4px}}
@page{{size:{page};margin:12mm}}
</style></head><body>{body}</body></html>"""


def find_pdf_browser() -> str | None:
    for env_name in ["MULTISEQ_BROWSER", "CHROME_BIN", "CHROMIUM_BIN", "BROWSER_PATH"]:
        value = os.environ.get(env_name)
        if value and Path(value).exists():
            return value
    candidates = [
        DEFAULT_EDGE,
        Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"),
        Path("/usr/bin/chromium"),
        Path("/usr/bin/chromium-browser"),
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    for name in ["chromium", "chromium-browser", "google-chrome", "google-chrome-stable", "msedge"]:
        resolved = shutil.which(name)
        if resolved:
            return resolved
    return None


def edge_print_pdf(html_path: Path, pdf_path: Path):
    browser = find_pdf_browser()
    if not browser:
        print(f"WARN: No Chromium/Chrome/Edge browser found; PDF not created for {html_path}", file=sys.stderr)
        return
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            browser,
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--disable-dev-shm-usage",
            f"--print-to-pdf={pdf_path}",
            html_path.as_uri(),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        timeout=90,
    )


def region_track_svg(hit: AmpliconHit, primer: PrimerPair, width: int = 900) -> str:
    if not hit.product_len:
        return ""
    flank = max(0, (len(hit.window_seq) - hit.product_len) // 2)
    total = len(hit.window_seq)
    left = 40
    right = 40
    track_w = width - left - right

    def x(offset: int) -> float:
        return left + (offset / max(1, total)) * track_w

    f_offset = flank
    r_offset = flank + hit.product_len - len(primer.r_seq)
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='110' viewBox='0 0 {width} 110'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<line x1='{left}' y1='50' x2='{width-right}' y2='50' stroke='#3366aa' stroke-width='2'/>",
        f"<rect x='{x(flank)}' y='38' width='{max(2, x(flank + hit.product_len)-x(flank))}' height='24' fill='#e8f0fe' stroke='#7aa5ff'/>",
        f"<rect x='{x(f_offset)}' y='30' width='{max(3, x(f_offset+len(primer.f_seq))-x(f_offset))}' height='40' fill='#18a058' opacity='0.85'/>",
        f"<rect x='{x(r_offset)}' y='30' width='{max(3, x(r_offset+len(primer.r_seq))-x(r_offset))}' height='40' fill='#d83b36' opacity='0.85'/>",
        f"<text x='{left}' y='22' font-size='12' font-family='Arial'>window {hit.window_start}-{hit.window_end}</text>",
        f"<text x='{x(flank)}' y='84' font-size='12' font-family='Arial'>PCR {hit.product_start}-{hit.product_end}, {hit.product_len} bp</text>",
        f"<text x='{x(f_offset)}' y='28' font-size='11' font-family='Arial'>F</text>",
        f"<text x='{x(r_offset)}' y='28' font-size='11' font-family='Arial'>R</text>",
        "</svg>",
    ]
    return "\n".join(parts)


def marked_sequence(hit: AmpliconHit, primer: PrimerPair, genome_len: int, line_width: int = 100) -> str:
    seq = hit.window_seq
    if not seq:
        return ""
    flank = max(0, (len(seq) - hit.product_len) // 2)
    f_range = range(flank, flank + len(primer.f_seq))
    r_start = flank + hit.product_len - len(primer.r_seq)
    r_range = range(r_start, r_start + len(primer.r_seq))
    product_range = range(flank, flank + hit.product_len)
    out: list[str] = []
    for start in range(0, len(seq), line_width):
        chunk = seq[start : start + line_width]
        pieces: list[str] = []
        for i, b in enumerate(chunk, start):
            esc = html.escape(b)
            if i in f_range:
                pieces.append(f"<span style='background:#b7ebc6;color:#064e3b;font-weight:700'>{esc}</span>")
            elif i in r_range:
                pieces.append(f"<span style='background:#ffd0d0;color:#7f1d1d;font-weight:700'>{esc}</span>")
            elif i in product_range:
                pieces.append(f"<span style='background:#e8f0fe'>{esc}</span>")
            else:
                pieces.append(esc)
        coord_start = coord_add(hit.window_start, start, genome_len)
        coord_end = coord_add(hit.window_start, start + len(chunk) - 1, genome_len)
        out.append(f"{coord_start:>6}-{coord_end:<6}  {''.join(pieces)}")
    return "\n".join(out)


def stage1_sample_reports(root: Path, primers: list[PrimerPair], records: list[FastaRecord], flank: int, max_product: int) -> dict[tuple[str, str], list[AmpliconHit]]:
    stage = root / "01_sample_amplicon_reports"
    html_dir = stage / "html"
    pdf_dir = stage / "pdf"
    csv_dir = stage / "tables"
    all_hits: dict[tuple[str, str], list[AmpliconHit]] = {}
    summary_rows: list[dict] = []

    for rec in records:
        body = [f"<h1>{html.escape(rec.sample)} 理论PCR扩增区域报告</h1>"]
        body.append(f"<p>Genome length: {len(rec.seq)} bp; source: {html.escape(str(rec.path))}</p>")
        for primer in primers:
            hits = find_amplicons(rec, primer, flank, max_product)
            all_hits[(rec.sample, primer.pair_name)] = hits
            body.append(f"<div class='page'><h2>{primer.order}. {html.escape(primer.pair_name)}</h2>")
            body.append(
                "<table><tbody>"
                f"<tr><th>F primer</th><td>{html.escape(primer.f_name)}: {primer.f_seq}</td></tr>"
                f"<tr><th>R primer</th><td>{html.escape(primer.r_name)}: {primer.r_seq}</td></tr>"
                "</tbody></table>"
            )
            for idx, hit in enumerate(hits, 1):
                row = hit_to_row(hit, primer, rec.path)
                summary_rows.append(row)
                cls = "ok" if hit.status == "FOUND" else "warn"
                body.append(f"<h3 class='{cls}'>Hit {idx}: {html.escape(hit.status)} {html.escape(hit.message)}</h3>")
                if hit.product_len:
                    body.append(region_track_svg(hit, primer))
                    body.append("<table><tbody>")
                    body.append(f"<tr><th>Strand</th><td>{hit.strand}</td><th>Product</th><td>{hit.product_start}-{hit.product_end}, {hit.product_len} bp</td></tr>")
                    body.append(f"<tr><th>Window</th><td>{hit.window_start}-{hit.window_end}</td><th>F/R coords</th><td>F {hit.f_start}-{hit.f_end}; R {hit.r_start}-{hit.r_end}</td></tr>")
                    body.append("</tbody></table>")
                    body.append(f"<div class='seq'>{marked_sequence(hit, primer, len(rec.seq))}</div>")
                else:
                    body.append("<p class='warn'>未找到理论扩增片段。</p>")
            body.append("</div>")
        html_path = html_dir / f"{safe_name(rec.sample)}_amplicons.html"
        pdf_path = pdf_dir / f"{safe_name(rec.sample)}_amplicons.pdf"
        html_path.parent.mkdir(parents=True, exist_ok=True)
        html_path.write_text(html_doc(f"{rec.sample} amplicons", "\n".join(body)), encoding="utf-8")
        edge_print_pdf(html_path, pdf_path)

    fields = [
        "sample",
        "source_fasta",
        "primer_order",
        "primer_pair",
        "status",
        "strand",
        "product_start",
        "product_end",
        "product_len",
        "window_start",
        "window_end",
        "f_start",
        "f_end",
        "r_start",
        "r_end",
        "message",
    ]
    write_csv(csv_dir / "all_sample_amplicon_locations.csv", summary_rows, fields)
    return all_hits


def hit_to_row(hit: AmpliconHit, primer: PrimerPair, source: Path) -> dict:
    return {
        "sample": hit.sample,
        "source_fasta": str(source),
        "primer_order": primer.order,
        "primer_pair": primer.pair_name,
        "status": hit.status,
        "strand": hit.strand,
        "product_start": hit.product_start or "",
        "product_end": hit.product_end or "",
        "product_len": hit.product_len or "",
        "window_start": hit.window_start or "",
        "window_end": hit.window_end or "",
        "f_start": hit.f_start or "",
        "f_end": hit.f_end or "",
        "r_start": hit.r_start or "",
        "r_end": hit.r_end or "",
        "message": hit.message,
    }


def nw(ref: str, seq: str) -> tuple[str, str]:
    m, n = len(ref), len(seq)
    gap, match, mismatch = -2, 2, -1
    prev = [j * gap for j in range(n + 1)]
    trace: list[list[str]] = []
    for i in range(1, m + 1):
        cur = [i * gap] + [0] * n
        drow = ["U"] * (n + 1)
        ri = ref[i - 1]
        for j in range(1, n + 1):
            vals = (
                prev[j - 1] + (match if ri == seq[j - 1] else mismatch),
                prev[j] + gap,
                cur[j - 1] + gap,
            )
            k = 0 if vals[0] >= vals[1] and vals[0] >= vals[2] else (1 if vals[1] >= vals[2] else 2)
            cur[j] = vals[k]
            drow[j] = "DUL"[k]
        trace.append(drow)
        prev = cur
    i, j = m, n
    ar: list[str] = []
    aq: list[str] = []
    while i > 0 or j > 0:
        d = trace[i - 1][j] if i > 0 and j >= 0 else "L"
        if i > 0 and j > 0 and d == "D":
            ar.append(ref[i - 1])
            aq.append(seq[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or d == "U"):
            ar.append(ref[i - 1])
            aq.append("-")
            i -= 1
        else:
            ar.append("-")
            aq.append(seq[j - 1])
            j -= 1
    return "".join(reversed(ar)), "".join(reversed(aq))


def msa_to_first(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    if not records:
        return []
    ref = records[0][1]
    max_ins = [0] * (len(ref) + 1)
    paired = []
    for name, seq in records:
        ar, aq = nw(ref, seq)
        inserts = [[] for _ in range(len(ref) + 1)]
        bases = ["-"] * len(ref)
        pos = 0
        for rb, qb in zip(ar, aq):
            if rb == "-":
                inserts[pos].append(qb)
            else:
                bases[pos] = qb
                pos += 1
        for i, ins in enumerate(inserts):
            max_ins[i] = max(max_ins[i], len(ins))
        paired.append((name, inserts, bases))
    out = []
    for name, inserts, bases in paired:
        seq_out: list[str] = []
        for i in range(len(ref) + 1):
            ins = inserts[i]
            seq_out.extend(ins + ["-"] * (max_ins[i] - len(ins)))
            if i < len(ref):
                seq_out.append(bases[i])
        out.append((name, "".join(seq_out)))
    return out


def variant_columns(aln: list[tuple[str, str]]) -> list[int]:
    if not aln:
        return []
    length = len(aln[0][1])
    cols = []
    for i in range(length):
        states = {seq[i] for _, seq in aln if seq[i] not in {"N", "?"}}
        if len(states) > 1:
            cols.append(i)
    return cols


def display_window_sequence(hit: AmpliconHit, display_flank: int | None = None) -> str:
    """Return the product plus flanks in the same F-to-R orientation as product_seq."""
    seq = revcomp(hit.window_seq) if hit.strand == "-" else hit.window_seq
    if display_flank is None:
        return seq
    available_flank = max(0, (len(seq) - hit.product_len) // 2)
    keep_flank = max(0, min(display_flank, available_flank))
    start = available_flank - keep_flank
    end = available_flank + hit.product_len + keep_flank
    return seq[start:end]


def display_window_coord(hit: AmpliconHit, offset: int, genome_len: int, display_flank: int | None = None) -> int:
    if display_flank is not None:
        available_flank = max(0, (len(display_window_sequence(hit)) - hit.product_len) // 2)
        # Recalculate using the full stored window so the display offset maps back
        # to the original genome coordinate.
        full_flank = max(0, (len(hit.window_seq) - hit.product_len) // 2)
        start_shift = full_flank - available_flank
        offset += start_shift
    if hit.strand == "-":
        return coord_add(hit.window_end, -offset, genome_len)
    return coord_add(hit.window_start, offset, genome_len)


def build_alignment_position_labels(aln: list[tuple[str, str]], ref_hit: AmpliconHit | None, genome_len: int, display_flank: int | None = None) -> dict[int, str]:
    labels: dict[int, str] = {}
    if not aln or not ref_hit or not genome_len:
        return labels
    offset = -1
    last_coord = ref_hit.window_start
    for col, base in enumerate(aln[0][1]):
        if base != "-":
            offset += 1
            last_coord = display_window_coord(ref_hit, offset, genome_len, display_flank)
            labels[col] = str(last_coord)
        else:
            labels[col] = f"gap_after_{last_coord}"
    return labels


def build_alignment_primer_marks(
    aln: list[tuple[str, str]],
    selected_hits: dict[str, AmpliconHit],
    primer: PrimerPair,
    display_flank: int | None = None,
) -> dict[tuple[int, int], str]:
    marks: dict[tuple[int, int], str] = {}
    for row, (sample, seq) in enumerate(aln):
        hit = selected_hits.get(sample)
        if not hit:
            continue
        flank = max(0, (len(display_window_sequence(hit, display_flank)) - hit.product_len) // 2)
        f_range = range(flank, flank + len(primer.f_seq))
        r_start = flank + hit.product_len - len(primer.r_seq)
        r_range = range(r_start, r_start + len(primer.r_seq))
        offset = -1
        for col, base in enumerate(seq):
            if base == "-":
                continue
            offset += 1
            if offset in f_range:
                marks[(row, col)] = "F"
            elif offset in r_range:
                marks[(row, col)] = "R"
    return marks


def alignment_pdf_svg(
    aln: list[tuple[str, str]],
    title: str,
    colors: dict[str, str],
    position_labels: dict[int, str] | None = None,
    primer_marks: dict[tuple[int, int], str] | None = None,
    chunk_size: int = 90,
) -> tuple[str, list[dict]]:
    var_cols = variant_columns(aln)
    if not aln:
        return "", []
    length = len(aln[0][1])
    windows = [(s, min(s + chunk_size, length)) for s in range(0, length, chunk_size)]

    row_h = 19
    cell = 12
    left = 155
    top = 55
    block_gap = 44
    max_wcols = max(e - s for s, e in windows)
    width = max(980, left + max_wcols * cell + 80)
    height = top + sum((len(aln) * row_h + 66) for _ in windows) + block_gap * len(windows) + 42
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='18' y='28' font-family='Arial, Microsoft YaHei' font-size='18' font-weight='700'>{html.escape(title)}</text>",
        "<text x='18' y='46' font-family='Arial, Microsoft YaHei' font-size='11' fill='#475467'>Complete amplicon window is shown in chunks: left flank + F primer + PCR product + R primer + right flank. Green boxes = F primer; red boxes = R primer; red outlines = variant columns.</text>",
    ]
    y = top
    variant_rows: list[dict] = []
    seen_variants: set[int] = set()
    var_set = set(var_cols)
    for wi, (s, e) in enumerate(windows, 1):
        parts.append(f"<text x='18' y='{y-12}' font-family='Arial' font-size='12'>Block {wi}: alignment columns {s+1}-{e}</text>")
        for r, (name, seq) in enumerate(aln):
            yy = y + r * row_h
            parts.append(f"<text x='12' y='{yy+14}' font-family='Arial' font-size='12'>{html.escape(name)}</text>")
            for j, pos in enumerate(range(s, e)):
                b = seq[pos]
                x = left + j * cell
                fill = colors.get(b, colors.get("-", "#d9d9d9"))
                text_color = "#ffffff" if b == "G" and fill.lower() in {"#111111", "black"} else "#111111"
                parts.append(f"<rect x='{x}' y='{yy}' width='{cell}' height='{row_h-2}' fill='{fill}' stroke='white'/>")
                mark = (primer_marks or {}).get((r, pos))
                if mark:
                    stroke = "#16a34a" if mark == "F" else "#dc2626"
                    parts.append(f"<rect x='{x+1}' y='{yy+1}' width='{cell-2}' height='{row_h-4}' fill='none' stroke='{stroke}' stroke-width='2'/>")
                parts.append(f"<text x='{x+cell/2}' y='{yy+13}' text-anchor='middle' font-family='Arial' font-size='10' font-weight='700' fill='{text_color}'>{html.escape(b)}</text>")
        for pos in var_cols:
            if s <= pos < e:
                x = left + (pos - s) * cell
                parts.append(f"<rect x='{x}' y='{y-2}' width='{cell}' height='{len(aln)*row_h+2}' fill='none' stroke='#e11d48' stroke-width='2'/>")
                states = "".join(sorted({seq[pos] for _, seq in aln}))
                typ = "InDel" if "-" in states else "SNV"
                label = (position_labels or {}).get(pos, str(pos + 1))
                if pos not in seen_variants:
                    variant_rows.append({"alignment_col": pos + 1, "position_label": label, "type": typ, "states": states})
                    seen_variants.add(pos)
        label_y = y + len(aln) * row_h + 18
        last_label_x = -999
        for pos in range(s, e):
            x = left + (pos - s) * cell + cell / 2
            should_label = pos in {s, e - 1} or pos in var_set or (pos - s) % 20 == 0
            if not should_label or x - last_label_x < 46:
                continue
            label = (position_labels or {}).get(pos, str(pos + 1))
            parts.append(f"<line x1='{x}' y1='{y+len(aln)*row_h+3}' x2='{x}' y2='{y+len(aln)*row_h+8}' stroke='#667085'/>")
            parts.append(f"<text x='{x}' y='{label_y}' text-anchor='end' transform='rotate(-45 {x} {label_y})' font-family='Arial' font-size='9' fill='#475467'>{html.escape(label)}</text>")
            last_label_x = x
        y += len(aln) * row_h + 66 + block_gap
    parts.append("</svg>")
    return "\n".join(parts), variant_rows


def stage2_alignment_reports(
    root: Path,
    primers: list[PrimerPair],
    records: list[FastaRecord],
    hits: dict[tuple[str, str], list[AmpliconHit]],
    align_samples: str,
    colors: dict[str, str],
    align_flank: int | None = None,
):
    stage = root / "02_multi_sequence_alignment"
    pdf_dir = stage / "pdf"
    html_dir = stage / "html"
    fasta_dir = stage / "alignment_fasta"
    table_dir = stage / "variant_tables"
    if align_samples.lower() == "all":
        sample_set = {r.sample for r in records}
    else:
        sample_set = {s.strip() for s in re.split(r"[,;，；\s]+", align_samples) if s.strip()}
    genome_len = {r.sample: len(r.seq) for r in records}
    summary: list[dict] = []

    for primer in primers:
        seqs: list[tuple[str, str]] = []
        selected_hits: dict[str, AmpliconHit] = {}
        missing: list[str] = []
        for rec in records:
            if rec.sample not in sample_set:
                continue
            candidates = [h for h in hits.get((rec.sample, primer.pair_name), []) if h.product_seq and h.status in {"FOUND", "WARN_MULTIPLE_HITS"}]
            if candidates:
                seqs.append((rec.sample, display_window_sequence(candidates[0], align_flank)))
                selected_hits[rec.sample] = candidates[0]
            else:
                missing.append(rec.sample)
        aln = msa_to_first(seqs)
        fasta_path = fasta_dir / f"{safe_name(primer.pair_name)}.alignment.fasta"
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        with fasta_path.open("w", encoding="ascii") as f:
            for name, seq in aln:
                f.write(f">{name}\n{textwrap.fill(seq, 80)}\n")
        position_labels: dict[int, str] = {}
        if aln:
            ref_sample = aln[0][0]
            ref_hit = selected_hits.get(ref_sample)
            ref_len = genome_len.get(ref_sample, 0)
            position_labels = build_alignment_position_labels(aln, ref_hit, ref_len, align_flank)
        primer_marks = build_alignment_primer_marks(aln, selected_hits, primer, align_flank)
        svg, variant_rows = alignment_pdf_svg(aln, f"{primer.order}. {primer.pair_name} 多序列比对", colors, position_labels, primer_marks)
        variant_path = table_dir / f"{safe_name(primer.pair_name)}.variant_sites.csv"
        write_csv(variant_path, variant_rows, ["alignment_col", "position_label", "type", "states"])
        miss_html = "".join(f"<li>{html.escape(m)}: not found</li>" for m in missing)
        body = (
            f"<h1>{html.escape(primer.pair_name)} 多序列比对</h1>"
            f"<p>Samples: {len(seqs)} aligned; missing: {len(missing)}</p>"
            f"<p>Displayed region: theoretical amplicon plus {align_flank if align_flank is not None else 'all available'} bp flanking sequence on both sides. "
            f"F primer is boxed in green; R primer is boxed in red.</p>"
            f"<div>{svg}</div>"
            f"<h2>Missing samples</h2><ul>{miss_html}</ul>"
        )
        html_path = html_dir / f"{safe_name(primer.pair_name)}.alignment.html"
        pdf_path = pdf_dir / f"{safe_name(primer.pair_name)}.alignment.pdf"
        html_path.parent.mkdir(parents=True, exist_ok=True)
        html_path.write_text(html_doc(f"{primer.pair_name} alignment", body, landscape=True), encoding="utf-8")
        edge_print_pdf(html_path, pdf_path)
        summary.append(
            {
                "primer_order": primer.order,
                "primer_pair": primer.pair_name,
                "aligned_samples": len(seqs),
                "missing_samples": ";".join(missing),
                "variant_sites": len(variant_rows),
                "pdf": str(pdf_path),
                "alignment_fasta": str(fasta_path),
                "variant_table": str(variant_path),
            }
        )
    write_csv(stage / "multi_sequence_alignment_summary.csv", summary, ["primer_order", "primer_pair", "aligned_samples", "missing_samples", "variant_sites", "pdf", "alignment_fasta", "variant_table"])


def parse_colors(args) -> dict[str, str]:
    colors = dict(BASE_COLORS)
    for key, value in [("A", args.color_a), ("T", args.color_t), ("C", args.color_c), ("G", args.color_g), ("-", args.color_gap)]:
        if value:
            colors[key] = value
    return colors


def input_index(root: Path, primers: list[PrimerPair], records: list[FastaRecord]):
    idx = root / "00_input_index"
    write_csv(
        idx / "primer_table_parsed.csv",
        [
            {
                "order": p.order,
                "primer_pair": p.pair_name,
                "forward_primer_name": p.f_name,
                "forward_primer_seq": p.f_seq,
                "reverse_primer_name": p.r_name,
                "reverse_primer_seq": p.r_seq,
            }
            for p in primers
        ],
        ["order", "primer_pair", "forward_primer_name", "forward_primer_seq", "reverse_primer_name", "reverse_primer_seq"],
    )
    write_csv(
        idx / "sample_fasta_index.csv",
        [{"sample": r.sample, "path": str(r.path), "header": r.header, "length": len(r.seq)} for r in records],
        ["sample", "path", "header", "length"],
    )


def zip_pdfs(root: Path):
    zip_dir = root / "04_pdf_only_zip"
    zip_dir.mkdir(parents=True, exist_ok=True)
    zip_path = zip_dir / "全部PDF.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for pdf in root.glob("**/*.pdf"):
            if "04_pdf_only_zip" in pdf.parts:
                continue
            z.write(pdf, pdf.relative_to(root))
    return zip_path


def classify_sanger_files(sanger_dir: Path, out_dir: Path):
    seq_ext = {".fasta", ".fa", ".fas", ".seq", ".txt"}
    peak_ext = {".ab1", ".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff"}
    rows = []
    for p in sorted(sanger_dir.rglob("*")):
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        kind = "sequence" if ext in seq_ext else ("peak_or_support" if ext in peak_ext else "other")
        direction = ""
        if re.search(r"(^|[_\-\[\(])F([_\-\]\)]|$)", p.stem, re.I):
            direction = "F"
        elif re.search(r"(^|[_\-\[\(])R([_\-\]\)]|$)", p.stem, re.I):
            direction = "R"
        rows.append({"file": p.name, "path": str(p), "kind": kind, "direction_guess": direction})
    write_csv(out_dir / "sanger_file_index.csv", rows, ["file", "path", "kind", "direction_guess"])


def parse_number_map(text: str | None) -> dict[str, str]:
    if not text:
        return {}
    out: dict[str, str] = {}
    for part in re.split(r"[,;，；\s]+", text.strip()):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif ":" in part:
            k, v = part.split(":", 1)
        else:
            continue
        out[k.strip()] = v.strip()
    return out


def primer_key(pair_name: str) -> str:
    m = re.match(r"(G\d+|G\d+\d*)", pair_name)
    if m:
        return m.group(1)
    return pair_name.split("-")[0].strip()


def parse_sanger_filename(path: Path) -> dict | None:
    name = path.name
    m = re.search(r"_\((G\d+)-([^)]+)\)_\[(G\d+)-([FR])\]", name, re.I)
    if not m:
        return None
    return {
        "group": m.group(1),
        "number": m.group(2),
        "tag_group": m.group(3),
        "direction": m.group(4).upper(),
    }


def read_sequence_file(path: Path) -> str:
    seq = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        seq.append(clean_seq(line))
    return "".join(seq)


def sanger_metrics(theory: str, read: str) -> tuple[str, str, dict, list[dict]]:
    at, ar = nw(theory, read)
    tpos = rpos = 0
    matches = mismatches = gaps = insertions = covered = 0
    first = last = ""
    diffs = []
    for col, (tb, rb) in enumerate(zip(at, ar), 1):
        t_label = ""
        r_label = ""
        if tb != "-":
            tpos += 1
            t_label = tpos
        if rb != "-":
            rpos += 1
            r_label = rpos
        if tb != "-" and rb != "-":
            covered += 1
            if first == "":
                first = tpos
            last = tpos
            if tb == rb:
                matches += 1
            else:
                mismatches += 1
                diffs.append({"alignment_col": col, "theory_pos": t_label, "read_pos": r_label, "type": "Mismatch", "theory": tb, "read": rb})
        elif tb != "-" and rb == "-":
            gaps += 1
            if first == "":
                first = tpos
            last = tpos
            diffs.append({"alignment_col": col, "theory_pos": t_label, "read_pos": "", "type": "Gap_in_Read", "theory": tb, "read": "-"})
        elif tb == "-" and rb != "-":
            insertions += 1
            diffs.append({"alignment_col": col, "theory_pos": f"gap_after_{tpos}", "read_pos": r_label, "type": "Insertion_vs_Theory", "theory": "-", "read": rb})
    identity = round(matches / max(1, covered) * 100, 2)
    coverage = round(covered / max(1, len(theory)) * 100, 2)
    stats = {
        "theory_length": len(theory),
        "read_length": len(read),
        "covered_theory_bases": covered,
        "coverage_percent": coverage,
        "identity_percent": identity,
        "matches": matches,
        "mismatches": mismatches,
        "gaps_in_read": gaps,
        "insertions_vs_theory": insertions,
        "first_covered_theory_pos": first,
        "last_covered_theory_pos": last,
    }
    return at, ar, stats, diffs


def choose_read_orientation(theory: str, seq: str, nominal_direction: str) -> tuple[str, str, str, str, dict, list[dict]]:
    candidates = [("forward", seq), ("revcomp", revcomp(seq))]
    best = None
    for label, oriented in candidates:
        at, ar, stats, diffs = sanger_metrics(theory, oriented)
        rank = (stats["identity_percent"], stats["coverage_percent"], -stats["gaps_in_read"] - stats["insertions_vs_theory"])
        if best is None or rank > best[0]:
            best = (rank, label, oriented, at, ar, stats, diffs)
    assert best is not None
    _, label, oriented, at, ar, stats, diffs = best
    return label, oriented, at, ar, stats, diffs


def read_overlay_svg(hit: AmpliconHit, primer: PrimerPair, read_stats: list[dict], width: int = 900) -> str:
    flank = max(0, (len(hit.window_seq) - hit.product_len) // 2)
    total = len(hit.window_seq)
    left, right = 55, 45
    track_w = width - left - right

    def x(offset: float) -> float:
        return left + (offset / max(1, total)) * track_w

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{120 + 24*len(read_stats)}' viewBox='0 0 {width} {120 + 24*len(read_stats)}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='12' y='20' font-family='Arial' font-size='12'>window {hit.window_start}-{hit.window_end}; PCR {hit.product_start}-{hit.product_end}</text>",
        f"<line x1='{left}' y1='52' x2='{width-right}' y2='52' stroke='#3366aa' stroke-width='2'/>",
        f"<rect x='{x(flank)}' y='40' width='{max(2, x(flank + hit.product_len)-x(flank))}' height='24' fill='#e8f0fe' stroke='#7aa5ff'/>",
    ]
    f_offset = flank
    r_offset = flank + hit.product_len - len(primer.r_seq)
    parts.append(f"<rect x='{x(f_offset)}' y='32' width='{max(3, x(f_offset+len(primer.f_seq))-x(f_offset))}' height='40' fill='#18a058' opacity='0.85'/>")
    parts.append(f"<rect x='{x(r_offset)}' y='32' width='{max(3, x(r_offset+len(primer.r_seq))-x(r_offset))}' height='40' fill='#d83b36' opacity='0.85'/>")
    y = 92
    for rs in read_stats:
        first = rs.get("first_covered_theory_pos") or 1
        last = rs.get("last_covered_theory_pos") or 1
        sx = x(flank + int(first) - 1)
        ex = x(flank + int(last))
        color = "#7c3aed" if rs.get("direction") == "R" else "#0f766e"
        parts.append(f"<text x='12' y='{y+12}' font-family='Arial' font-size='11'>{html.escape(rs.get('direction',''))} read</text>")
        parts.append(f"<rect x='{sx}' y='{y}' width='{max(2, ex-sx)}' height='14' fill='{color}' opacity='0.75'/>")
        parts.append(f"<text x='{min(width-180, ex+6)}' y='{y+12}' font-family='Arial' font-size='11'>cov {rs.get('coverage_percent')}%, id {rs.get('identity_percent')}%</text>")
        y += 24
    parts.append("</svg>")
    return "\n".join(parts)


def pairwise_alignment_svg(
    theory_aln: str,
    read_aln: str,
    title: str,
    coord_labels: list[str],
    colors: dict[str, str] | None = None,
    primer_marks: dict[int, str] | None = None,
    block: int = 80,
    width: int = 1060,
) -> str:
    colors = colors or BASE_COLORS
    cell = 12
    left = 72
    top = 42
    row_h = 26
    block_gap = 34
    nblocks = max(1, (len(theory_aln) + block - 1) // block)
    height = top + nblocks * (row_h * 3 + block_gap) + 30
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='12' y='24' font-family='Arial, Microsoft YaHei' font-size='16' font-weight='700'>{html.escape(title)}</text>",
    ]
    for bi, start in enumerate(range(0, len(theory_aln), block), 1):
        end = min(len(theory_aln), start + block)
        y = top + (bi - 1) * (row_h * 3 + block_gap)
        label_start = coord_labels[start] if start < len(coord_labels) else ""
        label_end = coord_labels[end - 1] if end - 1 < len(coord_labels) else ""
        parts.append(f"<text x='{left}' y='{y-8}' font-family='Arial' font-size='11'>Genome position: {html.escape(str(label_start))} - {html.escape(str(label_end))}</text>")
        rows = [("Theory", theory_aln[start:end]), ("Sanger", read_aln[start:end]), ("Diff", "".join("|" if a != b else " " for a, b in zip(theory_aln[start:end], read_aln[start:end])))]
        for ri, (label, seq) in enumerate(rows):
            yy = y + ri * row_h
            parts.append(f"<text x='10' y='{yy+16}' font-family='Arial' font-size='14'>{label}</text>")
            for j, b in enumerate(seq):
                x = left + j * cell
                if label == "Diff":
                    if b.strip():
                        parts.append(f"<text x='{x+cell/2}' y='{yy+17}' text-anchor='middle' font-family='Arial' font-size='16' fill='#cc1f1a'>|</text>")
                    continue
                fill = colors.get(b, colors.get("-", "#d9d9d9"))
                mark = (primer_marks or {}).get(start + j, "")
                if mark == "F":
                    fill = "#86efac"
                    stroke = "#078a3f"
                    stroke_w = 1
                elif mark == "R":
                    fill = "#fca5a5"
                    stroke = "#d11f1f"
                    stroke_w = 1
                else:
                    stroke = "white"
                    stroke_w = 1
                text_color = "#ffffff" if b == "G" and fill.lower() in {"#111111", "black"} else "#111111"
                opacity = "0.92" if mark else "0.28"
                parts.append(f"<rect x='{x}' y='{yy}' width='{cell}' height='{row_h-4}' fill='{fill}' opacity='{opacity}' stroke='{stroke}' stroke-width='{stroke_w}'/>")
                parts.append(f"<text x='{x+cell/2}' y='{yy+16}' text-anchor='middle' font-family='Consolas' font-size='13' fill='{text_color}'>{html.escape(b)}</text>")
    parts.append("</svg>")
    return "\n".join(parts)


def alignment_coord_labels(theory_aln: str, hit: AmpliconHit, genome_len: int) -> list[str]:
    labels: list[str] = []
    tpos = 0
    last_coord = hit.product_start
    for b in theory_aln:
        if b != "-":
            tpos += 1
            last_coord = coord_add(hit.product_start, tpos - 1, genome_len)
            labels.append(str(last_coord))
        else:
            labels.append(f"gap_after_{last_coord}")
    return labels


def window_alignment_for_sanger(product_theory_aln: str, read_aln: str, hit: AmpliconHit, genome_len: int) -> tuple[str, str, list[str]]:
    """Embed product-level Sanger alignment inside the ±flank theory window.

    Statistics are still calculated on the PCR product. Visualization shows the
    full theory window (product ±500 bp by default) and every real Sanger base
    available in the product alignment. Uncovered flank/product bases are gaps
    in the Sanger row.
    """
    flank = max(0, (len(hit.window_seq) - hit.product_len) // 2)
    left = hit.window_seq[:flank]
    right = hit.window_seq[flank + hit.product_len :]
    theory_display = left + product_theory_aln + right
    read_display = ("-" * len(left)) + read_aln + ("-" * len(right))
    labels: list[str] = []
    for i in range(len(left)):
        labels.append(str(coord_add(hit.window_start, i, genome_len)))
    labels.extend(alignment_coord_labels(product_theory_aln, hit, genome_len))
    right_start = coord_add(hit.product_start, hit.product_len, genome_len)
    for i in range(len(right)):
        labels.append(str(coord_add(right_start, i, genome_len)))
    return theory_display, read_display, labels


def primer_marks_for_display(display_theory: str, hit: AmpliconHit, primer: PrimerPair) -> dict[int, str]:
    """Return alignment-column marks for F/R primer regions in displayed theory."""
    flank = max(0, (len(hit.window_seq) - hit.product_len) // 2)
    right_len = max(0, len(hit.window_seq) - flank - hit.product_len)
    product_aln_end = len(display_theory) - right_len
    marks: dict[int, str] = {}
    product_pos = 0
    for col in range(flank, product_aln_end):
        if display_theory[col] == "-":
            continue
        product_pos += 1
        if 1 <= product_pos <= len(primer.f_seq):
            marks[col] = "F"
        elif hit.product_len - len(primer.r_seq) < product_pos <= hit.product_len:
            marks[col] = "R"
    return marks


def stage3_sanger_compare(
    root: Path,
    primers: list[PrimerPair],
    records: list[FastaRecord],
    hits: dict[tuple[str, str], list[AmpliconHit]],
    sanger_dir: Path,
    number_map: dict[str, str],
):
    stage = root / "03_sanger_vs_theory"
    pdf_dir = stage / "pdf"
    html_dir = stage / "html"
    table_dir = stage / "tables"
    classify_sanger_files(sanger_dir, table_dir)

    sample_lookup = {r.sample: r for r in records}
    primer_by_key = {primer_key(p.pair_name): p for p in primers}
    seq_files: dict[tuple[str, str, str], list[Path]] = {}
    support_files: dict[tuple[str, str, str], list[Path]] = {}
    observed_keys: set[tuple[str, str]] = set()
    for p in sorted(sanger_dir.rglob("*")):
        if not p.is_file():
            continue
        meta = parse_sanger_filename(p)
        if not meta:
            continue
        sample = number_map.get(meta["number"], meta["number"])
        key = (meta["group"], sample, meta["direction"])
        observed_keys.add((meta["group"], sample))
        if p.suffix.lower() in {".fasta", ".fa", ".fas", ".seq", ".txt"}:
            seq_files.setdefault(key, []).append(p)
        else:
            support_files.setdefault(key, []).append(p)

    summary_rows: list[dict] = []
    missing_rows: list[dict] = []
    for pkey, primer in primer_by_key.items():
        body = [f"<h1>{html.escape(primer.pair_name)} 胶回收测序结果 vs 理论扩增片段</h1>"]
        body.append("<p>仅展示真实存在的胶回收测序序列文件；未上传或测序失败的样本只记录在缺失表中。每个read按Theory/Sanger/Diff三行显示。</p>")
        samples_for_primer = sorted({sample for group, sample in observed_keys if group == pkey}, key=natural_key)
        if not samples_for_primer:
            body.append("<p class='warn'>没有识别到该引物对的胶回收测序序列文件。</p>")
        for sample in samples_for_primer:
            rec = sample_lookup.get(sample)
            if not rec:
                missing_rows.append({"primer_pair": primer.pair_name, "sample": sample, "direction": "sample", "status": "sample_genome_not_uploaded"})
                continue
            sample = rec.sample
            candidates = [h for h in hits.get((sample, primer.pair_name), []) if h.product_seq and h.status in {"FOUND", "WARN_MULTIPLE_HITS"}]
            if not candidates:
                missing_rows.append({"primer_pair": primer.pair_name, "sample": sample, "direction": "theory", "status": "theory_not_found"})
                continue
            hit = candidates[0]
            read_blocks: list[str] = []
            stats_rows: list[dict] = []
            for direction in ["F", "R"]:
                files = seq_files.get((pkey, sample, direction), [])
                if not files:
                    missing_rows.append({"primer_pair": primer.pair_name, "sample": sample, "direction": direction, "status": "sequence_file_missing"})
                    continue
                # Prefer FASTA over .seq when both exist for the same read.
                files = sorted(files, key=lambda x: 0 if x.suffix.lower() in {".fasta", ".fa", ".fas"} else 1)
                seq = read_sequence_file(files[0])
                if not seq:
                    missing_rows.append({"primer_pair": primer.pair_name, "sample": sample, "direction": direction, "status": "empty_sequence"})
                    continue
                orientation, oriented, at, ar, stats, diffs = choose_read_orientation(hit.product_seq, seq, direction)
                stats.update({"direction": direction, "orientation": orientation, "sequence_file": files[0].name})
                stats_rows.append(stats)
                display_theory, display_read, coords = window_alignment_for_sanger(at, ar, hit, len(rec.seq))
                primer_marks = primer_marks_for_display(display_theory, hit, primer)
                read_title = f"Sample {sample} {direction} read: {files[0].name}; orientation {orientation}; coverage {stats['coverage_percent']}%; identity {stats['identity_percent']}%"
                read_blocks.append(pairwise_alignment_svg(display_theory, display_read, read_title, coords, primer_marks=primer_marks))
                summary = {
                    "primer_pair": primer.pair_name,
                    "sample": sample,
                    "direction": direction,
                    "sequence_file": files[0].name,
                    "orientation_used": orientation,
                    **stats,
                    "support_files": ";".join(x.name for x in support_files.get((pkey, sample, direction), [])),
                }
                summary_rows.append(summary)
            if not read_blocks:
                continue
            body.append(f"<div class='page'><h2>Sample {html.escape(sample)}</h2>")
            body.append(region_track_svg(hit, primer))
            body.append("<table><thead><tr><th>Direction</th><th>File</th><th>Orientation</th><th>Coverage %</th><th>Identity %</th><th>Mismatches</th><th>Gaps</th><th>Insertions</th><th>Support files</th></tr></thead><tbody>")
            for rs in stats_rows:
                supports = ";".join(x.name for x in support_files.get((pkey, sample, rs["direction"]), []))
                body.append(
                    f"<tr><td>{rs['direction']}</td><td>{html.escape(rs['sequence_file'])}</td><td>{rs['orientation']}</td>"
                    f"<td>{rs['coverage_percent']}</td><td>{rs['identity_percent']}</td><td>{rs['mismatches']}</td>"
                    f"<td>{rs['gaps_in_read']}</td><td>{rs['insertions_vs_theory']}</td><td>{html.escape(supports)}</td></tr>"
                )
            body.append("</tbody></table>")
            body.extend(read_blocks)
            body.append("</div>")
        html_path = html_dir / f"{safe_name(primer.pair_name)}.sanger_vs_theory.html"
        pdf_path = pdf_dir / f"{safe_name(primer.pair_name)}.sanger_vs_theory.pdf"
        html_path.parent.mkdir(parents=True, exist_ok=True)
        html_path.write_text(html_doc(f"{primer.pair_name} sanger vs theory", "\n".join(body), landscape=True), encoding="utf-8")
        edge_print_pdf(html_path, pdf_path)

    if summary_rows:
        fields = list(summary_rows[0].keys())
        write_csv(table_dir / "sanger_vs_theory_summary.csv", summary_rows, fields)
    write_csv(table_dir / "sanger_missing_or_failed.csv", missing_rows, ["primer_pair", "sample", "direction", "status"])
    return summary_rows, missing_rows


def create_stage3_word_report(root: Path, summary_rows: list[dict], missing_rows: list[dict]):
    report_dir = root / "03_sanger_vs_theory" / "word_report"
    by_primer: dict[str, list[dict]] = {}
    for row in summary_rows:
        by_primer.setdefault(row.get("primer_pair", ""), []).append(row)
    overview = [
        f"Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.",
        f"Successful sequence-to-theory comparisons: {len(summary_rows)}.",
        f"Missing or failed records: {len(missing_rows)}.",
        "Only real uploaded sequencing sequence files are counted as comparisons. Samples without sequence files are listed as missing/failed, not treated as compared.",
        "The PDF visualizations remain the primary visual evidence; this Word report summarizes all real comparison statistics and flags.",
    ]
    primer_table = [["Primer pair", "Comparisons", "Mean coverage %", "Mean identity %", "Mismatches", "Gaps", "Insertions"]]
    notable = [["Primer pair", "Sample", "Direction", "Coverage %", "Identity %", "Issue"]]
    for primer, rows in sorted(by_primer.items(), key=lambda kv: natural_key(kv[0])):
        covs = [float(r.get("coverage_percent", 0) or 0) for r in rows]
        ids = [float(r.get("identity_percent", 0) or 0) for r in rows]
        mismatches = sum(int(r.get("mismatches", 0) or 0) for r in rows)
        gaps = sum(int(r.get("gaps_in_read", 0) or 0) for r in rows)
        insertions = sum(int(r.get("insertions_vs_theory", 0) or 0) for r in rows)
        primer_table.append([
            primer,
            str(len(rows)),
            f"{sum(covs)/len(covs):.2f}" if covs else "0",
            f"{sum(ids)/len(ids):.2f}" if ids else "0",
            str(mismatches),
            str(gaps),
            str(insertions),
        ])
        for r in rows:
            issues = []
            if float(r.get("coverage_percent", 0) or 0) < 80:
                issues.append("low coverage")
            if float(r.get("identity_percent", 0) or 0) < 95:
                issues.append("low identity")
            if int(r.get("mismatches", 0) or 0) > 10:
                issues.append("many mismatches")
            if issues:
                notable.append([primer, r.get("sample", ""), r.get("direction", ""), r.get("coverage_percent", ""), r.get("identity_percent", ""), "; ".join(issues)])
    missing_table = [["Primer pair", "Sample", "Direction", "Status"]]
    for r in missing_rows:
        missing_table.append([r.get("primer_pair", ""), r.get("sample", ""), r.get("direction", ""), r.get("status", "")])
    sections = [
        ("Overview", overview, None),
        ("Per-primer summary", ["Coverage and identity are calculated against the theoretical amplicon sequence for each sample/read."], primer_table),
        ("Notable comparisons", ["Rows below are automatically flagged by coverage, identity, or mismatch thresholds."], notable),
        ("Missing or failed files", ["Missing sequence files are treated as sequencing failure or absent uploaded files."], missing_table),
        ("Recommended follow-up", [
            "Check low-coverage and low-identity reads against the original chromatogram files.",
            "Confirm sample naming when a read is assigned to an unexpected primer pair or direction.",
            "For theory-not-found cases, inspect primer sequence, genome assembly orientation, and potential primer mismatches.",
        ], None),
    ]
    out = report_dir / "胶回收测序比对专业汇总报告.docx"
    make_docx(out, "胶回收测序结果与理论扩增片段比对专业汇总报告", sections)
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PCR amplicon multi-sequence alignment workflow")
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name in ["stage1", "align", "full", "sanger-index", "sanger-compare"]:
        p = sub.add_parser(name)
        p.add_argument("--primer-table", required=name != "sanger-index", help="Primer Excel table")
        p.add_argument("--genome-dir", help="Directory containing genome FASTA files")
        p.add_argument("--fasta", action="append", help="Genome FASTA file or directory; can repeat")
        p.add_argument("--out-dir", help="Run output directory; default creates desktop timestamp folder")
        p.add_argument("--flank", type=int, default=500)
        p.add_argument("--max-product", type=int, default=5000)
        p.add_argument("--align-flank", type=int, default=None, help="Flanking bp to show on each side in multi-sequence alignment; default uses --flank")
        p.add_argument("--align-samples", default="all", help="all, or comma-separated sample names")
        p.add_argument("--color-a", default=None)
        p.add_argument("--color-t", default=None)
        p.add_argument("--color-c", default=None)
        p.add_argument("--color-g", default=None)
        p.add_argument("--color-gap", default=None)
        p.add_argument("--zip-pdfs", action="store_true")
        p.add_argument("--sanger-dir", help="Gel-recovered sequencing result folder")
        p.add_argument("--number-map", help="Sequencing number to sample map, e.g. 1=11,2=51,3=1")
        p.add_argument("--word-report", action="store_true", help="After Sanger comparison, create a professional Word summary report")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    root = run_root(args.out_dir)
    if args.cmd == "sanger-index":
        if not args.sanger_dir:
            raise SystemExit("--sanger-dir is required for sanger-index")
        stage = root / "03_sanger_vs_theory"
        classify_sanger_files(Path(args.sanger_dir), stage)
        print(root)
        return 0

    primers = read_primers(Path(args.primer_table))
    records = read_fastas(args.fasta, args.genome_dir)
    if not primers:
        raise SystemExit("No primer pairs found in the primer table.")
    if not records:
        raise SystemExit("No genome FASTA records found.")
    input_index(root, primers, records)
    hits = stage1_sample_reports(root, primers, records, args.flank, args.max_product)
    if args.cmd in {"align", "full"}:
        stage2_alignment_reports(root, primers, records, hits, args.align_samples, parse_colors(args), args.align_flank if args.align_flank is not None else args.flank)
    if args.cmd in {"sanger-compare", "full"} and args.sanger_dir:
        sanger_summary, sanger_missing = stage3_sanger_compare(root, primers, records, hits, Path(args.sanger_dir), parse_number_map(args.number_map))
        if args.word_report:
            docx = create_stage3_word_report(root, sanger_summary, sanger_missing)
            print(f"Word report: {docx}")
    elif args.sanger_dir:
        classify_sanger_files(Path(args.sanger_dir), root / "03_sanger_vs_theory")
    if args.zip_pdfs:
        zip_path = zip_pdfs(root)
        print(f"PDF zip: {zip_path}")
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
