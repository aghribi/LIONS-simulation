"""PDF of the LIONS guide: docs/guide/*.md in filename order, rendered with
fpdf2.

Usage: uv run --with fpdf2 python generate_pdf.py
Output: docs/_build/pdf/LIONS.pdf
"""

import datetime
import glob
import os
import re

from fpdf import FPDF
from fpdf.enums import XPos, YPos

DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
GUIDE_DIR = os.path.join(DOCS_DIR, 'guide')
PDF_FILE = os.path.join(DOCS_DIR, '_build', 'pdf', 'LIONS.pdf')
TITLE = 'LIONS technical guide'


class DocPDF(FPDF):
    """PDF with page header and footer."""

    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(120, 120, 120)
            self.cell(0, 5, TITLE, align='C')
            self.ln(8)

    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', align='C')


def clean_text(text):
    """Strip inline markdown (code, bold, italics, links)."""
    text = re.sub(r'`([^`]+)`', r'\1', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'\1', text)
    text = re.sub(r'\*([^*]+)\*', r'\1', text)
    text = re.sub(r'\[([^\]]+)\]\([^)]+\)', r'\1', text)
    return text.strip()


def parse_table(lines, start_idx):
    """Parse a markdown table; return (rows, index of the next line)."""
    rows = []
    idx = start_idx
    while idx < len(lines):
        line = lines[idx].strip()
        if not line.startswith('|'):
            break
        if re.match(r'^\|[\s\-:|]+\|$', line):
            idx += 1
            continue
        cells = [c.strip() for c in line.split('|')[1:-1]]
        rows.append(cells)
        idx += 1
    return rows, idx


def cell_ln(pdf, w, h, txt, **kwargs):
    """cell() followed by a line break."""
    pdf.cell(w, h, txt, new_x=XPos.LMARGIN, new_y=YPos.NEXT, **kwargs)


def fit_text(pdf, txt, max_w):
    """Truncate txt to fit in max_w (mm)."""
    if pdf.get_string_width(txt) + 1 <= max_w:
        return txt
    while len(txt) > 0 and pdf.get_string_width(txt + '..') + 1 > max_w:
        txt = txt[:-1]
    return txt + '..' if txt else ''


def read_guide():
    """Lines of all guide pages, in filename order."""
    pages = sorted(glob.glob(os.path.join(GUIDE_DIR, '*.md')))
    if not pages:
        raise FileNotFoundError(f'no guide page in {GUIDE_DIR}')
    lines = []
    for page in pages:
        with open(page, encoding='utf-8') as f:
            lines += f.read().split('\n') + ['']
    return pages, lines


def generate_pdf():
    pages, lines = read_guide()

    pdf = DocPDF()
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=20)

    # Title page
    pdf.add_page()
    pdf.set_font('Helvetica', 'B', 28)
    pdf.ln(50)
    cell_ln(pdf, 0, 15, 'LIONS', align='C')
    pdf.set_font('Helvetica', '', 14)
    pdf.ln(5)
    cell_ln(pdf, 0, 10, 'Technical guide', align='C')
    pdf.ln(10)
    pdf.set_font('Helvetica', '', 11)
    cell_ln(
        pdf, 0, 8, 'Ion tracking in separated-sector cyclotrons', align='C'
    )
    pdf.ln(20)
    pdf.set_font('Helvetica', 'I', 10)
    cell_ln(pdf, 0, 8, 'Fortran code: P. Bertrand, GANIL', align='C')
    cell_ln(pdf, 0, 8, 'Python port and GUI: O. Kamalou, GANIL', align='C')
    pdf.ln(30)
    pdf.set_font('Helvetica', '', 9)
    cell_ln(
        pdf,
        0,
        6,
        f'Generated on {datetime.date.today().isoformat()}',
        align='C',
    )

    # Content
    i = 0
    in_code_block = False
    code_buffer = []

    while i < len(lines):
        line = lines[i]

        # Code blocks
        if line.strip().startswith('```'):
            if in_code_block:
                in_code_block = False
                pdf.set_font('Courier', '', 7)
                pdf.set_fill_color(240, 240, 240)
                pdf.set_x(pdf.l_margin)
                code_text = '\n'.join(code_buffer)
                if code_text.strip():
                    pdf.multi_cell(0, 3.5, code_text, border=1, fill=True)
                pdf.ln(3)
                code_buffer = []
                pdf.set_font('Helvetica', '', 10)
            else:
                in_code_block = True
                code_buffer = []
            i += 1
            continue

        if in_code_block:
            code_buffer.append(line)
            i += 1
            continue

        stripped = line.strip()

        # Blank lines
        if not stripped:
            i += 1
            continue

        # Rules ---
        if re.match(r'^-{3,}$', stripped):
            pdf.ln(3)
            pdf.set_draw_color(180, 180, 180)
            pdf.line(10, pdf.get_y(), 200, pdf.get_y())
            pdf.ln(3)
            i += 1
            continue

        # Headings: H1 (one per guide page) starts a new page
        level = len(stripped) - len(stripped.lstrip('#'))
        if 1 <= level <= 4 and stripped[level : level + 1] == ' ':
            title = clean_text(stripped[level + 1 :])
            if level == 1:
                pdf.add_page()
                pdf.set_font('Helvetica', 'B', 18)
                pdf.set_text_color(0, 51, 102)
                cell_ln(pdf, 0, 12, title)
                pdf.set_draw_color(0, 51, 102)
                pdf.line(10, pdf.get_y(), 200, pdf.get_y())
                pdf.ln(5)
            elif level == 2:
                pdf.ln(5)
                pdf.set_font('Helvetica', 'B', 13)
                pdf.set_text_color(0, 70, 130)
                cell_ln(pdf, 0, 8, title)
                pdf.ln(2)
            else:
                pdf.ln(3)
                pdf.set_font('Helvetica', 'B', 11)
                pdf.set_text_color(50, 50, 50)
                cell_ln(pdf, 0, 7, title)
                pdf.ln(1)
            pdf.set_text_color(0, 0, 0)
            pdf.set_font('Helvetica', '', 10)
            i += 1
            continue

        # Tables
        if stripped.startswith('|'):
            rows, end_idx = parse_table(lines, i)
            if rows:
                n_cols = len(rows[0])
                page_w = 190
                font_size = (
                    6
                    if n_cols > 6
                    else (7 if n_cols > 5 else (8 if n_cols > 3 else 9))
                )

                # Column widths from the content
                pdf.set_font('Helvetica', 'B', font_size)
                max_lens = [0] * n_cols
                for row in rows:
                    for ci in range(min(n_cols, len(row))):
                        txt = clean_text(row[ci])
                        sw = pdf.get_string_width(txt) + 2
                        if sw > max_lens[ci]:
                            max_lens[ci] = sw
                total = sum(max_lens) if sum(max_lens) > 0 else 1
                col_widths = [max(12, ml / total * page_w) for ml in max_lens]
                # Scale so that the widths add up to page_w
                s = sum(col_widths)
                col_widths = [w / s * page_w for w in col_widths]

                # Header row
                pdf.set_font('Helvetica', 'B', font_size)
                pdf.set_fill_color(0, 51, 102)
                pdf.set_text_color(255, 255, 255)
                for ci, cell in enumerate(rows[0]):
                    w = (
                        col_widths[ci]
                        if ci < len(col_widths)
                        else col_widths[-1]
                    )
                    txt = fit_text(pdf, clean_text(cell), w)
                    pdf.cell(w, 5, txt, border=1, fill=True)
                pdf.ln()
                pdf.set_text_color(0, 0, 0)

                # Data rows
                pdf.set_font('Helvetica', '', font_size)
                for ri, row in enumerate(rows[1:]):
                    fill = ri % 2 == 0
                    if fill:
                        pdf.set_fill_color(245, 245, 250)
                    for ci in range(n_cols):
                        w = (
                            col_widths[ci]
                            if ci < len(col_widths)
                            else col_widths[-1]
                        )
                        if ci < len(row):
                            txt = fit_text(pdf, clean_text(row[ci]), w)
                        else:
                            txt = ''
                        pdf.cell(w, 5, txt, border=1, fill=fill)
                    pdf.ln()

                pdf.ln(3)
                pdf.set_font('Helvetica', '', 10)
            i = end_idx
            continue

        # Lists
        if re.match(r'^[-*]\s', stripped) or re.match(r'^\d+\.\s', stripped):
            text = re.sub(r'^[-*]\s+', '', stripped)
            text = re.sub(r'^\d+\.\s+', '', text)
            text = clean_text(text)
            indent = len(line) - len(line.lstrip())
            if indent >= 4:
                prefix = '     - '
            elif stripped[0] in '-*':
                prefix = '  - '
            else:
                num = stripped.split('.')[0].strip()
                prefix = f'  {num}. '
            pdf.set_font('Helvetica', '', 9)
            pdf.set_x(pdf.l_margin)
            full_text = prefix + text
            if full_text.strip():
                pdf.multi_cell(0, 4.5, full_text)
            pdf.set_font('Helvetica', '', 10)
            i += 1
            continue

        # Paragraphs
        text = clean_text(stripped)
        if text:
            pdf.set_font('Helvetica', '', 10)
            pdf.set_x(pdf.l_margin)
            pdf.multi_cell(0, 5, text)
            pdf.ln(2)

        i += 1

    os.makedirs(os.path.dirname(PDF_FILE), exist_ok=True)
    pdf.output(PDF_FILE)
    print(f'{len(pages)} guide pages -> {PDF_FILE}')
    print(
        f'size: {os.path.getsize(PDF_FILE) / 1024:.0f} kB, pages: {pdf.page_no()}'
    )


if __name__ == '__main__':
    generate_pdf()
