#!/usr/bin/env python3
"""Report interativo (Shiny for Python) para as variantes do painel hotgenes.

Junta SNV/SV/CNV/STR (VCF) e DMRs (BED) num unico app com tabelas
filtraveis e um IGV.js embutido, sincronizado por clique de linha.
"""
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam

from shiny import App, Inputs, Outputs, Session, module, reactive, render, ui

# --------------------------------------------------------------------------
# Funcoes principais
# --------------------------------------------------------------------------

_BED_STANDARD_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
]


def bed2dataframe(path):
    """Le um arquivo BED (com ou sem header, com ou sem linha 'track') em um DataFrame.

    Sempre adiciona as colunas internas _chrom/_start/_end (0-based, half-open),
    usadas por bed_interest_filter. Essas colunas nao aparecem como opcao de
    coluna na UI.
    """
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("track") or line.startswith("browser") or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                # decide se essa primeira linha e um header (nomes de coluna)
                # ou ja e dado (BED "puro", sem header)
                is_data = True
                try:
                    int(fields[1])
                    int(fields[2])
                except (ValueError, IndexError):
                    is_data = False
                if is_data:
                    header = [
                        _BED_STANDARD_COLS[i] if i < len(_BED_STANDARD_COLS) else f"col{i+1}"
                        for i in range(len(fields))
                    ]
                    rows.append(fields)
                else:
                    header = fields
                continue
            rows.append(fields)

    df = pd.DataFrame(rows, columns=header)
    chrom_col, start_col, end_col = header[0], header[1], header[2]
    df["_chrom"] = df[chrom_col]
    df["_start"] = df[start_col].astype(int)
    df["_end"] = df[end_col].astype(int)
    return df


def _format_value(value, is_gt=False, phased=False):
    if value is None:
        return "."
    if isinstance(value, tuple):
        if is_gt:
            sep = "|" if phased else "/"
            return sep.join("." if v is None else str(v) for v in value)
        return ",".join("." if v is None else str(v) for v in value)
    return str(value)


def _header_record_meta(record):
    meta = {}
    for attr in ("id", "number", "type", "description"):
        value = getattr(record, attr, None)
        if value not in (None, ""):
            meta[attr.capitalize()] = value
    return meta


def vcf2dataframe(path):
    """Le um VCF (bgzip+tabix ou plano) em um DataFrame "achatado".

    Cada chave de INFO e de FORMAT vira uma coluna separada. Retorna
    (dataframe, header_meta), onde header_meta = {"INFO": {id: {...}}, "FORMAT": {id: {...}}}
    com os atributos do header do VCF (Description, Number, Type, Id).
    """
    vcf = pysam.VariantFile(str(path))

    info_ids = list(vcf.header.info.keys())
    format_ids = list(vcf.header.formats.keys())
    sample_names = list(vcf.header.samples)

    header_meta = {
        "INFO": {i: _header_record_meta(vcf.header.info[i]) for i in info_ids},
        "FORMAT": {f: _header_record_meta(vcf.header.formats[f]) for f in format_ids},
    }

    rows = []
    for rec in vcf.fetch():
        row = {
            "CHROM": rec.chrom,
            "POS": rec.pos,
            "ID": rec.id if rec.id else ".",
            "REF": rec.ref,
            "ALT": ",".join(rec.alts) if rec.alts else ".",
            "QUAL": rec.qual if rec.qual is not None else ".",
            "FILTER": ",".join(rec.filter.keys()) if list(rec.filter.keys()) else ".",
        }
        for key in info_ids:
            row[key] = _format_value(rec.info.get(key))

        for sample_name in sample_names:
            sample = rec.samples[sample_name]
            prefix = f"{sample_name}_" if len(sample_names) > 1 else ""
            for key in format_ids:
                col = f"{prefix}{key}"
                row.setdefault(col, None)
                val = sample.get(key)
                row[col] = _format_value(val, is_gt=(key == "GT"), phased=sample.phased)

        end = rec.info.get("END", rec.stop)
        row["_chrom"] = rec.chrom
        row["_start"] = rec.pos - 1
        row["_end"] = int(end) if end is not None else rec.pos
        rows.append(row)

    df = pd.DataFrame(rows)
    return df, header_meta


def bed_interest_filter(df, region_df):
    """Filtra `df` (de vcf2dataframe ou bed2dataframe) mantendo so as linhas
    que se sobrepoem (parcial ou totalmente) a alguma regiao de `region_df`
    (tipicamente gerado por bed2dataframe a partir de --region).
    """
    keep = pd.Series(False, index=df.index)
    for chrom, region_group in region_df.groupby("_chrom"):
        in_chrom = df["_chrom"] == chrom
        if not in_chrom.any():
            continue
        starts = df.loc[in_chrom, "_start"].to_numpy()
        ends = df.loc[in_chrom, "_end"].to_numpy()
        r_starts = region_group["_start"].to_numpy()
        r_ends = region_group["_end"].to_numpy()
        # overlap[i, j] = variante i sobrepoe regiao j
        overlap = (starts[:, None] < r_ends[None, :]) & (r_starts[None, :] < ends[:, None])
        any_overlap = overlap.any(axis=1)
        keep.loc[df.index[in_chrom]] = any_overlap
    return df.loc[keep].reset_index(drop=True)


# --------------------------------------------------------------------------
# Config (report.config)
# --------------------------------------------------------------------------

CONFIG_FILENAME = "report.config"
PATH_CONFIG_KEYS = ["snv_vcf", "sv_vcf", "cnv_vcf", "str_vcf", "dmr_bed", "region", "bam"]
REQUIRED_CONFIG_KEYS = ["sample"] + PATH_CONFIG_KEYS
CONFIG_DEFAULTS = {"host": "127.0.0.1", "port": "8000"}


def find_config_path():
    """Sem parametros nomeados: usa o argumento posicional se dado, senao
    procura o arquivo `report.config` ao lado deste script."""
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).resolve()
    return Path(__file__).resolve().parent / CONFIG_FILENAME


def load_config(config_path):
    if not config_path.is_file():
        raise SystemExit(
            f"Arquivo de config nao encontrado: {config_path}\n"
            f"Crie um `{CONFIG_FILENAME}` ao lado de {Path(__file__).name} (ou passe o caminho "
            f"como argumento posicional) com as chaves: {', '.join(REQUIRED_CONFIG_KEYS)}."
        )

    values = dict(CONFIG_DEFAULTS)
    for raw_line in config_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        values[key.strip()] = value.strip()

    missing = [k for k in REQUIRED_CONFIG_KEYS if not values.get(k)]
    if missing:
        raise SystemExit(f"Faltando no config ({config_path}): {', '.join(missing)}")

    for key in PATH_CONFIG_KEYS:
        p = Path(values[key])
        if not p.is_absolute():
            p = (config_path.parent / p).resolve()
        values[key] = str(p)

    values["port"] = int(values["port"])
    return SimpleNamespace(**values)


def parse_select_columns(spec):
    """'CHROM;POS;REF;ALT' -> ['CHROM', 'POS', 'REF', 'ALT']"""
    if not spec:
        return []
    return [c.strip() for c in spec.split(";") if c.strip()]


_COLOR_COLUMN_RE = re.compile(r"^([^\[\]]+)\[(.*)\]$")


def parse_color_columns(spec):
    """'ACMG_PRED[Pathogenic:#FF000080|VUS:#D3D3D3];GENE[...]' ->
    {'ACMG_PRED': {'Pathogenic': '#FF000080', 'VUS': '#D3D3D3'}, 'GENE': {...}}"""
    result = {}
    if not spec:
        return result
    for col_spec in spec.split(";"):
        col_spec = col_spec.strip()
        if not col_spec:
            continue
        m = _COLOR_COLUMN_RE.match(col_spec)
        if not m:
            continue
        col, body = m.group(1).strip(), m.group(2)
        value_colors = {}
        for pair in body.split("|"):
            pair = pair.strip()
            if not pair or ":" not in pair:
                continue
            value, color = pair.split(":", 1)
            value_colors[value.strip()] = color.strip()
        if value_colors:
            result[col] = value_colors
    return result


DEFAULT_STANDARD_PALETTE = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79", "#5254a3", "#6b6ecf", "#9c9ede",
    "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
]


def parse_color_palettes(spec):
    """'Standard[#1f77b4|#aec7e8|...];Warm[...]' -> {'Standard': ['#1f77b4', '#aec7e8', ...], ...}"""
    result = {}
    if not spec:
        return result
    for pal_spec in spec.split(";"):
        pal_spec = pal_spec.strip()
        if not pal_spec:
            continue
        m = _COLOR_COLUMN_RE.match(pal_spec)
        if not m:
            continue
        name, body = m.group(1).strip(), m.group(2)
        colors = [c.strip() for c in body.split("|") if c.strip()]
        if colors:
            result[name] = colors
    return result


def get_config_value(config_path, key):
    """Le uma unica chave direto do arquivo de config no disco (sem cache),
    para refletir o que o botao 'Save Filter' acabou de gravar na mesma sessao."""
    if not config_path.is_file():
        return ""
    for raw_line in config_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        if k.strip() == key:
            return v.strip()
    return ""


def set_config_value(config_path, key, value):
    """Atualiza (ou adiciona) `key = value` no arquivo de config, preservando o resto."""
    lines = config_path.read_text().splitlines() if config_path.is_file() else []
    new_line = f"{key} = {value}"
    for i, raw_line in enumerate(lines):
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, _ = line.partition("=")
        if k.strip() == key:
            lines[i] = new_line
            break
    else:
        lines.append(new_line)
    config_path.write_text("\n".join(lines) + "\n")


def filters_to_config_value(filters, display_cols):
    """[{'col': idx, 'value': ...}] -> JSON string com nome de coluna em vez de indice
    (portavel entre selecoes de coluna diferentes)."""
    out = []
    for f in filters:
        idx = f.get("col")
        if idx is None or not (0 <= idx < len(display_cols)):
            continue
        out.append({"col": display_cols[idx], "value": f.get("value")})
    return json.dumps(out)


def config_value_to_filters(value, display_cols):
    """Inverso de filters_to_config_value: nome de coluna -> indice na exibicao atual."""
    if not value:
        return []
    try:
        raw = json.loads(value)
    except (json.JSONDecodeError, TypeError):
        return []
    out = []
    for f in raw:
        name = f.get("col")
        if name not in display_cols:
            continue
        out.append({"col": display_cols.index(name), "value": f.get("value")})
    return out


# --------------------------------------------------------------------------
# UI / colunas
# --------------------------------------------------------------------------

STANDARD_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
INTERNAL_COLS = ("_chrom", "_start", "_end")
BASE_SELECTED = ["CHROM", "POS", "REF", "ALT"]
MAX_DEFAULT_SELECTED_BARS = 30

SECTIONS = ["snv", "sv", "cnv", "str", "dmr"]
SECTION_LABELS = {"snv": "SNV", "sv": "SV", "cnv": "CNV", "str": "STR", "dmr": "DMR"}


def vcf_choices_and_selected(df, header_meta, configured_selected):
    info_ids = [c for c in header_meta["INFO"] if c in df.columns]
    standard = [c for c in STANDARD_COLS if c in df.columns]
    fmt_cols = [
        c for c in df.columns
        if c not in standard and c not in info_ids and c not in INTERNAL_COLS
    ]
    choices = {
        "Standard": {c: c for c in standard},
        "INFO": {c: c for c in info_ids},
        "FORMAT": {c: c for c in fmt_cols},
    }
    configured_selected = [c for c in configured_selected if c in df.columns]
    selected = configured_selected or [c for c in BASE_SELECTED if c in df.columns]
    return choices, selected


def bed_choices_and_selected(df, configured_selected):
    cols = [c for c in df.columns if c not in INTERNAL_COLS]
    configured_selected = [c for c in configured_selected if c in cols]
    selected = configured_selected or cols
    return {"Standard": {c: c for c in cols}}, selected


# --------------------------------------------------------------------------
# Modulo reutilizavel: uma sessao (SNV/SV/CNV/STR/DMR)
# --------------------------------------------------------------------------

def _graphics_field(label_text, widget, extra_class=""):
    return ui.tags.div(
        ui.tags.span(label_text, class_="graphics-label"),
        widget,
        class_=f"graphics-row {extra_class}".strip(),
    )


def _column_bars_fields(column_id, bars_id, choices_by_group):
    return [
        _graphics_field(
            "Select a column",
            ui.input_select(column_id, None, choices=choices_by_group),
        ),
        _graphics_field(
            "Select bars",
            ui.input_selectize(bars_id, None, choices=[], multiple=True),
            extra_class="graphics-row-grow",
        ),
    ]


def _snv_graphics_ui(choices_by_group, palette_names, default_palette):
    return ui.tags.div(
        _graphics_field(
            "Select a chart type",
            ui.input_select(
                "chart_type", None,
                choices=["Custom Bar Chart", "Custom Pizza Chart"],
                selected="Custom Bar Chart",
            ),
        ),
        ui.panel_conditional(
            "input.chart_type == 'Custom Bar Chart'",
            *_column_bars_fields("chart_selected_column", "chart_selected_bars", choices_by_group),
            _graphics_field(
                "Color palette",
                ui.input_select("chart_color_palette", None, choices=palette_names, selected=default_palette),
            ),
            _graphics_field(
                "Bar orientation",
                ui.input_radio_buttons(
                    "bar_orientation", None,
                    choices=["Horizontal", "Vertical"], selected="Vertical", inline=True,
                ),
            ),
            ui.output_plot("bar_chart", height="500px"),
            ui.tags.div(
                ui.input_action_button("save_svg_btn", "Save as SVG"),
                ui.input_action_button("save_png_btn", "Save as PNG"),
                class_="filter-btn-row",
            ),
        ),
        ui.panel_conditional(
            "input.chart_type == 'Custom Pizza Chart'",
            *_column_bars_fields("pizza_selected_column", "pizza_selected_bars", choices_by_group),
            _graphics_field(
                "Color palette",
                ui.input_select("pizza_color_palette", None, choices=palette_names, selected=default_palette),
            ),
            ui.output_plot("pizza_chart", height="500px"),
            ui.tags.div(
                ui.input_action_button("save_pizza_svg_btn", "Save as SVG"),
                ui.input_action_button("save_pizza_png_btn", "Save as PNG"),
                class_="filter-btn-row",
            ),
        ),
        class_="section-table-col",
    )


@module.ui
def variant_section_ui(choices_by_group, default_selected, tooltip_meta, hm_key, palette_names, default_palette):
    graphics_panel = (
        _snv_graphics_ui(choices_by_group, palette_names, default_palette)
        if hm_key == "snv" else ui.p("Em breve.")
    )
    return ui.navset_pill(
        ui.nav_panel(
            "Table",
            ui.tags.div(
                ui.input_switch("bedfilter", "Filter by regions of interest", value=True),
                ui.input_selectize(
                    "select_columns", "Columns",
                    choices=choices_by_group, selected=default_selected, multiple=True,
                    width="100%",
                ),
                ui.tags.div(
                    ui.input_action_button("load_filter_btn", "Load Filter"),
                    ui.input_action_button("reset_filter_btn", "Reset Filter"),
                    ui.input_action_button("save_filter_btn", "Save Filter"),
                    ui.tags.span(class_="filter-btn-spacer"),
                    ui.input_action_button("save_columns_btn", "Save Columns"),
                    ui.input_action_button("load_columns_btn", "Load Columns"),
                    class_="filter-btn-row",
                ),
                ui.output_data_frame("table"),
                ui.tags.div(
                    ui.input_action_button("save_tsv_btn", "Save Tsv"),
                    ui.input_action_button("save_xlsx_btn", "Save Xlsx"),
                    class_="filter-btn-row",
                ),
                ui.tags.script(f"window.__headerMeta = window.__headerMeta || {{}}; "
                                f"window.__headerMeta[{json.dumps(hm_key)}] = {json.dumps(tooltip_meta)};"),
                class_="section-table-col",
                **{"data-hm-key": hm_key},
            ),
        ),
        ui.nav_panel("Statistics", ui.p("Em breve.")),
        ui.nav_panel("Graphics", graphics_panel),
    )


@module.server
def variant_section_server(input: Inputs, output: Outputs, session: Session, *, df, df_interest, header_meta, color_columns, section_key, config_path, sample, export_dir, color_palettes):
    @reactive.calc
    def active_df():
        return df_interest if input.bedfilter() else df

    def display_columns(d):
        cols = list(input.select_columns() or [])
        cols = [c for c in cols if c in d.columns]
        if not cols:
            cols = [c for c in BASE_SELECTED if c in d.columns] or list(d.columns[:4])
        return cols

    @render.data_frame
    def table():
        d = active_df()
        cols = display_columns(d)
        d = d[cols]

        styles = []
        for col, value_colors in color_columns.items():
            if col not in d.columns:
                continue
            for value, color in value_colors.items():
                mask = (d[col] == value).tolist()
                if any(mask):
                    styles.append({"location": "body", "rows": mask, "cols": [col], "style": {"background-color": color}})

        return render.DataGrid(d, filters=True, selection_mode="row", height="420px", width="100%", styles=styles)

    @reactive.calc
    def filtered_active_df():
        """active_df() com o filtro (e ordenacao) da tabela aplicado, mantendo todas
        as colunas. Cai de volta para active_df() quando nenhum filtro esta ativo."""
        return active_df().loc[table.data_view().index]

    @reactive.effect
    async def _navigate_on_click():
        sel = table.data_view(selected=True)
        if sel is None or sel.empty:
            return
        idx = sel.index[0]
        row = active_df().loc[idx]
        chrom, start, end = row["_chrom"], int(row["_start"]), int(row["_end"])
        pad = max(50, int((end - start) * 0.1))
        locus = f"{chrom}:{max(1, start - pad + 1)}-{end + pad}"
        await session.send_custom_message("igv_goto", {"locus": locus})

    config_key = f"{section_key}_column_filter"

    @reactive.effect
    @reactive.event(input.save_filter_btn)
    def _save_filter():
        cols = display_columns(active_df())
        filters = filters_to_config_value(list(table.filter()), cols)
        set_config_value(config_path, config_key, filters)

    @reactive.effect
    @reactive.event(input.load_filter_btn)
    async def _load_filter():
        cols = display_columns(active_df())
        raw = get_config_value(config_path, config_key)
        filters = config_value_to_filters(raw, cols)
        await table.update_filter(filters if filters else None)

    @reactive.effect
    @reactive.event(input.reset_filter_btn)
    async def _reset_filter():
        await table.update_filter(None)

    columns_config_key = f"{section_key}_select_columns"

    @reactive.effect
    @reactive.event(input.save_columns_btn)
    def _save_columns():
        cols = list(input.select_columns() or [])
        set_config_value(config_path, columns_config_key, ";".join(cols))

    @reactive.effect
    @reactive.event(input.load_columns_btn)
    def _load_columns():
        raw = get_config_value(config_path, columns_config_key)
        cols = parse_select_columns(raw)
        valid_cols = [c for c in cols if c in active_df().columns]
        ui.update_selectize("select_columns", selected=valid_cols)

    export_label = SECTION_LABELS[section_key]

    @reactive.effect
    @reactive.event(input.save_tsv_btn)
    def _save_tsv():
        path = export_dir / f"{sample}_{export_label}_filtered.tsv"
        table.data_view().to_csv(path, sep="\t", index=False)
        ui.notification_show(f"Salvo: {path.name}", duration=4)

    @reactive.effect
    @reactive.event(input.save_xlsx_btn)
    def _save_xlsx():
        path = export_dir / f"{sample}_{export_label}_filtered.xlsx"
        table.data_view().to_excel(path, index=False)
        ui.notification_show(f"Salvo: {path.name}", duration=4)

    if section_key == "snv":
        def _update_bars_choices(col, bars_id):
            if not col:
                return
            values = sorted(filtered_active_df()[col].astype(str).unique().tolist())
            default_selected = values[:MAX_DEFAULT_SELECTED_BARS]
            ui.update_selectize(bars_id, choices=values, selected=default_selected)

        @reactive.effect
        def _update_chart_bars_choices():
            _update_bars_choices(input.chart_selected_column(), "chart_selected_bars")

        @reactive.effect
        def _update_pizza_bars_choices():
            _update_bars_choices(input.pizza_selected_column(), "pizza_selected_bars")

        def _counts_for(col, bars_selected):
            counts = filtered_active_df()[col].astype(str).value_counts().sort_values(ascending=False)
            selected_bars = set(bars_selected or [])
            return counts[counts.index.isin(selected_bars)]

        def _palette_colors(palette_name, n):
            colors = color_palettes.get(palette_name) or DEFAULT_STANDARD_PALETTE
            return [colors[i % len(colors)] for i in range(n)]

        def _bar_figure():
            col = input.chart_selected_column()
            orientation = input.bar_orientation()
            counts = _counts_for(col, input.chart_selected_bars())
            colors = _palette_colors(input.chart_color_palette(), len(counts))

            fig, ax = plt.subplots(
                figsize=(8, max(4, 0.35 * len(counts))) if orientation == "Horizontal"
                else (max(8, 0.35 * len(counts)), 5)
            )
            if orientation == "Horizontal":
                ax.barh(counts.index, counts.values, color=colors)
                ax.invert_yaxis()
                ax.set_xlabel("Count")
                ax.set_ylabel(col)
            else:
                ax.bar(counts.index, counts.values, color=colors)
                ax.set_ylabel("Count")
                ax.set_xlabel(col)
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
            ax.set_title(f"{col} distribution")
            fig.tight_layout()
            return fig

        @render.plot
        def bar_chart():
            return _bar_figure()

        @reactive.effect
        @reactive.event(input.save_svg_btn)
        def _save_chart_svg():
            fig = _bar_figure()
            path = export_dir / f"{sample}_{export_label}_{input.chart_selected_column()}_barchart.svg"
            fig.savefig(path, format="svg", bbox_inches="tight")
            plt.close(fig)
            ui.notification_show(f"Salvo: {path.name}", duration=4)

        @reactive.effect
        @reactive.event(input.save_png_btn)
        def _save_chart_png():
            fig = _bar_figure()
            path = export_dir / f"{sample}_{export_label}_{input.chart_selected_column()}_barchart.png"
            fig.savefig(path, format="png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            ui.notification_show(f"Salvo: {path.name}", duration=4)

        def _pie_figure():
            col = input.pizza_selected_column()
            counts = _counts_for(col, input.pizza_selected_bars())

            fig, ax = plt.subplots(figsize=(8, 8))
            if counts.empty:
                ax.text(0.5, 0.5, "No bars selected", ha="center", va="center")
                ax.axis("off")
            else:
                colors = _palette_colors(input.pizza_color_palette(), len(counts))
                ax.pie(counts.values, labels=counts.index, autopct="%1.1f%%", colors=colors, startangle=90)
            ax.set_title(f"{col} distribution")
            fig.tight_layout()
            return fig

        @render.plot
        def pizza_chart():
            return _pie_figure()

        @reactive.effect
        @reactive.event(input.save_pizza_svg_btn)
        def _save_pizza_svg():
            fig = _pie_figure()
            path = export_dir / f"{sample}_{export_label}_{input.pizza_selected_column()}_piechart.svg"
            fig.savefig(path, format="svg", bbox_inches="tight")
            plt.close(fig)
            ui.notification_show(f"Salvo: {path.name}", duration=4)

        @reactive.effect
        @reactive.event(input.save_pizza_png_btn)
        def _save_pizza_png():
            fig = _pie_figure()
            path = export_dir / f"{sample}_{export_label}_{input.pizza_selected_column()}_piechart.png"
            fig.savefig(path, format="png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            ui.notification_show(f"Salvo: {path.name}", duration=4)


# --------------------------------------------------------------------------
# IGV card
# --------------------------------------------------------------------------

def mount_static(paths):
    """Mapeia cada diretorio-pai unico entre `paths` para um alias /dataN,
    retorna (static_assets_dict, url_for(path))."""
    dirs = {}
    static_assets = {}
    for p in paths:
        d = str(Path(p).resolve().parent)
        if d not in dirs:
            alias = f"data{len(dirs)}"
            dirs[d] = alias
            static_assets[f"/{alias}"] = d

    def url_for(p):
        d = str(Path(p).resolve().parent)
        return f"/{dirs[d]}/{Path(p).name}"

    return static_assets, url_for


SECTION_LAYOUT_CSS = ui.tags.style("""
#variants_accordion .accordion-body {
  padding-top: 0;
}
.section-table-col {
  display: flex;
  flex-direction: column;
  align-items: stretch;
  width: 100%;
  gap: 8px;
  margin-top: 16px;
}
.section-table-col .shiny-input-container {
  width: 100% !important;
  max-width: 100% !important;
}
.filter-btn-row {
  display: flex;
  gap: 8px;
}
.filter-btn-row .shiny-input-container {
  width: auto !important;
}
.filter-btn-spacer {
  width: 20px;
}
.section-table-col .shiny-data-grid {
  width: 100% !important;
  text-align: left;
}
.section-table-col .shiny-data-grid table {
  width: 100% !important;
}
.graphics-row {
  display: flex;
  align-items: center;
  gap: 8px;
}
.graphics-row .graphics-label {
  font-weight: 500;
  white-space: nowrap;
}
.graphics-row .shiny-input-container {
  width: auto !important;
  max-width: none !important;
  margin-bottom: 0;
}
.graphics-row-grow .shiny-input-container {
  flex: 1 1 auto;
}

.navbar-brand .omics-navbar-logo {
  height: 32px;
  width: 32px;
  object-fit: contain;
  margin-right: 8px;
}
""")


def navbar_title(logo_url):
    return ui.tags.span(
        ui.tags.img(src=logo_url, alt="logo", class_="omics-navbar-logo"),
        "OMICS Report",
    )


HEADER_TOOLTIP_SCRIPT = ui.tags.script("""
(function () {
  if (window.__hmTooltipInstalled) return;
  window.__hmTooltipInstalled = true;
  const tip = document.createElement('div');
  tip.style.cssText = 'position:fixed;z-index:9999;background:#222;color:#fff;' +
    'padding:6px 8px;border-radius:4px;font-size:12px;white-space:pre;' +
    'pointer-events:none;display:none;max-width:340px;line-height:1.4;';
  document.body.appendChild(tip);

  document.addEventListener('mouseover', function (e) {
    const th = e.target.closest('th');
    if (!th) { tip.style.display = 'none'; return; }
    const wrap = th.closest('[data-hm-key]');
    const meta = wrap && window.__headerMeta ? window.__headerMeta[wrap.getAttribute('data-hm-key')] : null;
    const entry = meta ? meta[th.textContent.trim()] : null;
    if (!entry) { tip.style.display = 'none'; return; }
    tip.textContent = Object.entries(entry).map(([k, v]) => k + ': ' + v).join('\\n');
    tip.style.display = 'block';
  });
  document.addEventListener('mousemove', function (e) {
    if (tip.style.display === 'block') {
      tip.style.left = (e.clientX + 12) + 'px';
      tip.style.top = (e.clientY + 12) + 'px';
    }
  });
})();
""")


def igv_card_ui(track_specs):
    tracks_json = json.dumps(track_specs)
    return ui.card(
        ui.card_header("IGV Browser"),
        ui.tags.div(id="igvDiv"),
        ui.tags.script(f"""
        window.__igvTracks = {tracks_json};
        function __initIgv() {{
          if (typeof igv === 'undefined') {{ setTimeout(__initIgv, 100); return; }}
          igv.createBrowser(document.getElementById('igvDiv'), {{
            genome: 'hg38',
            tracks: window.__igvTracks
          }}).then(function(browser) {{
            window.igvBrowser = browser;
            Shiny.addCustomMessageHandler('igv_goto', function(msg) {{
              window.igvBrowser.search(msg.locus);
            }});
          }});
        }}
        __initIgv();
        """),
    )


IGV_JS_SRC = "https://cdn.jsdelivr.net/npm/igv@3/dist/igv.min.js"

APP_THEME = ui.Theme(preset="shiny").add_defaults(primary="#874D9C")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def build_app(args, config_path):
    region_df = bed2dataframe(args.region)

    vcf_paths = {"snv": args.snv_vcf, "sv": args.sv_vcf, "cnv": args.cnv_vcf, "str": args.str_vcf}
    dfs, dfs_interest, metas = {}, {}, {}
    for key, path in vcf_paths.items():
        df, meta = vcf2dataframe(path)
        dfs[key] = df
        metas[key] = meta
        dfs_interest[key] = bed_interest_filter(df, region_df)

    dmr_df = bed2dataframe(args.dmr_bed)
    dfs["dmr"] = dmr_df
    metas["dmr"] = {"INFO": {}, "FORMAT": {}}
    dfs_interest["dmr"] = bed_interest_filter(dmr_df, region_df)

    logo_path = Path(__file__).resolve().parent / "logo.png"
    static_assets, url_for = mount_static(
        [args.snv_vcf, args.sv_vcf, args.cnv_vcf, args.str_vcf, args.bam, logo_path]
    )

    track_specs = [
        {
            "name": "haplotagged BAM", "type": "alignment", "format": "bam",
            "url": url_for(args.bam), "indexURL": url_for(args.bam) + ".bai",
            "colorBy": "tag", "tag": "HP",
            "visibilityWindow": 5000, "samplingWindowSize": 100, "samplingDepth": 300,
        },
    ]
    for key, label in [("snv", "SNV"), ("sv", "SV"), ("cnv", "CNV"), ("str", "STR")]:
        path = vcf_paths[key]
        track_specs.append({
            "name": label, "type": "variant", "format": "vcf",
            "url": url_for(path), "indexURL": url_for(path) + ".tbi",
        })

    configured_selected = {
        key: parse_select_columns(getattr(args, f"{key}_select_columns", ""))
        for key in SECTIONS
    }
    color_columns = {
        key: parse_color_columns(getattr(args, f"{key}_color_columns", ""))
        for key in SECTIONS
    }
    color_palettes = parse_color_palettes(getattr(args, "color_palettes", "")) or {"Standard": DEFAULT_STANDARD_PALETTE}
    palette_names = list(color_palettes.keys())
    default_palette = "Standard" if "Standard" in color_palettes else palette_names[0]

    accordion_panels = []
    for key in SECTIONS:
        if key == "dmr":
            choices, selected = bed_choices_and_selected(dfs["dmr"], configured_selected[key])
        else:
            choices, selected = vcf_choices_and_selected(dfs[key], metas[key], configured_selected[key])
        tooltip_meta = {**metas[key]["INFO"], **metas[key]["FORMAT"]}
        accordion_panels.append(
            ui.accordion_panel(
                SECTION_LABELS[key],
                variant_section_ui(key, choices, selected, tooltip_meta, key, palette_names, default_palette),
            )
        )

    app_ui = ui.page_navbar(
        ui.nav_panel(
            "Variants",
            ui.card(
                ui.card_header("Variants"),
                ui.accordion(*accordion_panels, id="variants_accordion", open=["SNV"]),
            ),
            igv_card_ui(track_specs),
        ),
        ui.nav_panel("Run", ui.p("Em breve.")),
        ui.nav_panel("Help", ui.p("Em breve.")),
        ui.nav_panel("About", ui.p("Em breve.")),
        ui.nav_spacer(),
        ui.nav_control(ui.input_dark_mode()),
        title=navbar_title(url_for(logo_path)),
        header=ui.TagList(ui.tags.script(src=IGV_JS_SRC), SECTION_LAYOUT_CSS, HEADER_TOOLTIP_SCRIPT),
        id="main_navbar",
        theme=APP_THEME,
    )

    export_dir = Path(__file__).resolve().parent

    def server(input: Inputs, output: Outputs, session: Session):
        for key in SECTIONS:
            variant_section_server(
                key, df=dfs[key], df_interest=dfs_interest[key], header_meta=metas[key],
                color_columns=color_columns[key], section_key=key, config_path=config_path,
                sample=args.sample, export_dir=export_dir, color_palettes=color_palettes,
            )

    return App(app_ui, server, static_assets=static_assets)


if __name__ == "__main__":
    config_path = find_config_path()
    args = load_config(config_path)
    app = build_app(args, config_path)
    app.run(host=args.host, port=args.port)
